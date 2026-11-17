use crate::{checked_len, linalg, CoreError};
use rayon::prelude::*;

const TILE: usize = 256;

/// Project a symmetric covariance kernel for several bases and decay scales.
pub(crate) fn project(
    n: usize,
    bases: &[&[f64]],
    widths: &[usize],
    scales: &[f64],
    fill_tile: impl Fn(usize, usize, usize, usize, &mut [f64]) + Sync,
) -> Result<Vec<Vec<f64>>, CoreError> {
    if n == 0 || bases.is_empty() || bases.len() != widths.len() || scales.is_empty() {
        return Err(CoreError::invalid_input(
            "expected positive covariance-grid dimensions",
        ));
    }
    if scales
        .iter()
        .any(|&scale| !scale.is_finite() || scale < 0.0)
    {
        return Err(CoreError::invalid_input(
            "covariance scales must be finite and nonnegative",
        ));
    }
    let mut offsets = Vec::with_capacity(widths.len());
    let mut scale_width = 0usize;
    for (&basis, &width) in bases.iter().zip(widths) {
        if width == 0 || basis.len() != checked_len(&[n, width])? {
            return Err(CoreError::invalid_input(
                "covariance-grid basis dimensions do not match",
            ));
        }
        if basis.iter().any(|value| !value.is_finite()) {
            return Err(CoreError::invalid_input(
                "covariance-grid basis must be finite",
            ));
        }
        offsets.push(scale_width);
        scale_width = scale_width
            .checked_add(checked_len(&[width, width])?)
            .ok_or(CoreError::invalid_input("covariance grid is too large"))?;
    }
    let size = checked_len(&[scales.len(), scale_width])?;
    let blocks = n.div_ceil(TILE);
    let mut tiles = Vec::new();
    for i in 0..blocks {
        for j in i..blocks {
            tiles.push((i * TILE, j * TILE));
        }
    }
    let grid = Projection {
        n,
        bases,
        widths,
        offsets: &offsets,
        scale_width,
        fill_tile: &fill_tile,
    };
    let mut result = vec![0.0; size];
    if crate::threads::is_serial() {
        grid.evaluate(&tiles, scales, &mut result);
    } else {
        let workers = rayon::current_num_threads().max(1);
        if tiles.len() >= workers || scales.len() == 1 {
            let tiles_per_shard = tiles.len().div_ceil(workers);
            let partials: Vec<_> = tiles
                .par_chunks(tiles_per_shard)
                .map(|chunk| {
                    let mut partial = vec![0.0; size];
                    grid.evaluate(chunk, scales, &mut partial);
                    partial
                })
                .collect();
            for partial in partials {
                for (total, value) in result.iter_mut().zip(partial) {
                    *total += value;
                }
            }
        } else {
            result
                .par_chunks_mut(scale_width)
                .zip(scales.par_iter())
                .for_each(|(out, scale)| {
                    grid.evaluate(&tiles, std::slice::from_ref(scale), out);
                });
        }
    }

    Ok(widths
        .iter()
        .zip(offsets)
        .map(|(&width, start)| {
            let matrix_size = width * width;
            result
                .chunks_exact(scale_width)
                .flat_map(|scale| scale[start..start + matrix_size].iter().copied())
                .collect()
        })
        .collect())
}

struct Projection<'a, F> {
    n: usize,
    bases: &'a [&'a [f64]],
    widths: &'a [usize],
    offsets: &'a [usize],
    scale_width: usize,
    fill_tile: &'a F,
}

impl<F: Fn(usize, usize, usize, usize, &mut [f64])> Projection<'_, F> {
    fn evaluate(&self, tiles: &[(usize, usize)], scales: &[f64], output: &mut [f64]) {
        let width_max = *self.widths.iter().max().expect("at least one basis");
        let block_size = self.n.min(TILE);
        let mut distance = vec![0.0; block_size * block_size];
        let mut kernel = vec![0.0; block_size * block_size];
        let mut product = vec![0.0; block_size * width_max];
        let mut contribution = vec![0.0; width_max * width_max];

        for &(i, j) in tiles {
            let ni = (self.n - i).min(TILE);
            let nj = (self.n - j).min(TILE);
            let distance = &mut distance[..ni * nj];
            let kernel = &mut kernel[..ni * nj];
            (self.fill_tile)(i, j, ni, nj, distance);
            for (&scale, out) in scales.iter().zip(output.chunks_exact_mut(self.scale_width)) {
                linalg::exp_scaled(distance, -scale, kernel);
                for ((&basis, &q), &offset) in self.bases.iter().zip(self.widths).zip(self.offsets)
                {
                    linalg::multiply(
                        kernel,
                        &basis[j * q..(j + nj) * q],
                        &mut product,
                        ni,
                        nj,
                        q,
                        false,
                    );
                    linalg::multiply(
                        &basis[i * q..(i + ni) * q],
                        &product,
                        &mut contribution,
                        q,
                        ni,
                        q,
                        true,
                    );
                    for a in 0..q {
                        for b in a..q {
                            let ab = a * q + b;
                            let ba = b * q + a;
                            let value = if i == j {
                                0.5 * (contribution[ab] + contribution[ba])
                            } else {
                                contribution[ab] + contribution[ba]
                            };
                            out[offset + ab] += value;
                            if ab != ba {
                                out[offset + ba] += value;
                            }
                        }
                    }
                }
            }
        }
    }
}
