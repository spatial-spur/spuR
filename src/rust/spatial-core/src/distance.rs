//! Geographic and Euclidean distances in the host's requested units.

use crate::{checked_len, CoreError};
use rayon::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Metric {
    Euclidean,
    /// Great-circle distance as fraction of 2 * pi * R (equivalently asin(...) / pi).
    /// Input columns are `(lat, lon)` in degrees.
    GreatCircleFractionLatLon,
    /// Same as `GreatCircleFractionLatLon`, with input columns `(lon, lat)` in degrees.
    GreatCircleFractionLonLat,
    /// Great-circle distance in the units of `radius`. Input columns are `(lat, lon)` in degrees.
    GreatCircle {
        radius: f64,
    },
}

#[derive(Clone, Copy)]
struct GeographicPoint {
    lat: f64,
    lon: f64,
    cos_lat: f64,
}

impl GeographicPoint {
    fn new(lat: f64, lon: f64) -> Self {
        let lat = lat.to_radians();
        Self {
            lat,
            lon: lon.to_radians(),
            cos_lat: lat.cos(),
        }
    }

    fn half_angle(self, other: Self) -> f64 {
        let lat = ((self.lat - other.lat) * 0.5).sin();
        let lon = ((self.lon - other.lon) * 0.5).sin();
        (lat * lat + self.cos_lat * other.cos_lat * lon * lon)
            .clamp(0.0, 1.0)
            .sqrt()
            .asin()
    }
}

pub fn haversine(lat1: f64, lon1: f64, lat2: f64, lon2: f64, radius: f64) -> f64 {
    2.0 * radius * GeographicPoint::new(lat1, lon1).half_angle(GeographicPoint::new(lat2, lon2))
}

pub(crate) struct Distances<'a> {
    coordinates: &'a [f64],
    width: usize,
    geographic: Option<Vec<GeographicPoint>>,
    multiplier: f64,
}

impl<'a> Distances<'a> {
    pub(crate) fn new(
        coordinates: &'a [f64],
        n: usize,
        width: usize,
        metric: Metric,
    ) -> Result<Self, CoreError> {
        if width == 0 || coordinates.len() != checked_len(&[n, width])? {
            return Err(CoreError::invalid_input(
                "coordinate dimensions do not match",
            ));
        }
        if coordinates.iter().any(|v| !v.is_finite()) {
            return Err(CoreError::invalid_input("coordinates must be finite"));
        }
        let geographic = match metric {
            Metric::Euclidean => None,
            Metric::GreatCircleFractionLatLon | Metric::GreatCircle { .. } => {
                if width != 2 {
                    return Err(CoreError::invalid_input(
                        "geographic coordinates require two columns",
                    ));
                }
                Some(
                    coordinates
                        .chunks_exact(2)
                        .map(|row| GeographicPoint::new(row[0], row[1]))
                        .collect(),
                )
            }
            Metric::GreatCircleFractionLonLat => {
                if width != 2 {
                    return Err(CoreError::invalid_input(
                        "geographic coordinates require two columns",
                    ));
                }
                Some(
                    coordinates
                        .chunks_exact(2)
                        .map(|row| GeographicPoint::new(row[1], row[0]))
                        .collect(),
                )
            }
        };
        let multiplier = match metric {
            Metric::Euclidean => 1.0,
            Metric::GreatCircleFractionLatLon | Metric::GreatCircleFractionLonLat => {
                std::f64::consts::FRAC_1_PI
            }
            Metric::GreatCircle { radius } if radius.is_finite() && radius > 0.0 => 2.0 * radius,
            Metric::GreatCircle { .. } => {
                return Err(CoreError::invalid_input(
                    "Earth radius must be finite and positive",
                ))
            }
        };
        Ok(Self {
            coordinates,
            width,
            geographic,
            multiplier,
        })
    }

    pub(crate) fn pair(&self, i: usize, j: usize) -> f64 {
        if let Some(points) = &self.geographic {
            self.multiplier * points[i].half_angle(points[j])
        } else {
            let a = &self.coordinates[i * self.width..(i + 1) * self.width];
            let b = &self.coordinates[j * self.width..(j + 1) * self.width];
            let squares = a
                .iter()
                .zip(b)
                .map(|(left, right)| (left - right).powi(2))
                .sum::<f64>();
            if squares.is_finite() {
                squares.sqrt()
            } else {
                a.iter()
                    .zip(b)
                    .fold(0.0_f64, |norm, (left, right)| norm.hypot(left - right))
            }
        }
    }
}

/// The full symmetric distance matrix in row-major order.
pub fn distance_matrix(
    coordinates: &[f64],
    n: usize,
    width: usize,
    metric: Metric,
) -> Result<Vec<f64>, CoreError> {
    let distances = Distances::new(coordinates, n, width, metric)?;
    let mut matrix = vec![0.0; checked_len(&[n, n])?];
    if n == 0 {
        return Ok(matrix);
    }
    if crate::threads::is_serial() || n < 64 {
        for (i, row) in matrix.chunks_exact_mut(n).enumerate() {
            for (j, value) in row[..i].iter_mut().enumerate() {
                *value = distances.pair(i, j);
            }
        }
    } else {
        matrix.par_chunks_mut(n).enumerate().for_each(|(i, row)| {
            for (j, value) in row[..i].iter_mut().enumerate() {
                *value = distances.pair(i, j);
            }
        });
    }
    for i in 0..n {
        for j in 0..i {
            matrix[j * n + i] = matrix[i * n + j];
        }
    }
    if matrix.iter().any(|d| !d.is_finite()) {
        return Err(CoreError::invalid_input(
            "pairwise distance exceeds the floating-point range",
        ));
    }
    Ok(matrix)
}

/// Extract strict lower triangle from a symmetric `n x n` matrix.
pub(crate) fn lower_triangle(distances: &[f64], n: usize) -> Vec<f64> {
    let pairs = n.saturating_mul(n.saturating_sub(1)) / 2;
    let mut out = Vec::with_capacity(pairs);
    for (i, row) in distances.chunks_exact(n).enumerate() {
        out.extend_from_slice(&row[..i]);
    }
    out
}
