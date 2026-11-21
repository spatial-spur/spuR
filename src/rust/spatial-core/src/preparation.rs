//! Regression-side preparation for exact and Nyström SCPC. All matrices are row-major.

use crate::{
    checked_len, linalg, nystrom::normalize_coords, run_approx_scpc, run_exact_scpc,
    ApproxScpcRequest, CoreError, ExactScpcRequest, Metric, ScpcOutput,
};

const RANK_TOL: f64 = 1e-10;

/// Inputs extracted from the fitted host model. All indices are zero-based.
#[derive(Clone, Copy)]
pub struct ScpcPreparation<'a> {
    pub n_sites: usize,
    pub n_obs: usize,
    pub n_scores: usize,
    /// The first k coefficients being reported.
    pub coefficients: &'a [f64],
    /// First k rows of the fitted model's bread (k × n_scores).
    pub bread: &'a [f64],
    /// Fitted score contribution for each observation (n_obs × n_scores).
    pub scores: &'a [f64],
    /// Actual regressor columns, including the fitted intercept, before FE absorption.
    /// For IV these are structural regressors. None selects unconditional inference.
    pub design: Option<&'a [f64]>,
    pub design_cols: usize,
    /// Full instrument design, including exogenous regressors and the fitted intercept.
    pub instruments: Option<&'a [f64]>,
    pub instrument_cols: usize,
    /// Optional fitted regressor matrix in the score columns' order (n_obs × n_scores).
    /// Otherwise OLS uses the supplied design and IV computes P_Z X in the core.
    pub direction_design: Option<&'a [f64]>,
    pub cluster_ids: Option<&'a [usize]>,
    /// One-way fixed effects; the same groups are applied later to spatial directions.
    pub fixed_effects: Option<&'a [usize]>,
}

#[derive(Debug)]
pub struct PreparedScpcInputs {
    pub wj: Vec<f64>,
    pub xj: Option<Vec<f64>>,
    pub xjs: Option<Vec<f64>>,
    /// Left and right projection factors, already transposed as required for IV.
    pub proj_u: Option<Vec<f64>>,
    pub proj_v: Option<Vec<f64>>,
    pub proj_rank: usize,
}

fn check_matrix(a: &[f64], rows: usize, cols: usize) -> Result<(), CoreError> {
    if a.len() != checked_len(&[rows, cols])? || a.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input(
            "SCPC model matrix has invalid dimensions or non-finite entries",
        ));
    }
    Ok(())
}

impl ScpcPreparation<'_> {
    fn validate(&self) -> Result<(), CoreError> {
        let k = self.coefficients.len();
        if self.n_sites < 2
            || self.n_obs < self.n_sites
            || k == 0
            || k > self.n_scores
            || [
                self.n_obs,
                self.n_scores,
                self.design_cols,
                self.instrument_cols,
            ]
            .into_iter()
            .any(|v| v > i32::MAX as usize)
        {
            return Err(CoreError::invalid_input("invalid SCPC model dimensions"));
        }
        check_matrix(self.coefficients, 1, k)?;
        check_matrix(self.bread, k, self.n_scores)?;
        check_matrix(self.scores, self.n_obs, self.n_scores)?;
        if let Some(ids) = self.cluster_ids {
            if ids.len() != self.n_obs || ids.iter().any(|&g| g >= self.n_sites) {
                return Err(CoreError::invalid_input("invalid SCPC cluster indices"));
            }
            let mut seen = vec![false; self.n_sites];
            for &g in ids {
                seen[g] = true;
            }
            if seen.contains(&false) {
                return Err(CoreError::invalid_input(
                    "SCPC cluster indices must cover the spatial sites",
                ));
            }
        } else if self.n_obs != self.n_sites {
            return Err(CoreError::invalid_input(
                "SCPC requires cluster indices when observations exceed sites",
            ));
        }
        if let Some(ids) = self.fixed_effects {
            if ids.len() != self.n_obs || ids.iter().any(|&g| g >= self.n_obs) {
                return Err(CoreError::invalid_input(
                    "invalid SCPC fixed-effect indices",
                ));
            }
        }
        if let Some(a) = self.design {
            if self.design_cols == 0
                || (self.direction_design.is_none() && self.design_cols != self.n_scores)
            {
                return Err(CoreError::invalid_input(
                    "SCPC design columns do not match the fitted scores",
                ));
            }
            check_matrix(a, self.n_obs, self.design_cols)?;
        } else if self.design_cols != 0
            || self.instruments.is_some()
            || self.direction_design.is_some()
        {
            return Err(CoreError::invalid_input(
                "conditional SCPC inputs require a model design",
            ));
        }
        if let Some(z) = self.instruments {
            if self.instrument_cols == 0 {
                return Err(CoreError::invalid_input(
                    "IV residualizer is rank deficient: empty instrument design",
                ));
            }
            check_matrix(z, self.n_obs, self.instrument_cols)?;
        } else if self.instrument_cols != 0 {
            return Err(CoreError::invalid_input(
                "instrument columns require an instrument matrix",
            ));
        }
        if let Some(x) = self.direction_design {
            check_matrix(x, self.n_obs, self.n_scores)?;
        }
        Ok(())
    }
}

/// Forward spacing within a binary64 binade, including the finite endpoint.
fn spacing(magnitude: f64) -> f64 {
    if magnitude == f64::MAX {
        magnitude - f64::from_bits(magnitude.to_bits() - 1)
    } else {
        f64::from_bits(magnitude.to_bits() + 1) - magnitude
    }
}

/// Remove a one-way fixed effect in group-scaled coordinates so translations do
/// not hide representable within variation and sums cannot overflow.
fn within(a: &mut [f64], cols: usize, ids: &[usize], screen: bool) -> Result<(), CoreError> {
    let groups = ids.iter().copied().max().unwrap_or(0) + 1;
    let mut first = vec![usize::MAX; groups];
    let mut count = vec![0usize; groups];
    for (i, &g) in ids.iter().enumerate() {
        if count[g] == 0 {
            first[g] = i;
        }
        count[g] += 1;
    }
    for col in 0..cols {
        let mut largest = vec![0.0_f64; groups];
        let mut scale = vec![0.0_f64; groups];
        for (i, &g) in ids.iter().enumerate() {
            let x = a[i * cols + col];
            largest[g] = largest[g].max(x.abs());
            scale[g] = scale[g].max((x - a[first[g] * cols + col]).abs());
        }
        let mut anchor = vec![0.0; groups];
        for g in 0..groups {
            if scale[g] < largest[g] {
                anchor[g] = a[first[g] * cols + col];
            } else {
                scale[g] = largest[g];
            }
            if scale[g] == 0.0 {
                scale[g] = 1.0;
            }
        }
        let mut mean = vec![0.0; groups];
        let mut correction = vec![0.0; groups];
        for (i, &g) in ids.iter().enumerate() {
            let x = &mut a[i * cols + col];
            *x = (*x - anchor[g]) / scale[g];
            let add = *x - correction[g];
            let sum = mean[g] + add;
            correction[g] = (sum - mean[g]) - add;
            mean[g] = sum;
        }
        for g in 0..groups {
            if count[g] > 0 {
                mean[g] /= count[g] as f64;
            }
        }
        let mut amplitude = vec![0.0_f64; groups];
        for (i, &g) in ids.iter().enumerate() {
            let x = &mut a[i * cols + col];
            *x -= mean[g];
            amplitude[g] = amplitude[g].max(x.abs());
        }
        for (i, &g) in ids.iter().enumerate() {
            let x = &mut a[i * cols + col];
            if screen && 2.0 * amplitude[g] < spacing(largest[g]) / scale[g] {
                *x = 0.0;
            } else {
                *x *= scale[g];
            }
            if !x.is_finite() {
                return Err(CoreError::invalid_input(
                    "fixed-effect demeaning produced a non-finite entry",
                ));
            }
        }
    }
    Ok(())
}

fn equilibrate(a: &mut [f64], cols: usize) -> Vec<f64> {
    let mut scales = vec![0.0_f64; cols];
    for row in a.chunks_exact(cols) {
        for (s, x) in scales.iter_mut().zip(row) {
            *s = s.max(x.abs());
        }
    }
    for row in a.chunks_exact_mut(cols) {
        for (x, &s) in row.iter_mut().zip(&scales) {
            if s > 0.0 {
                *x /= s;
            }
        }
    }
    scales
}

/// Thin SVD A = U diag(s) V'. U and V' are row-major, with V' omitted when unused.
struct ThinSvd {
    u: Vec<f64>,
    s: Vec<f64>,
    vt: Vec<f64>,
}

#[cfg(all(target_vendor = "apple", not(feature = "portable")))]
fn thin_svd(a: &[f64], rows: usize, cols: usize, want_right: bool) -> Result<ThinSvd, CoreError> {
    #[link(name = "Accelerate", kind = "framework")]
    extern "C" {
        fn dgesvd_(
            jobu: *const u8,
            jobvt: *const u8,
            m: *const i32,
            n: *const i32,
            a: *mut f64,
            lda: *const i32,
            s: *mut f64,
            u: *mut f64,
            ldu: *const i32,
            vt: *mut f64,
            ldvt: *const i32,
            work: *mut f64,
            lwork: *const i32,
            info: *mut i32,
        );
    }
    // LAPACK sees the row-major matrix as A'. Its column-major VT and U are
    // already the required row-major U and V' of A, respectively.
    let m = cols as i32;
    let n = rows as i32;
    let rank = rows.min(cols);
    let mut data = a.to_vec();
    let mut s = vec![0.0; rank];
    let mut u = vec![0.0; rows * rank];
    let mut vt = vec![0.0; if want_right { cols * rank } else { 1 }];
    let job = if want_right { b'S' } else { b'N' };
    let ldu = if want_right { m } else { 1 };
    let mut decompose = |work: &mut [f64], lwork: i32| -> Result<(), CoreError> {
        let mut info = 0;
        // All input and output buffers cover the dimensions above. A one-word
        // work buffer is the LAPACK workspace query when lwork is -1.
        unsafe {
            dgesvd_(
                &job,
                &b'S',
                &m,
                &n,
                data.as_mut_ptr(),
                &m,
                s.as_mut_ptr(),
                vt.as_mut_ptr(),
                &ldu,
                u.as_mut_ptr(),
                &(rank as i32),
                work.as_mut_ptr(),
                &lwork,
                &mut info,
            );
        }
        if info != 0 {
            return Err(CoreError::internal("SCPC model SVD failed to converge"));
        }
        Ok(())
    };
    let mut query = [0.0];
    decompose(&mut query, -1)?;
    if !(1.0..=i32::MAX as f64).contains(&query[0]) {
        return Err(CoreError::invalid_input(
            "SCPC model SVD workspace is too large",
        ));
    }
    let lwork = query[0].ceil() as i32;
    decompose(&mut vec![0.0; lwork as usize], lwork)?;
    if !want_right {
        vt.clear();
    }
    Ok(ThinSvd { u, s, vt })
}

#[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
fn thin_svd(a: &[f64], rows: usize, cols: usize, want_right: bool) -> Result<ThinSvd, CoreError> {
    use faer::linalg::svd::{svd, svd_scratch, ComputeSvdVectors};
    let rank = rows.min(cols);
    let mut u = vec![0.0; rows * rank];
    let mut vt = if want_right {
        vec![0.0; cols * rank]
    } else {
        Vec::new()
    };
    let mut diag = faer::diag::Diag::<f64>::zeros(rank);
    let right = if want_right {
        ComputeSvdVectors::Thin
    } else {
        ComputeSvdVectors::No
    };
    let par = linalg::faer_par();
    let req = svd_scratch::<f64>(
        rows,
        cols,
        ComputeSvdVectors::Thin,
        right,
        par,
        faer::prelude::default(),
    );
    let mut mem = faer::dyn_stack::MemBuffer::new(req);
    svd(
        faer::MatRef::from_row_major_slice(a, rows, cols),
        diag.as_mut(),
        Some(faer::MatMut::from_row_major_slice_mut(&mut u, rows, rank)),
        want_right.then(|| faer::MatMut::from_column_major_slice_mut(&mut vt, cols, rank)),
        par,
        faer::dyn_stack::MemStack::new(&mut mem),
        faer::prelude::default(),
    )
    .map_err(|_| CoreError::internal("SCPC model SVD failed to converge"))?;
    let s = (0..rank).map(|i| diag.column_vector()[i]).collect();
    Ok(ThinSvd { u, s, vt })
}

fn product(
    a: &[f64],
    b: &[f64],
    rows: usize,
    inner: usize,
    cols: usize,
    transpose_a: bool,
) -> Result<Vec<f64>, CoreError> {
    let mut result = vec![0.0; checked_len(&[rows, cols])?];
    linalg::multiply(a, b, &mut result, rows, inner, cols, transpose_a);
    Ok(result)
}

fn transpose(a: &[f64], rows: usize, cols: usize) -> Vec<f64> {
    let mut result = vec![0.0; a.len()];
    for i in 0..rows {
        for j in 0..cols {
            result[j * rows + i] = a[i * cols + j];
        }
    }
    result
}

fn directions(
    xj: &[f64],
    k: usize,
    sites: usize,
    cluster: Option<&[usize]>,
) -> Result<Vec<f64>, CoreError> {
    let Some(ids) = cluster else {
        return Ok(xj
            .iter()
            .map(|&v| if v < 0.0 { -1.0 } else { 1.0 })
            .collect());
    };
    let mut scales = vec![0.0_f64; checked_len(&[sites, k])?];
    for (row, &g) in xj.chunks_exact(k).zip(ids) {
        for j in 0..k {
            scales[g * k + j] = scales[g * k + j].max(row[j].abs());
        }
    }
    let mut ss = vec![0.0_f64; scales.len()];
    let mut result = xj.to_vec();
    for (row, &g) in result.chunks_exact_mut(k).zip(ids) {
        for j in 0..k {
            if scales[g * k + j] > 0.0 {
                row[j] /= scales[g * k + j];
                ss[g * k + j] += row[j] * row[j];
            }
        }
    }
    for s in &mut ss {
        *s = s.sqrt();
    }
    for (row, &g) in result.chunks_exact_mut(k).zip(ids) {
        for j in 0..k {
            if ss[g * k + j] > 0.0 {
                row[j] /= ss[g * k + j];
            }
        }
    }
    Ok(result)
}

/// Construct model influence values, conditional directions and projection factors.
pub fn prepare_scpc_inputs(req: ScpcPreparation<'_>) -> Result<PreparedScpcInputs, CoreError> {
    crate::threads::check_native_process()?;
    req.validate()?;
    let k = req.coefficients.len();
    let mut bread_t = transpose(req.bread, k, req.n_scores);
    let grouped;
    let scores = if let Some(ids) = req.cluster_ids {
        let mut sums = vec![0.0; checked_len(&[req.n_sites, req.n_scores])?];
        for (&group, row) in ids.iter().zip(req.scores.chunks_exact(req.n_scores)) {
            for (col, &value) in row.iter().enumerate() {
                sums[group * req.n_scores + col] += value;
            }
        }
        grouped = sums;
        grouped.as_slice()
    } else {
        req.scores
    };
    let mut wj = product(scores, &bread_t, req.n_sites, req.n_scores, k, false)?;
    let n = req.n_sites as f64;
    for row in wj.chunks_exact_mut(k) {
        for (v, coef) in row.iter_mut().zip(req.coefficients) {
            *v = n * *v + coef;
        }
    }
    if wj.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input(
            "SCPC model influence values are non-finite",
        ));
    }
    let mut out = PreparedScpcInputs {
        wj,
        xj: None,
        xjs: None,
        proj_u: None,
        proj_v: None,
        proj_rank: 0,
    };
    let Some(raw) = req.design else {
        return Ok(out);
    };
    let rows = req.n_obs;
    let p = req.design_cols;
    let mut x = raw.to_vec();
    if let Some(ids) = req.fixed_effects {
        within(&mut x, p, ids, true)?;
    }

    let mut xj = if let Some(fitted) = req.direction_design {
        if let Some(ids) = req.fixed_effects {
            let mut fitted = fitted.to_vec();
            within(&mut fitted, req.n_scores, ids, false)?;
            product(&fitted, &bread_t, rows, req.n_scores, k, false)?
        } else {
            product(fitted, &bread_t, rows, req.n_scores, k, false)?
        }
    } else if req.instruments.is_none() {
        product(&x, &bread_t, rows, p, k, false)?
    } else {
        Vec::new()
    };

    let scales = equilibrate(&mut x, p);
    if let Some(raw_z) = req.instruments {
        let deficient = CoreError::invalid_input(
            "IV residualizer is rank deficient: projected design is not full rank",
        );
        if rows < p || scales.contains(&0.0) {
            return Err(deficient);
        }
        let mut z = raw_z.to_vec();
        if let Some(ids) = req.fixed_effects {
            within(&mut z, req.instrument_cols, ids, true)?;
        }
        equilibrate(&mut z, req.instrument_cols);
        let zsvd = thin_svd(&z, rows, req.instrument_cols, false)?;
        let rank = zsvd
            .s
            .iter()
            .take_while(|&&s| s > RANK_TOL * zsvd.s[0])
            .count();
        if rank < p {
            return Err(deficient);
        }
        let qz = if rank == zsvd.s.len() {
            zsvd.u
        } else {
            zsvd.u
                .chunks_exact(zsvd.s.len())
                .flat_map(|row| row[..rank].iter().copied())
                .collect()
        };
        // C = Qz' X. Its SVD avoids squaring the condition number in X' Pz X.
        let c = product(&qz, &x, rank, rows, p, true)?;
        let mut csvd = thin_svd(&c, rank, p, true)?;
        if csvd.s[p - 1] <= RANK_TOL * csvd.s[0].max(1.0) {
            return Err(deficient);
        }
        for row in csvd.u.chunks_exact_mut(p) {
            for (u, s) in row.iter_mut().zip(&csvd.s) {
                *u /= s;
            }
        }
        let c_dual = product(&csvd.u, &csvd.vt, rank, p, p, false)?;
        out.proj_u = Some(product(&qz, &c_dual, rows, rank, p, false)?);
        out.proj_v = Some(x); // Xi_s' M_IV' Xi R requires the transposed IV projection.
        if req.direction_design.is_none() {
            for (row, scale) in bread_t.chunks_exact_mut(k).zip(scales) {
                for b in row {
                    *b *= scale;
                }
            }
            let qz_scores = product(&c, &bread_t, rank, p, k, false)?;
            xj = product(&qz, &qz_scores, rows, rank, k, false)?;
        }
    } else {
        let deficient = CoreError::invalid_input("OLS residualizer is rank deficient");
        if rows < p || scales.contains(&0.0) {
            return Err(deficient);
        }
        let svd = thin_svd(&x, rows, p, false)?;
        if svd.s[p - 1] <= RANK_TOL * svd.s[0] {
            return Err(deficient);
        }
        out.proj_u = Some(svd.u.clone());
        out.proj_v = Some(svd.u);
    }
    for value in &mut xj {
        *value *= n;
    }
    if xj.iter().any(|v| !v.is_finite()) {
        return Err(CoreError::invalid_input(
            "SCPC model score directions are non-finite",
        ));
    }
    out.xjs = Some(directions(&xj, k, req.n_sites, req.cluster_ids)?);
    out.xj = Some(xj);
    out.proj_rank = p;
    Ok(out)
}

/// Prepare and run the host model as one operation. `seed=None` uses exact SCPC.
/// For Nyström, coordinates and every aligned site/observation vector are put
/// in the same canonical order before entering the prepared-input API.
#[allow(clippy::too_many_arguments)]
pub fn run_model_scpc(
    prep: ScpcPreparation<'_>,
    coords: &[f64],
    d: usize,
    metric: Metric,
    avc: f64,
    cvs: bool,
    nthreads: Option<usize>,
    seed: Option<u32>,
) -> Result<ScpcOutput, CoreError> {
    let n = prep.n_sites;
    let k = prep.coefficients.len();
    let uncond = prep.design.is_none();
    if nthreads == Some(0) {
        return Err(CoreError::invalid_input("thread count must be at least 1"));
    }
    if d == 0 || coords.len() != checked_len(&[n, d])? {
        return Err(CoreError::invalid_input(
            "invalid SCPC coordinate dimensions",
        ));
    }
    let mut inputs = prepare_scpc_inputs(prep)?;
    let Some(seed) = seed else {
        return run_exact_scpc(&ExactScpcRequest {
            coords,
            n,
            d,
            metric,
            avc,
            uncond,
            cvs,
            coef: prep.coefficients,
            wj: &inputs.wj,
            cluster_id: prep.cluster_ids,
            fe_id: prep.fixed_effects,
            xj: inputs.xj.as_deref(),
            xjs: inputs.xjs.as_deref(),
            proj_u: inputs.proj_u.as_deref(),
            proj_v: inputs.proj_v.as_deref(),
            proj_rank: inputs.proj_rank,
            cv_tolerance: 0.001,
            nthreads,
        });
    };
    let (coords, perm) = normalize_coords(coords, n, d, metric)?;
    let permute = |a: &[f64], width: usize| -> Vec<f64> {
        perm.iter()
            .flat_map(|&old| a[old * width..(old + 1) * width].iter().copied())
            .collect()
    };
    inputs.wj = permute(&inputs.wj, k);
    let mut fe = prep.fixed_effects.map(std::borrow::Cow::Borrowed);
    let clusters = if let Some(ids) = prep.cluster_ids {
        let mut inverse = vec![0; n];
        for (i, &old) in perm.iter().enumerate() {
            inverse[old] = i;
        }
        Some(ids.iter().map(|&g| inverse[g]).collect::<Vec<_>>())
    } else {
        for (v, width) in [
            (&mut inputs.xj, k),
            (&mut inputs.xjs, k),
            (&mut inputs.proj_u, inputs.proj_rank),
            (&mut inputs.proj_v, inputs.proj_rank),
        ] {
            if let Some(a) = v {
                *a = permute(a, width);
            }
        }
        if let Some(ids) = prep.fixed_effects {
            fe = Some(std::borrow::Cow::Owned(
                perm.iter().map(|&i| ids[i]).collect(),
            ));
        }
        None
    };
    run_approx_scpc(&ApproxScpcRequest {
        coords: &coords,
        n,
        d,
        metric,
        avc,
        uncond,
        cvs,
        coef: prep.coefficients,
        wj: &inputs.wj,
        cluster_id: clusters.as_deref(),
        fe_id: fe.as_deref(),
        xj: inputs.xj.as_deref(),
        xjs: inputs.xjs.as_deref(),
        proj_u: inputs.proj_u.as_deref(),
        proj_v: inputs.proj_v.as_deref(),
        proj_rank: inputs.proj_rank,
        cv_tolerance: 0.001,
        nthreads,
        seed,
        landmarks: 1000,
        samples: 20,
        pair_count: 1_000_000,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn selected_ols_columns_share_one_scale_invariant_rank_cutoff() {
        let x = [-1., 1., -1., 1.];
        let z = [-1., -1., 1., 1.];
        for scale in [1e-200, 1.0, 1e200] {
            for (gap, identified) in [(1e-9, true), (1e-11, false), (1e-13, false)] {
                let design: Vec<f64> = x
                    .iter()
                    .zip(z)
                    .flat_map(|(&x, z)| [x, (x + gap * z) * scale])
                    .collect();
                let req = ScpcPreparation {
                    n_sites: 4,
                    n_obs: 4,
                    n_scores: 2,
                    coefficients: &[0.0],
                    bread: &[1., 0.],
                    scores: &[0.; 8],
                    design: Some(&design),
                    design_cols: 2,
                    instruments: None,
                    instrument_cols: 0,
                    direction_design: None,
                    cluster_ids: None,
                    fixed_effects: None,
                };
                assert_eq!(
                    prepare_scpc_inputs(req).is_ok(),
                    identified,
                    "scale={scale}, gap={gap}"
                );
            }
        }
    }

    #[test]
    fn grouping_preserves_representable_contrasts_at_the_largest_binade() {
        let max = f64::MAX;
        let hi = f64::from_bits(max.to_bits() - 1);
        let lo = f64::from_bits(max.to_bits() - 2);
        let ulp = max - hi;
        let mut values = [max, 0.0, hi, 1.0, lo, 2.0];
        within(&mut values, 1, &[0, 1, 0, 1, 0, 1], true).unwrap();
        assert_eq!(values, [ulp, -1., 0., 0., -ulp, 1.]);
        let anchor = 2.0_f64.powi(52);
        let mut one_step = [anchor, anchor + 1., anchor, anchor + 1.];
        within(&mut one_step, 1, &[0, 0, 1, 1], true).unwrap();
        assert_eq!(one_step, [-0.5, 0.5, -0.5, 0.5]);
        let tiny = f64::from_bits(1);
        assert_eq!(spacing(tiny), tiny);
        let mut subnormal = [0., tiny, 2. * tiny];
        within(&mut subnormal, 1, &[0; 3], true).unwrap();
        assert_eq!(subnormal, [-tiny, 0., tiny]);
    }

    #[test]
    fn clustered_norms_handle_extreme_units_and_zero_has_positive_unclustered_sign() {
        let x = [-3.0, 0.0, 4.0, -0.0];
        assert_eq!(directions(&x, 1, 4, None).unwrap(), [-1.0, 1.0, 1.0, 1.0]);
        for scale in [1e-200, 1.0, 1e200] {
            let values = x.map(|v| v * scale);
            let normalized = directions(&values, 1, 2, Some(&[0, 0, 0, 1])).unwrap();
            for (a, b) in normalized.iter().zip([-0.6, 0., 0.8, 0.]) {
                assert!((a - b).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn iv_uses_the_fitted_score_design_with_redundant_instruments_and_fixed_effects() {
        let fe = [0, 0, 0, 1, 1, 1];
        let x = [3., 5., 7., 6., 8., 10.];
        let z: Vec<_> = x.iter().flat_map(|v| [*v, 0., *v * 2.]).collect();
        let fitted = [11., 13., 9., -8., -2., -5.];
        let req = ScpcPreparation {
            n_sites: 3,
            n_obs: 6,
            n_scores: 1,
            coefficients: &[0.5],
            bread: &[2.],
            scores: &[2., -3., 4., -4., 3., -2.],
            design: Some(&x),
            design_cols: 1,
            instruments: Some(&z),
            instrument_cols: 3,
            direction_design: Some(&fitted),
            cluster_ids: Some(&[0, 0, 1, 1, 2, 2]),
            fixed_effects: Some(&fe),
        };
        let out = prepare_scpc_inputs(req).unwrap();
        assert_eq!(out.wj, [-5.5, 0.5, 6.5]);
        assert_eq!(out.xj.unwrap(), [0., 12., -12., -18., 18., 0.]);
        let u = out.proj_u.unwrap();
        let v = out.proj_v.unwrap();
        assert!((u.iter().zip(&v).map(|(a, b)| a * b).sum::<f64>() - 1.).abs() < 1e-12);
        // Removing the optional fitted matrix makes the core project the
        // structural X through Z; in this example X is already in that span.
        let projected = prepare_scpc_inputs(ScpcPreparation {
            direction_design: None,
            ..req
        })
        .unwrap()
        .xj
        .unwrap();
        for (a, b) in projected.iter().zip([-12., 0., 12., -12., 0., 12.]) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn iv_rejects_roundoff_only_projections_even_when_their_condition_number_is_one() {
        let x = [1., 0., 0., 1., 0., 0., 0., 0.];
        for projected in [1e-16, 1e-310] {
            let z = [projected, 0., 0., projected, 1., 0., 0., 1.];
            let req = ScpcPreparation {
                n_sites: 4,
                n_obs: 4,
                n_scores: 2,
                coefficients: &[0.],
                bread: &[1., 0.],
                scores: &[0.; 8],
                design: Some(&x),
                design_cols: 2,
                instruments: Some(&z),
                instrument_cols: 2,
                direction_design: None,
                cluster_ids: None,
                fixed_effects: None,
            };
            assert!(prepare_scpc_inputs(req)
                .unwrap_err()
                .message()
                .contains("rank deficient"));
        }
    }
}
