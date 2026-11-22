use std::f64::consts::FRAC_2_PI;

use super::*;
use crate::student_t::betainc;
use crate::ErrorKind;

fn solve_linear(a: &mut [f64], b: &mut [f64], p: usize) {
    for col in 0..p {
        let mut pivot = col;
        let mut max_val = a[col * p + col].abs();
        for row in (col + 1)..p {
            if a[row * p + col].abs() > max_val {
                max_val = a[row * p + col].abs();
                pivot = row;
            }
        }
        if pivot != col {
            for k in 0..p {
                a.swap(col * p + k, pivot * p + k);
            }
            b.swap(col, pivot);
        }
        let diag = a[col * p + col];
        for k in col..p {
            a[col * p + k] /= diag;
        }
        b[col] /= diag;
        for row in 0..p {
            if row == col {
                continue;
            }
            let factor = a[row * p + col];
            for k in col..p {
                a[row * p + k] -= factor * a[col * p + k];
            }
            b[row] -= factor * b[col];
        }
    }
}

fn invert_symmetric(mat: &[f64], p: usize) -> Vec<f64> {
    let mut inv = vec![0.0; p * p];
    for col in 0..p {
        let mut a = mat.to_vec();
        let mut rhs = vec![0.0; p];
        rhs[col] = 1.0;
        solve_linear(&mut a, &mut rhs, p);
        for row in 0..p {
            inv[row * p + col] = rhs[row];
        }
    }
    inv
}

fn modified_gram_schmidt(x: &[f64], n: usize, p: usize) -> Vec<f64> {
    let mut q = vec![0.0; n * p];
    let mut rank = 0usize;
    for col in 0..p {
        let mut v: Vec<f64> = (0..n).map(|i| x[i * p + col]).collect();
        let orig_norm = v.iter().map(|val| val * val).sum::<f64>().sqrt();
        if orig_norm == 0.0 {
            continue;
        }
        for val in &mut v {
            *val /= orig_norm;
        }
        for _ in 0..2 {
            for r in 0..rank {
                let dot: f64 = (0..n).map(|i| q[i * p + r] * v[i]).sum();
                for i in 0..n {
                    v[i] -= dot * q[i * p + r];
                }
            }
        }
        let norm = v.iter().map(|val| val * val).sum::<f64>().sqrt();
        if norm > 1e-12 {
            for i in 0..n {
                q[i * p + rank] = v[i] / norm;
            }
            rank += 1;
        }
    }
    if rank < p {
        let mut compact = vec![0.0; n * rank];
        for i in 0..n {
            compact[i * rank..(i + 1) * rank].copy_from_slice(&q[i * p..i * p + rank]);
        }
        compact
    } else {
        q
    }
}

fn load_bench_csv(name: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let path = format!("{}/tests/fixtures/{}", env!("CARGO_MANIFEST_DIR"), name);
    let text = std::fs::read_to_string(path).expect("fixture csv");
    let mut y = Vec::new();
    let mut x = Vec::new();
    let mut coords = Vec::new();
    for line in text.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<f64> = line
            .split(',')
            .map(|s| s.trim().parse::<f64>().unwrap())
            .collect();
        y.push(cols[1]);
        x.extend_from_slice(&[1.0, cols[2], cols[3]]);
        coords.extend_from_slice(&[cols[6], cols[7]]);
    }
    (y, x, coords)
}

fn prepare_bench_ols(
    y: &[f64],
    x: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, usize) {
    let n = y.len();
    let p = 3usize;
    let mut xtx = vec![0.0; p * p];
    let mut xty = vec![0.0; p];
    for i in 0..n {
        for r in 0..p {
            xty[r] += x[i * p + r] * y[i];
            for c in 0..p {
                xtx[r * p + c] += x[i * p + r] * x[i * p + c];
            }
        }
    }
    let inv_xtx = invert_symmetric(&xtx, p);
    let mut coef = xty.clone();
    let mut a = xtx.clone();
    solve_linear(&mut a, &mut coef, p);

    let mut wj = vec![0.0; n * p];
    let mut xj = vec![0.0; n * p];
    let mut xjs = vec![0.0; n * p];
    for i in 0..n {
        let yhat: f64 = (0..p).map(|k| x[i * p + k] * coef[k]).sum();
        let resid = y[i] - yhat;
        for j in 0..p {
            let x_dir: f64 = (0..p)
                .map(|k| (n as f64) * inv_xtx[j * p + k] * x[i * p + k])
                .sum();
            xj[i * p + j] = x_dir;
            xjs[i * p + j] = if x_dir > 0.0 {
                1.0
            } else if x_dir < 0.0 {
                -1.0
            } else {
                0.0
            };
            wj[i * p + j] = x_dir * resid + coef[j];
        }
    }
    let q_mat = modified_gram_schmidt(x, n, p);
    let proj_rank = q_mat.len() / n;
    (coef, wj, xj, xjs, q_mat, proj_rank)
}

#[test]
fn conditional_basis_weights_apply_regressor_before_projection() {
    let x = [-2.0, -1.0, 1.0, 2.0];
    let sign = [-1.0, -1.0, 1.0, 1.0];
    let norm = 10.0_f64.sqrt();
    let w: Vec<_> = x.iter().flat_map(|value| [0.5, -value / norm]).collect();
    let projection: Vec<_> = x.iter().flat_map(|value| [0.5, value / norm]).collect();
    let a = 1.5 / norm;
    let expected = [1.0, a, 0.5, -a, 0.5, a, 1.0, -a];
    for scale in [1.0, 1e160, 1e-200] {
        let weights = build_conditional_basis(
            &w,
            4,
            1,
            None,
            None,
            &x.map(|v| v * scale),
            &sign,
            1,
            0,
            Some(&projection),
            Some(&projection),
            2,
        )
        .unwrap();
        for (got, expected) in weights.iter().zip(expected) {
            assert!((got / scale - expected).abs() < 1e-12);
        }
        let weights = normalize_covariance_basis(weights);
        assert!(weights.iter().all(|v| v.is_finite() && v.abs() <= 1.0));
    }
}

#[test]
fn score_projection_restores_units_after_extreme_rescaling() {
    let wfin = [0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5];
    let wj = [1.0, 1.0, 1.0, -1.0];
    let (se, t) = compute_uncond_stats(&wfin, &wj, 4, 1, &[0.5]).unwrap();
    for scale in [0.7 * f64::MAX, 1e-300] {
        let scaled = wj.map(|v| v * scale);
        let (s, u) = compute_uncond_stats(&wfin, &scaled, 4, 1, &[0.5 * scale]).unwrap();
        assert_eq!(s[0] / scale, se[0]);
        assert_eq!(u, t);
    }
}

#[test]
fn scpc_scale_grid_includes_the_finite_endpoint_when_needed() {
    assert_eq!(build_scale_grid(1.0, 1.3), [1.0, 1.2, 1.3]);
    assert_eq!(build_scale_grid(1.0, 0.8), [1.0, 1.2]);
}

#[test]
fn rejection_probability_respects_exchangeability_and_strict_singular_boundary() {
    let mut ws = ProbabilityWorkspace::new(1);
    for r in [1.0 - 5e-15, -1.0 + f64::EPSILON, 0.9] {
        assert!((ws.evaluate(&[1.0, r, r, 1.0], 2, 1.0).unwrap() - 0.5).abs() < 1e-15);
    }
    for r in [-1.0, 1.0] {
        let omega = [1.0, r, r, 1.0];
        assert_eq!(ws.evaluate(&omega, 2, 1.0).unwrap(), 0.0);
        assert_eq!(ws.evaluate(&omega, 2, 1.0 - 1e-15).unwrap(), 1.0);
        assert_eq!(ws.evaluate(&omega, 2, 1.0 + 1e-15).unwrap(), 0.0);
    }
    assert!(ws.evaluate(&[0.0, 1.0, 1.0, 0.0], 2, 1.0).is_err());
    assert_eq!(ws.evaluate(&[0.0; 4], 2, 0.0).unwrap(), 0.0);
    for scale in [2.0_f64.powi(-900), 1.0, 2.0_f64.powi(900)] {
        let a = 2401.0 * scale;
        let b = 49.0 * scale;
        assert_eq!(ws.evaluate(&[a, b, b, scale], 2, 49.0).unwrap(), 0.0);
        let near = f64::from_bits(b.to_bits() - 1);
        assert_eq!(ws.evaluate(&[a, near, near, scale], 2, 49.0).unwrap(), 0.5);
    }
    let mut ws = ProbabilityWorkspace::new(2);
    let omega = [25., 15., 20., 15., 9., 12., 20., 12., 16.];
    assert_eq!(ws.evaluate(&omega, 3, 1.0).unwrap(), 0.0);
    assert_eq!(ws.evaluate(&omega, 3, 1.0 - 1e-15).unwrap(), 1.0);
    let mut nearly_rank_one = omega;
    nearly_rank_one[0] = f64::from_bits(25.0_f64.to_bits() + 1);
    let p = ws.evaluate(&nearly_rank_one, 3, 1.0).unwrap();
    assert!((p - 0.5).abs() < 2e-7, "positive direction lost: p={p}");
    let r = 1.0 - 5e-15;
    let padded = [1.0, 0.0, r, 0.0, 0.0, 0.0, r, 0.0, 1.0];
    assert_eq!(ws.evaluate(&padded, 3, 1.0).unwrap(), 0.5);
}

#[test]
fn identity_covariance_probability_matches_student_t_and_near_zero_tail() {
    for q in [1usize, 2, 5, 10] {
        let width = q + 1;
        let mut identity = vec![0.0; width * width];
        for i in 0..width {
            identity[i * width + i] = 1.0;
        }
        let mut ws = ProbabilityWorkspace::new(q);
        for t in [1e-6, 0.1, 1.0, 2.5, 10.0] {
            let cv_norm = t / (q as f64).sqrt();
            let got = ws.evaluate(&identity, width, cv_norm).unwrap();
            let expected = if t * t < q as f64 {
                1.0 - betainc(0.5, 0.5 * q as f64, (t * t) / (q as f64 + t * t))
            } else {
                betainc(0.5 * q as f64, 0.5, (q as f64) / (q as f64 + t * t))
            };
            assert!(
                (got - expected).abs() < 1e-10,
                "q={q} t={t}: got={got}, expected={expected}"
            );
        }
    }

    // Near-zero test statistic reproducer from rust_core_handoff.md:
    // coords=[[0.0],[1.0]], coef=[1e-6], wj=[[1],[-1]], avc=0.03, uncond=true
    let coords = [0.0, 1.0];
    let coef = [1e-6];
    let wj = [1.0, -1.0];
    let out = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n: 2,
        d: 1,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: true,
        cvs: false,
        coef: &coef,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: None,
        xjs: None,
        proj_u: None,
        proj_v: None,
        proj_rank: 0,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    assert_eq!(out.q, 1);
    let rho = 0.03_f64;
    let expected_p = 1.0 - FRAC_2_PI * (1e-6 * ((1.0 - rho) / (1.0 + rho)).sqrt()).atan();
    assert!(
        (out.scpcstats[3] - expected_p).abs() < 1e-12,
        "got p={}, expected={expected_p}",
        out.scpcstats[3]
    );
}

#[test]
fn rejection_probability_is_invariant_to_large_finite_covariance_scale() {
    let mut ws = ProbabilityWorkspace::new(2);
    let scale = 1.4e308;
    let omega = [scale, 0.0, 0.0, 0.0, scale, 0.0, 0.0, 0.0, scale];
    let expected = 1.0 - 1.0 / 2.0_f64.sqrt();
    let actual = ws.evaluate(&omega, 3, 1.0).unwrap();
    assert!((actual - expected).abs() < 1e-12);
}

#[test]
fn exact_scpc_matches_release_score_statistics_and_is_scale_invariant() {
    let (y, x, coords) = load_bench_csv("bench_200.csv");
    let n = y.len();
    let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);

    let out_euclid = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: false,
        cvs: true,
        coef: &coef,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank,
        cv_tolerance: 0.001,
        nthreads: Some(2),
    })
    .unwrap();

    assert_eq!(out_euclid.q, 19);
    assert!(
        (out_euclid.c0 - 3.500409719601247).abs() < 1e-3,
        "c0={}",
        out_euclid.c0
    );
    let dist_euclid = distance_matrix(&coords, n, 2, Metric::Euclidean).unwrap();
    let dvec_euclid = lower_triangle(&dist_euclid, n);
    assert!(
        (crate::calibration::avg_correlation(out_euclid.c0, &dvec_euclid) - 0.03).abs() < 1e-12
    );
    assert!(
        (out_euclid.cvfin - 2.1008365544083087).abs() < 1e-3,
        "cvfin={}",
        out_euclid.cvfin
    );
    let expected_stats = [
        [1.4834167407175964, 0.11454779987636342, 12.95019845268713],
        [1.5473280763431876, 0.09650902048915827, 16.032989128897157],
        [-0.6589528089605742, 0.10738682321940547, -6.13625386435211],
    ];
    for (row, expected) in out_euclid.scpcstats.chunks_exact(6).zip(expected_stats) {
        for (got, exp) in row.iter().zip(expected) {
            assert!(
                (got - exp).abs() <= 1e-3 * exp.abs().max(1.0),
                "got={got}, exp={exp}"
            );
        }
    }
    // Check 5% CV reuse consistency between scpcstats CI and scpccvs[, 2]
    let cvs = out_euclid.scpccvs.as_ref().unwrap();
    for j in 0..3 {
        let se = out_euclid.scpcstats[j * 6 + 1];
        let half_w = out_euclid.scpcstats[j * 6 + 5] - out_euclid.scpcstats[j * 6];
        assert!((half_w - cvs[j * 4 + 2] * se).abs() < 1e-12);
    }

    // Single-threaded serial execution on caller thread matches multi-threaded exactly
    let out_one = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: false,
        cvs: false,
        coef: &coef,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    assert_eq!(out_one.q, out_euclid.q);
    for (a, b) in out_one.scpcstats.iter().zip(out_euclid.scpcstats.iter()) {
        assert!((a - b).abs() < 1e-12);
    }

    // Positive outcome rescaling (`y * 1e-12`) preserves q, t, p and scales SE, CI proportionally
    let y_scale = 1e-12;
    let y_scaled: Vec<f64> = y.iter().map(|&v| v * y_scale).collect();
    let (coef_ys, wj_ys, _, _, _, _) = prepare_bench_ols(&y_scaled, &x);
    let out_ys = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: false,
        cvs: true,
        coef: &coef_ys,
        wj: &wj_ys,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    assert_eq!(out_ys.q, out_euclid.q);
    for j in 0..3 {
        assert!(
            (out_ys.scpcstats[j * 6 + 1] - out_euclid.scpcstats[j * 6 + 1] * y_scale).abs() < 1e-22
        );
        assert!((out_ys.scpcstats[j * 6 + 2] - out_euclid.scpcstats[j * 6 + 2]).abs() < 1e-10);
        assert!((out_ys.scpcstats[j * 6 + 3] - out_euclid.scpcstats[j * 6 + 3]).abs() < 1e-10);
        assert!(
            (out_ys.scpcstats[j * 6 + 4] - out_euclid.scpcstats[j * 6 + 4] * y_scale).abs() < 1e-22
        );
        assert!(
            (out_ys.scpcstats[j * 6 + 5] - out_euclid.scpcstats[j * 6 + 5] * y_scale).abs() < 1e-22
        );
    }

    // Coordinate rescaling (`coords * 1e-6`) preserves q, se, t, p and scales c0 by 1e6
    let c_scale = 1e-6;
    let coords_scaled: Vec<f64> = coords.iter().map(|&v| v * c_scale).collect();
    let out_cs = run_exact_scpc(&ExactScpcRequest {
        coords: &coords_scaled,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: false,
        cvs: false,
        coef: &coef,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank,
        cv_tolerance: 0.001,
        nthreads: Some(2),
    })
    .unwrap();
    assert_eq!(out_cs.q, out_euclid.q);
    assert!((out_cs.c0 * c_scale - out_euclid.c0).abs() / out_euclid.c0 < 1e-10);
    for (a, b) in out_cs.scpcstats.iter().zip(out_euclid.scpcstats.iter()) {
        assert!((a - b).abs() < 1e-10, "a={a}, b={b}");
    }

    // Exact zero wj produces structured DegenerateVariance error
    let wj_zero = vec![0.0; n * 3];
    let err = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: true,
        cvs: false,
        coef: &coef,
        wj: &wj_zero,
        cluster_id: None,
        fe_id: None,
        xj: None,
        xjs: None,
        proj_u: None,
        proj_v: None,
        proj_rank: 0,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap_err();
    assert_eq!(err.kind(), ErrorKind::DegenerateVariance);

    // Latlong conditional (`s_1 = lon, s_2 = lat`)
    let out_geo = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::GreatCircleFractionLonLat,
        avc: 0.03,
        uncond: false,
        cvs: true,
        coef: &coef,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank,
        cv_tolerance: 0.001,
        nthreads: Some(2),
    })
    .unwrap();
    assert_eq!(out_geo.q, 18);
}

#[test]
fn exact_scpc_matches_release_score_statistics_bench_1000() {
    let (y, x, coords) = load_bench_csv("bench_1000.csv");
    let n = y.len();
    let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);

    let out = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: false,
        cvs: true,
        coef: &coef,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank,
        cv_tolerance: 0.001,
        nthreads: Some(4),
    })
    .unwrap();

    assert_eq!(out.q, 8);
    let dist_1000 = distance_matrix(&coords, n, 2, Metric::Euclidean).unwrap();
    let dvec_1000 = lower_triangle(&dist_1000, n);
    assert!((crate::calibration::avg_correlation(out.c0, &dvec_1000) - 0.03).abs() < 1e-12);
    assert!((out.c0 - 0.43896032977784594).abs() < 1e-3, "c0={}", out.c0);
    assert!(
        (out.cvfin - 2.365574447704166).abs() < 1e-3,
        "cvfin={}",
        out.cvfin
    );
    let expected_stats = [
        [1.5780823405600297, 0.13216765018777865, 11.940004519396016],
        [1.745699446088806, 0.0697326838936868, 25.03416401913148],
        [-0.9039531747843945, 0.037101498337623266, -24.3643306951765],
    ];
    for (row, expected) in out.scpcstats.chunks_exact(6).zip(expected_stats) {
        assert!((row[0] - expected[0]).abs() < 1e-12);
        for (got, exp) in row.iter().zip(expected) {
            assert!(
                (got - exp).abs() <= 2e-3 * exp.abs().max(1.0),
                "got={got}, exp={exp}"
            );
        }
    }
}

#[test]
fn cluster_id_five_site_cutoffs_and_singleton_equivalence() {
    // Five-site clustered Gaussian cutoffs from test_conditional_ols.do:
    // s = [0,1,2,4,6], 2 obs/cluster, x_g = {x_0 - 0.25, x_0 + 0.25}, y = x + (s == 0)
    let s_sites = [0.0, 1.0, 2.0, 4.0, 6.0];
    let x0 = [-1.0, -4.0, 2.0, 2.0, 1.0];
    let mut s_obs = Vec::new();
    let mut x_obs = Vec::new();
    let mut y_obs = Vec::new();
    let mut cl_id = Vec::new();
    for (g, (&sg, &xg)) in s_sites.iter().zip(x0.iter()).enumerate() {
        for delta in [-0.25_f64, 0.25_f64] {
            let xval = xg + delta;
            s_obs.push(sg);
            x_obs.push(xval);
            y_obs.push(xval + if sg == 0.0 { 1.0 } else { 0.0 });
            cl_id.push(g);
        }
    }
    let g_sites = s_sites.len();
    let m_obs = y_obs.len();

    // Case 1: With intercept (y ~ x), cluster-aggregated scores, n_eff = G = 5
    let mut x_int = Vec::with_capacity(m_obs * 2);
    for &xv in &x_obs {
        x_int.extend_from_slice(&[1.0, xv]);
    }
    let mut xtx = [0.0; 4];
    let mut xty = [0.0; 2];
    for i in 0..m_obs {
        for r in 0..2 {
            xty[r] += x_int[i * 2 + r] * y_obs[i];
            for c in 0..2 {
                xtx[r * 2 + c] += x_int[i * 2 + r] * x_int[i * 2 + c];
            }
        }
    }
    let bread = invert_symmetric(&xtx, 2);
    let beta: Vec<f64> = (0..2)
        .map(|r| (0..2).map(|c| bread[r * 2 + c] * xty[c]).sum())
        .collect();
    let mut s_cl = vec![0.0; g_sites * 2];
    for i in 0..m_obs {
        let resid = y_obs[i] - (beta[0] + beta[1] * x_obs[i]);
        let g = cl_id[i];
        s_cl[g * 2] += x_int[i * 2] * resid;
        s_cl[g * 2 + 1] += x_int[i * 2 + 1] * resid;
    }
    let mut wj_int = vec![0.0; g_sites];
    for g in 0..g_sites {
        let mut dot = 0.0;
        for c in 0..2 {
            dot += bread[2 + c] * s_cl[g * 2 + c];
        }
        wj_int[g] = (g_sites as f64) * dot + beta[1];
    }
    let mut xj_int = vec![0.0; m_obs];
    let mut sq_sum_int = vec![0.0; g_sites];
    for i in 0..m_obs {
        let val = (g_sites as f64) * (bread[2] + bread[3] * x_obs[i]);
        xj_int[i] = val;
        sq_sum_int[cl_id[i]] += val * val;
    }
    let xjs_int: Vec<f64> = (0..m_obs)
        .map(|i| {
            let d = sq_sum_int[cl_id[i]].sqrt();
            if d > 0.0 {
                xj_int[i] / d
            } else {
                0.0
            }
        })
        .collect();
    let q_int = modified_gram_schmidt(&x_int, m_obs, 2);
    let out_int = run_exact_scpc(&ExactScpcRequest {
        coords: &s_sites,
        n: g_sites,
        d: 1,
        metric: Metric::Euclidean,
        avc: 0.5,
        uncond: false,
        cvs: true,
        coef: &[beta[1]],
        wj: &wj_int,
        cluster_id: Some(&cl_id),
        fe_id: None,
        xj: Some(&xj_int),
        xjs: Some(&xjs_int),
        proj_u: Some(&q_int),
        proj_v: Some(&q_int),
        proj_rank: 2,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    let cv5_int = out_int.scpccvs.as_ref().unwrap()[2];
    assert!(
        (cv5_int - 17.380977).abs() < 0.002,
        "clustered int cv5={cv5_int}"
    );

    // Case 2: No intercept (y ~ x - 1)
    let denom_noint: f64 = x_obs.iter().map(|v| v * v).sum();
    let beta_noint: f64 = x_obs
        .iter()
        .zip(y_obs.iter())
        .map(|(xv, yv)| xv * yv)
        .sum::<f64>()
        / denom_noint;
    let mut s_cl_no = vec![0.0; g_sites];
    for i in 0..m_obs {
        let resid = y_obs[i] - beta_noint * x_obs[i];
        s_cl_no[cl_id[i]] += x_obs[i] * resid;
    }
    let wj_no: Vec<f64> = s_cl_no
        .iter()
        .map(|sg| (g_sites as f64) * (sg / denom_noint) + beta_noint)
        .collect();
    let xj_no: Vec<f64> = x_obs
        .iter()
        .map(|xv| (g_sites as f64) * (xv / denom_noint))
        .collect();
    let mut sq_sum_no = vec![0.0; g_sites];
    for i in 0..m_obs {
        sq_sum_no[cl_id[i]] += xj_no[i] * xj_no[i];
    }
    let xjs_no: Vec<f64> = (0..m_obs)
        .map(|i| {
            let d = sq_sum_no[cl_id[i]].sqrt();
            if d > 0.0 {
                xj_no[i] / d
            } else {
                0.0
            }
        })
        .collect();
    let q_no = modified_gram_schmidt(&x_obs, m_obs, 1);
    let out_no = run_exact_scpc(&ExactScpcRequest {
        coords: &s_sites,
        n: g_sites,
        d: 1,
        metric: Metric::Euclidean,
        avc: 0.5,
        uncond: false,
        cvs: true,
        coef: &[beta_noint],
        wj: &wj_no,
        cluster_id: Some(&cl_id),
        fe_id: None,
        xj: Some(&xj_no),
        xjs: Some(&xjs_no),
        proj_u: Some(&q_no),
        proj_v: Some(&q_no),
        proj_rank: 1,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    let cv5_no = out_no.scpccvs.as_ref().unwrap()[2];
    assert!(
        (cv5_no - 9.902023).abs() < 0.002,
        "clustered noint cv5={cv5_no}"
    );

    // Case 3: Singleton clusters (cluster_id = 0..n-1) match unclustered cluster_id = None to < 1e-12
    let (y, x, coords) = load_bench_csv("bench_200.csv");
    let n = y.len();
    let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);
    let singleton: Vec<usize> = (0..n).collect();
    let out_none = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: false,
        cvs: true,
        coef: &coef,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    let out_sing = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: false,
        cvs: true,
        coef: &coef,
        wj: &wj,
        cluster_id: Some(&singleton),
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    for (a, b) in out_none.scpcstats.iter().zip(out_sing.scpcstats.iter()) {
        assert!((a - b).abs() < 1e-12, "a={a}, b={b}");
    }
    for (a, b) in out_none
        .scpccvs
        .unwrap()
        .iter()
        .zip(out_sing.scpccvs.unwrap().iter())
    {
        assert!((a - b).abs() < 1e-12, "cv a={a}, b={b}");
    }
}

#[test]
fn exact_scpc_iv_and_one_way_fe_match_explicit_dummies() {
    let (_, x_all, _) = load_bench_csv("bench_200.csv");
    let n_obs = 60usize;
    let g_sites = 30usize;
    let n_fe = 5usize;
    // Staggered clusters (pairs of obs) and fixed effects (crossing cluster boundaries)
    let cl_id: Vec<usize> = (0..n_obs).map(|i| i / 2).collect();
    let fe_id: Vec<usize> = (0..n_obs).map(|i| i % n_fe).collect();
    let mut s_cl = Vec::with_capacity(g_sites * 2);
    for g in 0..g_sites {
        s_cl.push((g as f64) * 0.13 + ((g % 5) as f64) * 0.21);
        s_cl.push(((g % 7) as f64) * 0.17 - (g as f64) * 0.04);
    }
    let s_obs: Vec<f64> = cl_id
        .iter()
        .enumerate()
        .flat_map(|(i, &g)| [s_cl[g * 2] + (i as f64) * 0.01, s_cl[g * 2 + 1]])
        .collect();

    // Build exogenous w, excluded instrument z, endogenous x1, and outcome y
    let mut w = Vec::with_capacity(n_obs);
    let mut z = Vec::with_capacity(n_obs);
    let mut x1 = Vec::with_capacity(n_obs);
    let mut y = Vec::with_capacity(n_obs);
    for i in 0..n_obs {
        let wi = x_all[i * 3 + 1];
        let zi = x_all[i * 3 + 2] + 0.3 * (i as f64).sin();
        let ui = 0.4 * ((i as f64) * 1.1).cos() + 0.15 * (fe_id[i] as f64);
        let xi = 0.85 * zi + 0.35 * wi + 0.55 * ui + 0.2 * ((fe_id[i] as f64) - 2.0);
        let yi = 1.2 * xi + 0.6 * wi + 0.9 * (fe_id[i] as f64) + ui;
        w.push(wi);
        z.push(zi);
        x1.push(xi);
        y.push(yi);
    }

    let demean = |col: &[f64]| -> Vec<f64> {
        let mut sums = vec![0.0; n_fe];
        let mut cnts = vec![0usize; n_fe];
        for (i, &val) in col.iter().enumerate() {
            sums[fe_id[i]] += val;
            cnts[fe_id[i]] += 1;
        }
        col.iter()
            .enumerate()
            .map(|(i, &val)| val - sums[fe_id[i]] / (cnts[fe_id[i]] as f64))
            .collect()
    };

    // --- Part 1: Unclustered & Clustered 1-way FE OLS (y ~ x1 + w | fe) ---
    let x_dm_c0 = demean(&x1);
    let x_dm_c1 = demean(&w);
    let y_dm = demean(&y);
    let mut x_dm = vec![0.0; n_obs * 2];
    let mut x_full_dum = vec![0.0; n_obs * (2 + n_fe)];
    for i in 0..n_obs {
        x_dm[i * 2] = x_dm_c0[i];
        x_dm[i * 2 + 1] = x_dm_c1[i];
        x_full_dum[i * (2 + n_fe)] = x1[i];
        x_full_dum[i * (2 + n_fe) + 1] = w[i];
        x_full_dum[i * (2 + n_fe) + 2 + fe_id[i]] = 1.0;
    }
    let solve_ols =
        |x_mat: &[f64], y_vec: &[f64], p_cols: usize| -> (Vec<f64>, Vec<f64>, Vec<f64>) {
            let mut xtx = vec![0.0; p_cols * p_cols];
            let mut xty = vec![0.0; p_cols];
            for i in 0..n_obs {
                for r in 0..p_cols {
                    xty[r] += x_mat[i * p_cols + r] * y_vec[i];
                    for c in 0..p_cols {
                        xtx[r * p_cols + c] += x_mat[i * p_cols + r] * x_mat[i * p_cols + c];
                    }
                }
            }
            let inv = invert_symmetric(&xtx, p_cols);
            let mut beta = xty.clone();
            let mut a = xtx;
            solve_linear(&mut a, &mut beta, p_cols);
            let bread: Vec<f64> = inv.iter().map(|v| v * (n_obs as f64)).collect();
            let mut res = vec![0.0; n_obs];
            for i in 0..n_obs {
                let fit: f64 = (0..p_cols).map(|c| x_mat[i * p_cols + c] * beta[c]).sum();
                res[i] = y_vec[i] - fit;
            }
            (beta, bread, res)
        };
    let (beta_dm, bread_dm, res_dm) = solve_ols(&x_dm, &y_dm, 2);
    let mut wj_ols_u = vec![0.0; n_obs * 2];
    let mut xj_ols = vec![0.0; n_obs * 2];
    for i in 0..n_obs {
        for j in 0..2 {
            let bdot: f64 = (0..2).map(|c| bread_dm[j * 2 + c] * x_dm[i * 2 + c]).sum();
            wj_ols_u[i * 2 + j] = (n_obs as f64) * res_dm[i] * bdot + beta_dm[j];
            xj_ols[i * 2 + j] = (n_obs as f64) * bdot;
        }
    }
    let xjs_ols_u: Vec<f64> = xj_ols.iter().map(|&v| v.signum()).collect();
    let q_dm = modified_gram_schmidt(&x_dm, n_obs, 2);
    let q_dum = modified_gram_schmidt(&x_full_dum, n_obs, 2 + n_fe);

    let out_fe_ols = run_exact_scpc(&ExactScpcRequest {
        coords: &s_obs,
        n: n_obs,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.1,
        uncond: false,
        cvs: true,
        coef: &beta_dm,
        wj: &wj_ols_u,
        cluster_id: None,
        fe_id: Some(&fe_id),
        xj: Some(&xj_ols),
        xjs: Some(&xjs_ols_u),
        proj_u: Some(&q_dm),
        proj_v: Some(&q_dm),
        proj_rank: 2,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    let out_dum_ols = run_exact_scpc(&ExactScpcRequest {
        coords: &s_obs,
        n: n_obs,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.1,
        uncond: false,
        cvs: true,
        coef: &beta_dm,
        wj: &wj_ols_u,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj_ols),
        xjs: Some(&xjs_ols_u),
        proj_u: Some(&q_dum),
        proj_v: Some(&q_dum),
        proj_rank: 2 + n_fe,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();
    assert_eq!(out_fe_ols.q, out_dum_ols.q);
    for (a, b) in out_fe_ols
        .scpcstats
        .iter()
        .zip(out_dum_ols.scpcstats.iter())
    {
        assert!((a - b).abs() < 1e-11, "uncl fe ols stats a={a}, b={b}");
    }
    for (a, b) in out_fe_ols
        .scpccvs
        .unwrap()
        .iter()
        .zip(out_dum_ols.scpccvs.unwrap().iter())
    {
        assert!((a - b).abs() < 1e-11, "uncl fe ols cvs a={a}, b={b}");
    }

    // Clustered 1-way FE OLS
    let mut s_cl_ols = vec![0.0; g_sites * 2];
    for i in 0..n_obs {
        let g = cl_id[i];
        for c in 0..2 {
            s_cl_ols[g * 2 + c] += x_dm[i * 2 + c] * res_dm[i];
        }
    }
    let mut wj_ols_c = vec![0.0; g_sites * 2];
    for g in 0..g_sites {
        for j in 0..2 {
            let mut dot = 0.0;
            for c in 0..2 {
                dot += (bread_dm[j * 2 + c] / (n_obs as f64)) * s_cl_ols[g * 2 + c];
            }
            wj_ols_c[g * 2 + j] = (g_sites as f64) * dot + beta_dm[j];
        }
    }
    let xj_ols_c: Vec<f64> = (0..n_obs)
        .flat_map(|i| {
            let x_dm = &x_dm;
            let bread_dm = &bread_dm;
            (0..2).map(move |j| {
                let bdot: f64 = (0..2)
                    .map(|c| (bread_dm[j * 2 + c] / (n_obs as f64)) * x_dm[i * 2 + c])
                    .sum();
                (g_sites as f64) * bdot
            })
        })
        .collect();
    let mut sq_ols_c = vec![0.0; g_sites * 2];
    for i in 0..n_obs {
        for j in 0..2 {
            sq_ols_c[cl_id[i] * 2 + j] += xj_ols_c[i * 2 + j].powi(2);
        }
    }
    let xjs_ols_c: Vec<f64> = (0..n_obs)
        .flat_map(|i| {
            let cl_id = &cl_id;
            let xj_ols_c = &xj_ols_c;
            let sq_ols_c = &sq_ols_c;
            (0..2).map(move |j| {
                let d = sq_ols_c[cl_id[i] * 2 + j].sqrt();
                if d > 0.0 {
                    xj_ols_c[i * 2 + j] / d
                } else {
                    0.0
                }
            })
        })
        .collect();

    let out_fe_ols_c = run_exact_scpc(&ExactScpcRequest {
        coords: &s_cl,
        n: g_sites,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.1,
        uncond: false,
        cvs: true,
        coef: &beta_dm,
        wj: &wj_ols_c,
        cluster_id: Some(&cl_id),
        fe_id: Some(&fe_id),
        xj: Some(&xj_ols_c),
        xjs: Some(&xjs_ols_c),
        proj_u: Some(&q_dm),
        proj_v: Some(&q_dm),
        proj_rank: 2,
        cv_tolerance: 0.001,
        nthreads: Some(2),
    })
    .unwrap();
    let out_dum_ols_c = run_exact_scpc(&ExactScpcRequest {
        coords: &s_cl,
        n: g_sites,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.1,
        uncond: false,
        cvs: true,
        coef: &beta_dm,
        wj: &wj_ols_c,
        cluster_id: Some(&cl_id),
        fe_id: None,
        xj: Some(&xj_ols_c),
        xjs: Some(&xjs_ols_c),
        proj_u: Some(&q_dum),
        proj_v: Some(&q_dum),
        proj_rank: 2 + n_fe,
        cv_tolerance: 0.001,
        nthreads: Some(2),
    })
    .unwrap();
    for (a, b) in out_fe_ols_c
        .scpcstats
        .iter()
        .zip(out_dum_ols_c.scpcstats.iter())
    {
        assert!((a - b).abs() < 1e-11, "cl fe ols stats a={a}, b={b}");
    }

    // --- Part 2: Unclustered & Clustered 1-way FE IV (y ~ w | fe | x1 ~ z) ---
    let z_dm_c0 = demean(&z);
    let mut z_dm = vec![0.0; n_obs * 2];
    let mut z_full_dum = vec![0.0; n_obs * (2 + n_fe)];
    for i in 0..n_obs {
        z_dm[i * 2] = z_dm_c0[i];
        z_dm[i * 2 + 1] = x_dm_c1[i];
        z_full_dum[i * (2 + n_fe)] = z[i];
        z_full_dum[i * (2 + n_fe) + 1] = w[i];
        z_full_dum[i * (2 + n_fe) + 2 + fe_id[i]] = 1.0;
    }
    let build_iv_proj =
        |x_mat: &[f64], z_mat: &[f64], p_cols: usize| -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
            let qz = modified_gram_schmidt(z_mat, n_obs, p_cols);
            let rz = qz.len() / n_obs;
            let mut qzt_x = vec![0.0; rz * p_cols];
            linalg::multiply(&qz, x_mat, &mut qzt_x, rz, n_obs, p_cols, true);
            let mut pzx = vec![0.0; n_obs * p_cols];
            linalg::multiply(&qz, &qzt_x, &mut pzx, n_obs, rz, p_cols, false);
            let mut a = vec![0.0; p_cols * p_cols];
            linalg::multiply(x_mat, &pzx, &mut a, p_cols, n_obs, p_cols, true);
            let a_inv = invert_symmetric(&a, p_cols);
            let mut proj_v = vec![0.0; n_obs * p_cols];
            linalg::multiply(&pzx, &a_inv, &mut proj_v, n_obs, p_cols, p_cols, false);
            (pzx, a_inv, x_mat.to_vec(), proj_v)
        };
    let (pzx_dm, a_inv_dm, u_iv_dm, v_iv_dm) = build_iv_proj(&x_dm, &z_dm, 2);
    let (_, _, u_iv_dum, v_iv_dum) = build_iv_proj(&x_full_dum, &z_full_dum, 2 + n_fe);

    // Compute IV coefficients & influence functions
    let mut pzy = vec![0.0; 2];
    for i in 0..n_obs {
        pzy[0] += pzx_dm[i * 2] * y_dm[i];
        pzy[1] += pzx_dm[i * 2 + 1] * y_dm[i];
    }
    let mut beta_iv = vec![0.0; 2];
    for j in 0..2 {
        beta_iv[j] = a_inv_dm[j * 2] * pzy[0] + a_inv_dm[j * 2 + 1] * pzy[1];
    }
    let mut res_iv = vec![0.0; n_obs];
    for i in 0..n_obs {
        res_iv[i] = y_dm[i] - (x_dm[i * 2] * beta_iv[0] + x_dm[i * 2 + 1] * beta_iv[1]);
    }
    let mut wj_iv_c = vec![0.0; g_sites * 2];
    let mut s_cl_iv = vec![0.0; g_sites * 2];
    for i in 0..n_obs {
        let g = cl_id[i];
        for c in 0..2 {
            s_cl_iv[g * 2 + c] += pzx_dm[i * 2 + c] * res_iv[i];
        }
    }
    for g in 0..g_sites {
        for j in 0..2 {
            let dot: f64 = (0..2)
                .map(|c| a_inv_dm[j * 2 + c] * s_cl_iv[g * 2 + c])
                .sum();
            wj_iv_c[g * 2 + j] = (g_sites as f64) * dot + beta_iv[j];
        }
    }
    let xj_iv_c: Vec<f64> = (0..n_obs)
        .flat_map(|i| {
            let pzx_dm = &pzx_dm;
            let a_inv_dm = &a_inv_dm;
            (0..2).map(move |j| {
                let dot: f64 = (0..2)
                    .map(|c| a_inv_dm[j * 2 + c] * pzx_dm[i * 2 + c])
                    .sum();
                (g_sites as f64) * dot
            })
        })
        .collect();
    let mut sq_iv_c = vec![0.0; g_sites * 2];
    for i in 0..n_obs {
        for j in 0..2 {
            sq_iv_c[cl_id[i] * 2 + j] += xj_iv_c[i * 2 + j].powi(2);
        }
    }
    let xjs_iv_c: Vec<f64> = (0..n_obs)
        .flat_map(|i| {
            let cl_id = &cl_id;
            let xj_iv_c = &xj_iv_c;
            let sq_iv_c = &sq_iv_c;
            (0..2).map(move |j| {
                let d = sq_iv_c[cl_id[i] * 2 + j].sqrt();
                if d > 0.0 {
                    xj_iv_c[i * 2 + j] / d
                } else {
                    0.0
                }
            })
        })
        .collect();

    // IV denominator weights are Xi_s' M' Xi R where M' = I - V U'.
    let out_fe_iv_c = run_exact_scpc(&ExactScpcRequest {
        coords: &s_cl,
        n: g_sites,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.1,
        uncond: false,
        cvs: true,
        coef: &beta_iv,
        wj: &wj_iv_c,
        cluster_id: Some(&cl_id),
        fe_id: Some(&fe_id),
        xj: Some(&xj_iv_c),
        xjs: Some(&xjs_iv_c),
        proj_u: Some(&v_iv_dm),
        proj_v: Some(&u_iv_dm),
        proj_rank: 2,
        cv_tolerance: 0.001,
        nthreads: Some(2),
    })
    .unwrap();
    let out_dum_iv_c = run_exact_scpc(&ExactScpcRequest {
        coords: &s_cl,
        n: g_sites,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.1,
        uncond: false,
        cvs: true,
        coef: &beta_iv,
        wj: &wj_iv_c,
        cluster_id: Some(&cl_id),
        fe_id: None,
        xj: Some(&xj_iv_c),
        xjs: Some(&xjs_iv_c),
        proj_u: Some(&v_iv_dum),
        proj_v: Some(&u_iv_dum),
        proj_rank: 2 + n_fe,
        cv_tolerance: 0.001,
        nthreads: Some(2),
    })
    .unwrap();
    assert_eq!(out_fe_iv_c.q, out_dum_iv_c.q);
    for (a, b) in out_fe_iv_c
        .scpcstats
        .iter()
        .zip(out_dum_iv_c.scpcstats.iter())
    {
        assert!((a - b).abs() < 1e-11, "cl fe iv stats a={a}, b={b}");
    }
    for (a, b) in out_fe_iv_c
        .scpccvs
        .unwrap()
        .iter()
        .zip(out_dum_iv_c.scpccvs.unwrap().iter())
    {
        assert!((a - b).abs() < 1e-11, "cl fe iv cvs a={a}, b={b}");
    }
}

fn permute_rows(src: &[f64], perm: &[usize], cols: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(perm.len() * cols);
    for &src_idx in perm {
        out.extend_from_slice(&src[src_idx * cols..(src_idx + 1) * cols]);
    }
    out
}

#[test]
fn nystrom_conditional_results_follow_the_coefficient_across_subsets_and_order() {
    let (y, x, coords) = load_bench_csv("bench_200.csv");
    let n = y.len();
    let (coef, scores, xj, signs, projection, rank) = prepare_bench_ols(&y, &x);
    let (coords, perm) =
        crate::nystrom::normalize_coords(&coords, n, 2, Metric::Euclidean).unwrap();
    let projection = permute_rows(&projection, &perm, rank);
    let original_width = coef.len();
    let fit = |columns: &[usize]| {
        let select = |a: &[f64]| -> Vec<f64> {
            perm.iter()
                .flat_map(|&i| columns.iter().map(move |&j| a[i * original_width + j]))
                .collect()
        };
        let coef_subset: Vec<f64> = columns.iter().map(|&j| coef[j]).collect();
        run_approx_scpc(&ApproxScpcRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef_subset,
            wj: &select(&scores),
            cluster_id: None,
            fe_id: None,
            xj: Some(&select(&xj)),
            xjs: Some(&select(&signs)),
            proj_u: Some(&projection),
            proj_v: Some(&projection),
            proj_rank: rank,
            cv_tolerance: 0.001,
            nthreads: Some(1),
            seed: 1,
            landmarks: 100,
            samples: 3,
            pair_count: 200_000,
        })
        .unwrap()
    };
    let full = fit(&[0, 1, 2]);
    for order in [&[2, 0, 1][..], &[1][..]] {
        let sub = fit(order);
        assert_eq!(sub.q, full.q);
        for (i, &j) in order.iter().enumerate() {
            for (width, a, b) in [
                (6, &sub.scpcstats, &full.scpcstats),
                (
                    4,
                    sub.scpccvs.as_ref().unwrap(),
                    full.scpccvs.as_ref().unwrap(),
                ),
            ] {
                for (left, right) in a[i * width..(i + 1) * width]
                    .iter()
                    .zip(&b[j * width..(j + 1) * width])
                {
                    assert!(
                        (left - right).abs() < 1e-11,
                        "column {j}: {left} vs {right}"
                    );
                }
            }
        }
    }
}

#[test]
fn approx_scpc_matches_release_and_is_thread_count_invariant() {
    let (y, x, coords) = load_bench_csv("bench_1000.csv");
    let n = y.len();
    let p = 3usize;
    let (coef, wj, xj, xjs, q_mat, proj_rank) = prepare_bench_ols(&y, &x);

    for metric in [Metric::Euclidean, Metric::GreatCircleFractionLonLat] {
        let (coords_norm, perm) = crate::nystrom::normalize_coords(&coords, n, 2, metric).unwrap();
        let wj_p = permute_rows(&wj, &perm, p);
        let xj_p = permute_rows(&xj, &perm, p);
        let xjs_p = permute_rows(&xjs, &perm, p);
        let q_p = permute_rows(&q_mat, &perm, proj_rank);

        let out_1 = run_approx_scpc(&ApproxScpcRequest {
            coords: &coords_norm,
            n,
            d: 2,
            metric,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj_p,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj_p),
            xjs: Some(&xjs_p),
            proj_u: Some(&q_p),
            proj_v: Some(&q_p),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(1),
            seed: 1,
            landmarks: 1000,
            samples: 20,
            pair_count: 1_000_000,
        })
        .unwrap();

        let out_4 = run_approx_scpc(&ApproxScpcRequest {
            coords: &coords_norm,
            n,
            d: 2,
            metric,
            avc: 0.03,
            uncond: false,
            cvs: true,
            coef: &coef,
            wj: &wj_p,
            cluster_id: None,
            fe_id: None,
            xj: Some(&xj_p),
            xjs: Some(&xjs_p),
            proj_u: Some(&q_p),
            proj_v: Some(&q_p),
            proj_rank,
            cv_tolerance: 0.001,
            nthreads: Some(4),
            seed: 1,
            landmarks: 1000,
            samples: 20,
            pair_count: 1_000_000,
        })
        .unwrap();

        assert!(!out_1.tied);
        assert!(out_1.q >= 1);
        assert!(out_1.c0 > 0.0 && out_1.cmax >= out_1.c0 && out_1.cvfin > 0.0);
        assert_eq!(out_1.q, out_4.q);
        assert!((out_1.c0 - out_4.c0).abs() < 1e-12);
        assert!((out_1.cmax - out_4.cmax).abs() < 1e-12);
        assert!((out_1.cvfin - out_4.cvfin).abs() < 1e-12);
        for (row1, row4) in out_1
            .scpcstats
            .chunks_exact(6)
            .zip(out_4.scpcstats.chunks_exact(6))
        {
            assert!(row1[1] > 0.0 && row1[4] < row1[0] && row1[0] < row1[5]);
            assert!(row1[3] >= 0.0 && row1[3] <= 1.0);
            for (&a, &b) in row1.iter().zip(row4.iter()) {
                assert!((a - b).abs() < 1e-12, "approx stats mismatch a={a}, b={b}");
            }
        }
        for (&a, &b) in out_1
            .scpccvs
            .as_ref()
            .unwrap()
            .iter()
            .zip(out_4.scpccvs.as_ref().unwrap().iter())
        {
            assert!((a - b).abs() < 1e-12, "approx cvs mismatch a={a}, b={b}");
        }
        for c in 0..=out_1.q {
            let width = out_1.q + 1;
            let column = out_1
                .wfin
                .chunks_exact(width)
                .zip(out_4.wfin.chunks_exact(width));
            let sign = column
                .clone()
                .map(|(a, b)| a[c] * b[c])
                .sum::<f64>()
                .signum();
            for (a, b) in column {
                assert!(
                    (a[c] - sign * b[c]).abs() < 1e-12,
                    "approx basis direction {c}"
                );
            }
        }
    }

    // Unconditional Nystrom run on Euclidean coords is also thread-count invariant.
    let (coords_norm, perm) =
        crate::nystrom::normalize_coords(&coords, n, 2, Metric::Euclidean).unwrap();
    let wj_p = permute_rows(&wj, &perm, p);
    let out_u1 = run_approx_scpc(&ApproxScpcRequest {
        coords: &coords_norm,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: true,
        cvs: true,
        coef: &coef,
        wj: &wj_p,
        cluster_id: None,
        fe_id: None,
        xj: None,
        xjs: None,
        proj_u: None,
        proj_v: None,
        proj_rank: 0,
        cv_tolerance: 0.001,
        nthreads: Some(1),
        seed: 1,
        landmarks: 1000,
        samples: 20,
        pair_count: 1_000_000,
    })
    .unwrap();
    let out_u4 = run_approx_scpc(&ApproxScpcRequest {
        coords: &coords_norm,
        n,
        d: 2,
        metric: Metric::Euclidean,
        avc: 0.03,
        uncond: true,
        cvs: true,
        coef: &coef,
        wj: &wj_p,
        cluster_id: None,
        fe_id: None,
        xj: None,
        xjs: None,
        proj_u: None,
        proj_v: None,
        proj_rank: 0,
        cv_tolerance: 0.001,
        nthreads: Some(4),
        seed: 1,
        landmarks: 1000,
        samples: 20,
        pair_count: 1_000_000,
    })
    .unwrap();
    assert_eq!(out_u1.q, out_u4.q);
    for (&a, &b) in out_u1.scpcstats.iter().zip(out_u4.scpcstats.iter()) {
        assert!((a - b).abs() < 1e-12);
    }
    for (&a, &b) in out_u1
        .scpccvs
        .as_ref()
        .unwrap()
        .iter()
        .zip(out_u4.scpccvs.as_ref().unwrap().iter())
    {
        assert!((a - b).abs() < 1e-12);
    }
}

#[test]
fn iv_conditional_weights_use_transposed_residual_maker_and_match_ols_when_z_equals_x() {
    // Reproducer 1 from rust_core_handoff.md: s=[0,1,2,4,6], x=[-1,-4,2,2,1], z=x, y=x+(s==0).
    let n = 5usize;
    let coords = [0.0, 1.0, 2.0, 4.0, 6.0];
    let x = [-1.0, -4.0, 2.0, 2.0, 1.0];
    let y: Vec<f64> = x
        .iter()
        .enumerate()
        .map(|(i, &xi)| xi + if i == 0 { 1.0 } else { 0.0 })
        .collect();

    // Design [1, x], fit OLS y ~ 1 + x
    let mut xtx = [0.0; 4];
    let mut xty = [0.0; 2];
    for i in 0..n {
        xtx[0] += 1.0;
        xtx[1] += x[i];
        xtx[2] += x[i];
        xtx[3] += x[i] * x[i];
        xty[0] += y[i];
        xty[1] += x[i] * y[i];
    }
    let det = xtx[0] * xtx[3] - xtx[1] * xtx[2];
    let inv = [xtx[3] / det, -xtx[1] / det, -xtx[2] / det, xtx[0] / det];
    let beta = [
        inv[0] * xty[0] + inv[1] * xty[1],
        inv[2] * xty[0] + inv[3] * xty[1],
    ];
    let bread = inv;
    let mut wj = vec![0.0; n * 2];
    let mut xj = vec![0.0; n * 2];
    let mut xjs = vec![0.0; n * 2];
    let mut x_full = vec![0.0; n * 2];
    for i in 0..n {
        let resid = y[i] - (beta[0] + beta[1] * x[i]);
        let row = [1.0, x[i]];
        x_full[i * 2] = 1.0;
        x_full[i * 2 + 1] = x[i];
        for j in 0..2 {
            let d_ij = (n as f64) * (bread[j * 2] * row[0] + bread[j * 2 + 1] * row[1]);
            wj[i * 2 + j] = d_ij * resid + beta[j];
            xj[i * 2 + j] = d_ij;
            xjs[i * 2 + j] = if d_ij > 0.0 {
                1.0
            } else if d_ij < 0.0 {
                -1.0
            } else {
                0.0
            };
        }
    }
    let q_mat = modified_gram_schmidt(&x_full, n, 2);
    let out_ols = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 1,
        metric: Metric::Euclidean,
        avc: 0.5,
        uncond: false,
        cvs: true,
        coef: &beta,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&q_mat),
        proj_v: Some(&q_mat),
        proj_rank: 2,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();

    // For IV with z=x: V = P_Z X (X' P_Z X)^{-1} = X (X' X)^{-1}, U = X.
    let mut v_iv = vec![0.0; n * 2];
    for i in 0..n {
        v_iv[i * 2] = inv[0] + inv[1] * x[i];
        v_iv[i * 2 + 1] = inv[2] + inv[3] * x[i];
    }
    let out_iv = run_exact_scpc(&ExactScpcRequest {
        coords: &coords,
        n,
        d: 1,
        metric: Metric::Euclidean,
        avc: 0.5,
        uncond: false,
        cvs: true,
        coef: &beta,
        wj: &wj,
        cluster_id: None,
        fe_id: None,
        xj: Some(&xj),
        xjs: Some(&xjs),
        proj_u: Some(&v_iv),
        proj_v: Some(&x_full),
        proj_rank: 2,
        cv_tolerance: 0.001,
        nthreads: Some(1),
    })
    .unwrap();

    assert!((out_ols.scpcstats[6 + 3] - 0.0241534748).abs() < 1e-8);
    assert!((out_ols.scpccvs.as_ref().unwrap()[4 + 2] - 17.62678224).abs() < 1e-4);
    for (&a, &b) in out_ols.scpcstats.iter().zip(&out_iv.scpcstats) {
        assert!((a - b).abs() < 1e-11, "OLS vs IV stat mismatch: {a} vs {b}");
    }
    for (&a, &b) in out_ols
        .scpccvs
        .as_ref()
        .unwrap()
        .iter()
        .zip(out_iv.scpccvs.as_ref().unwrap())
    {
        assert!((a - b).abs() < 1e-11, "OLS vs IV cv mismatch: {a} vs {b}");
    }

    // Independent quadrature for this covariance gives 0.07143375445 at
    // the previous IV cutoff. The corrected cutoff must control size.
    let conditional = build_conditional_basis(
        &out_iv.wfin,
        n,
        out_iv.q,
        None,
        None,
        &xj,
        &xjs,
        2,
        1,
        Some(&v_iv),
        Some(&x_full),
        2,
    )
    .unwrap();
    let dist = distance_matrix(&coords, n, 1, Metric::Euclidean).unwrap();
    let omega = covariance_grids(&dist, n, &[&conditional], &[out_iv.q + 1], &[out_iv.c0])
        .unwrap()
        .remove(0);
    let mut probability = ProbabilityWorkspace::new(out_iv.q);
    let root_q = (out_iv.q as f64).sqrt();
    let old_p = probability
        .evaluate(&omega, out_iv.q + 1, 14.50861818 / root_q)
        .unwrap();
    assert!(
        (old_p - 0.07143375445).abs() < 1e-8,
        "old rejection probability {old_p}"
    );
    let corrected_p = probability
        .evaluate(
            &omega,
            out_iv.q + 1,
            out_iv.scpccvs.as_ref().unwrap()[6] / root_q,
        )
        .unwrap();
    assert!(
        corrected_p <= 0.05,
        "corrected rejection probability {corrected_p}"
    );

    // Genuinely oblique IV projection: M = I - U V', check build_conditional_basis == D_xjs M' D_xj B_q
    // U = [1, x], V = P_Z X (X' P_Z X)^{-1} where X is oblique to Z
    let z_ob = [-2.0, -1.0, 0.5, 3.0, -0.5];
    // Center z_ob and x to keep exact orthogonality with [1]
    let mx = x.iter().sum::<f64>() / (n as f64);
    let mz = z_ob.iter().sum::<f64>() / (n as f64);
    let xc: Vec<f64> = x.iter().map(|v| v - mx).collect();
    let zc: Vec<f64> = z_ob.iter().map(|v| v - mz).collect();
    let ztz: f64 = zc.iter().map(|v| v * v).sum();
    let ztx: f64 = zc.iter().zip(&xc).map(|(a, b)| a * b).sum();
    let pzx1: Vec<f64> = zc.iter().map(|&zi| zi * (ztx / ztz)).collect();
    let xt_pzx: f64 = xc.iter().zip(&pzx1).map(|(a, b)| a * b).sum();
    let mut u_ob = vec![0.0; n * 2];
    let mut v_ob = vec![0.0; n * 2];
    let mut w_raw = vec![0.0; n * 2];
    for i in 0..n {
        u_ob[i * 2] = 1.0;
        u_ob[i * 2 + 1] = xc[i];
        v_ob[i * 2] = 1.0 / (n as f64);
        v_ob[i * 2 + 1] = pzx1[i] / xt_pzx;
        w_raw[i * 2] = 1.0 / (n as f64).sqrt();
        w_raw[i * 2 + 1] = ((i as f64) * 1.1).sin();
    }
    // Build expected D_xjs (I - V U') D_xj b_1 directly
    let basis_got = build_conditional_basis(
        &w_raw,
        n,
        1,
        None,
        None,
        &xj,
        &xjs,
        2,
        1,
        Some(&v_ob),
        Some(&u_ob),
        2,
    )
    .unwrap();
    let mut u_t_rx = [0.0; 2];
    for i in 0..n {
        let rx_i = w_raw[i * 2 + 1] * xj[i * 2 + 1];
        u_t_rx[0] += u_ob[i * 2] * rx_i;
        u_t_rx[1] += u_ob[i * 2 + 1] * rx_i;
    }
    for i in 0..n {
        let rx_i = w_raw[i * 2 + 1] * xj[i * 2 + 1];
        let m_t_rx_i = rx_i - (v_ob[i * 2] * u_t_rx[0] + v_ob[i * 2 + 1] * u_t_rx[1]);
        let exp_i = m_t_rx_i * xjs[i * 2 + 1];
        assert!((basis_got[i * 2 + 1] - exp_i).abs() < 1e-12);
    }
}

#[test]
fn probability_workspace_handles_equicorrelated_zero_cv_and_extreme_scales() {
    // 1. Equicorrelated 5x5 covariance (rho = 0.5) across extreme cv values.
    let n5 = 5usize;
    let mut eq = vec![0.5_f64; n5 * n5];
    for i in 0..n5 {
        eq[i * n5 + i] = 1.0;
    }
    let mut ws5 = ProbabilityWorkspace::new(n5 - 1);
    assert_eq!(ws5.evaluate(&eq, n5, 0.0).unwrap(), 1.0);

    // Exact SCPC integral references for 5x5 equicorrelated covariance (rho = 0.5).
    // Analytical DPR1 eigenvalues: mu0=2.5 (w0=0.4), mu1..3=0.5 (w=0), s=0.6.
    let prob_refs = [
        (1.0_f64, 9.98925194663871224e-02_f64),
        (1e-6_f64, 1.0_f64 - 1.28265778655618448e-06_f64),
    ];
    for (cv, expected) in prob_refs {
        let got = ws5.evaluate(&eq, n5, cv).unwrap();
        let err = if cv < 1.0 {
            ((1.0 - got) - (1.0 - expected)).abs() / (1.0 - expected)
        } else {
            (got - expected).abs() / expected
        };
        assert!(
            err < 1e-8,
            "cv={cv}: got={got}, expected={expected}, rel={err}"
        );
    }
    // Check eigenvalue ratios directly across 24 orders of magnitude (1e-12 to 1e12).
    let mut eigensolver = crate::linalg::EigenWorkspace::new(n5);
    let mut r_buf = Vec::new();
    let mut l_buf = Vec::new();
    for cv in [1e-12_f64, 1e-10, 1e-8, 1.0, 1e8, 1e10, 1e12] {
        assert!(eigensolver
            .rejection_eigen_ratios(&eq, n5, cv, &mut r_buf, &mut l_buf)
            .unwrap()
            .is_none());
        assert_eq!(r_buf.len(), 4);
        let a2 = cv * cv;
        let (s_eq, w_eq, mu_eq) = (0.6_f64, 0.4_f64, 2.5_f64);
        // Solve quadratic for lambda_pos and eta analytically:
        let b_lam = a2 * mu_eq - 1.0;
        let c_lam = -s_eq * a2 * mu_eq;
        let lam_ref = if a2 * mu_eq < 1e-6 {
            1.0 - w_eq * a2 * mu_eq
        } else if a2 * mu_eq > 1e10 {
            s_eq + w_eq / (a2 * mu_eq)
        } else {
            let disc = (b_lam * b_lam - 4.0 * c_lam).sqrt();
            if b_lam <= 0.0 {
                0.5 * (-b_lam + disc)
            } else {
                (-2.0 * c_lam) / (b_lam + disc)
            }
        };
        let eta_ref = if a2 < 1e-8 {
            s_eq * mu_eq
        } else if a2 > 1e10 {
            mu_eq - w_eq / a2
        } else {
            let b_e = 1.0 - a2 * mu_eq;
            let c_e = -s_eq * mu_eq;
            let d_e = (b_e * b_e - 4.0 * a2 * c_e).sqrt();
            if b_e <= 0.0 {
                (a2 * mu_eq - 1.0 + d_e) / (2.0 * a2)
            } else {
                (-2.0 * c_e) / (b_e + d_e)
            }
        };
        let log_r0_ref = 2.0 * cv.ln() + eta_ref.ln() - lam_ref.ln();
        let log_rd_ref = 2.0 * cv.ln() + (0.5_f64).ln() - lam_ref.ln();
        let mut l_sorted = l_buf.clone();
        l_sorted.sort_by(|a, b| a.total_cmp(b));
        for i in 0..3 {
            assert!(
                (l_sorted[i] - log_rd_ref).abs() < 1e-12,
                "cv={cv} rd mismatch"
            );
        }
        assert!(
            (l_sorted[3] - log_r0_ref).abs() < 1e-12,
            "cv={cv} r0 mismatch"
        );
    }
    // Very large tail cv=1e8 still integrates cleanly (> 0 or tail bounded).
    assert!(ws5.evaluate(&eq, n5, 1e8).unwrap() < 1e-12);

    // 2. Correlated 10x10 SPD matrix at cv == 0.0 -> exact 1.0.
    let n10 = 10usize;
    let mut spd10 = vec![0.0_f64; n10 * n10];
    for i in 0..n10 {
        for j in 0..n10 {
            spd10[i * n10 + j] = (-0.25 * ((i as f64) - (j as f64)).abs()).exp();
        }
    }
    let mut ws10 = ProbabilityWorkspace::new(n10 - 1);
    assert_eq!(ws10.evaluate(&spd10, n10, 0.0).unwrap(), 1.0);

    // 3. Extreme diagonal scale disparity diag(1e200, 1e-120, 1e-120) at cv = 1e160 -> 1 - 1/sqrt(2).
    let mut ws3 = ProbabilityWorkspace::new(2);
    let p_diag = ws3
        .evaluate(
            &[1e200, 0.0, 0.0, 0.0, 1e-120, 0.0, 0.0, 0.0, 1e-120],
            3,
            1e160,
        )
        .unwrap();
    let expected_diag = 1.0 - std::f64::consts::FRAC_1_SQRT_2;
    assert!((p_diag - expected_diag).abs() < 1e-12, "p_diag={p_diag}");

    // 4. Rank-1 extreme scale disparity v = (2^400, 2^-400, 2^-400).
    let hi = 2.0_f64.powi(400);
    let lo = 2.0_f64.powi(-400);
    let v = [hi, lo, lo];
    let mut r1 = [0.0_f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            r1[i * 3 + j] = v[i] * v[j];
        }
    }
    let thresh = (hi / lo) / 2.0_f64.sqrt();
    assert_eq!(ws3.evaluate(&r1, 3, 0.5 * thresh).unwrap(), 1.0);
    assert_eq!(ws3.evaluate(&r1, 3, 2.0 * thresh).unwrap(), 0.0);

    // 5. build_scale_grid must not panic when cmax / c0 > f64::MAX.
    let grid = build_scale_grid(1e-200, 1e200);
    assert!(grid.len() > 2 && grid.iter().all(|&v| v.is_finite()));
}
