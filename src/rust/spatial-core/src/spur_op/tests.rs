use super::halflife::{
    accepts_halflife_null, add_halflife_likelihood_terms, half_life_covariances, half_life_mixture,
    half_life_search_grid, half_power, halflife_cancellation_bound, halflife_tile_size,
    invert_half_life, log_sum_exp, quadratic_form_from_lower, MAX_HALFLIFE_TILE_BYTES,
};
use super::*;
use crate::ErrorKind;

fn load_bench_200_spur() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/bench_200.csv"
    ))
    .expect("fixture csv");
    let mut y = Vec::new();
    let mut x_in = Vec::new();
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
        x_in.extend_from_slice(&[1.0, cols[2], cols[3]]);
        // (s_1, s_2)
        coords.extend_from_slice(&[cols[6], cols[7]]);
    }
    (y, x_in, coords)
}

// Build an explicit transformation matrix as an independent reference for
// the streamed NN/ISO operations, including their parallel path.
fn dense_neighbor_reference(req: &SpurTransformRequest<'_>) -> Vec<f64> {
    let n = req.n;
    let distances = distance_matrix(req.coords.unwrap(), n, req.d, req.metric).unwrap();
    let mut contrast = vec![0.0; n * n];
    for (i, (dst, src)) in contrast
        .chunks_exact_mut(n)
        .zip(distances.chunks_exact(n))
        .enumerate()
    {
        let radius = match req.kind {
            SpurTransformKind::Nn => src
                .iter()
                .copied()
                .filter(|&d| d > 0.0)
                .fold(f64::INFINITY, f64::min),
            SpurTransformKind::Iso { radius } => radius,
            _ => unreachable!(),
        };
        let selected = |j: usize, d: f64| match req.kind {
            SpurTransformKind::Nn => d > 0.0 && d == radius,
            _ => i != j && d < radius,
        };
        let count = src
            .iter()
            .enumerate()
            .filter(|&(j, &d)| selected(j, d))
            .count();
        if count > 0 {
            for (j, (cell, &d)) in dst.iter_mut().zip(src).enumerate() {
                if selected(j, d) {
                    *cell = -1.0 / count as f64;
                }
            }
            dst[i] = 1.0;
        }
    }
    let mut out = vec![0.0; req.y.len()];
    linalg::multiply(&contrast, req.y, &mut out, n, n, req.n_y, false);
    out
}

#[test]
fn streamed_neighbor_transforms_match_dense_reference() {
    let coords = [0.0, 0.0, 1.0, 2.0, 6.0];
    let y = [1.0, 3.0, 5.0, 7.0, 9.0];
    let nn = SpurTransformKind::Nn;
    let iso = SpurTransformKind::Iso { radius: 1.0 };
    let base = SpurTransformRequest {
        coords: Some(&coords),
        n: coords.len(),
        d: 1,
        metric: Metric::Euclidean,
        kind: nn,
        cluster_id: None,
        y: &y,
        n_y: 1,
        nthreads: Some(1),
    };
    let compare = |req: SpurTransformRequest<'_>| {
        let expected = dense_neighbor_reference(&req);
        let values = run_spur_transform(&req).unwrap();
        for (&actual, expected) in values.iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-12);
        }
        values
    };
    for (kind, expected) in [
        (nn, [-4.0, -2.0, 4.0 / 3.0, 2.0, 2.0]),
        (iso, [-2.0, 2.0, 0.0, 0.0, 0.0]),
    ] {
        let values = compare(SpurTransformRequest { kind, ..base });
        for (actual, expected) in values.into_iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-12);
        }
    }
    let err = run_spur_transform(&SpurTransformRequest {
        coords: Some(&[5.0; 3]),
        n: 3,
        y: &[1.0, 3.0, 5.0],
        ..base
    })
    .unwrap_err();
    assert_eq!(err.kind(), ErrorKind::InvalidInput);
    assert_eq!(
        compare(SpurTransformRequest {
            coords: Some(&[0.0, 1e-200, 1e160]),
            n: 3,
            y: &[1.0, 3.0, 5.0],
            ..base
        }),
        vec![-2.0, 2.0, 3.0]
    );
    let coords_many: Vec<f64> = (0..80).flat_map(|i| [i as f64, (i % 5) as f64]).collect();
    let y_many: Vec<f64> = (0..160).map(|i| (i as f64).sin()).collect();
    for kind in [nn, SpurTransformKind::Iso { radius: 0.01 }] {
        compare(SpurTransformRequest {
            coords: Some(&coords_many),
            n: 80,
            d: 2,
            metric: Metric::GreatCircleFractionLatLon,
            kind,
            y: &y_many,
            n_y: 2,
            nthreads: Some(4),
            ..base
        });
    }
}

#[test]
fn cluster_transform_allows_one_observation() {
    let result = run_spur_transform(&SpurTransformRequest {
        coords: None,
        n: 1,
        d: 0,
        metric: Metric::Euclidean,
        kind: SpurTransformKind::Cluster,
        cluster_id: Some(&[0]),
        y: &[2.0, -3.0],
        n_y: 2,
        nthreads: Some(1),
    })
    .unwrap();
    assert_eq!(result, [0.0; 2]);
}

#[test]
fn spur_basis_orientation_uses_the_first_nearly_maximal_entry() {
    for (offset, first_is_positive) in [(5e-13, true), (2e-12, false)] {
        let a = -1.0;
        let b = 1.0 + offset;
        let covariance = [a * a, a * b, a * b, b * b];
        let geom = SpurGeometry::from_coords(&[0.0, 1.0], 2, 1, Metric::Euclidean).unwrap();
        let (basis, info) = select_leading_basis(&covariance, &geom, 1).unwrap();
        assert!(!info.tied);
        assert_eq!(basis[0].is_sign_positive(), first_is_positive);
        assert_eq!(basis[1].is_sign_negative(), first_is_positive);
    }
}

#[test]
fn halflife_worker_scratch_and_likelihood_terms_are_stable() {
    let n_alternatives = 165;
    for q in [4, 64, 500, 1000, 1200] {
        let pairs = q * (q + 1) / 2;
        let (cols, packed) = halflife_tile_size(q, pairs, n_alternatives, 10_000);
        assert!(cols <= TILE_COLS);
        assert_eq!(packed, q <= 64);
        if packed {
            assert!(pairs * n_alternatives * size_of::<f64>() <= MAX_HALFLIFE_TILE_BYTES);
        }
        let tile_doubles = cols
            * (2 * q
                + if packed {
                    pairs + n_alternatives + 2
                } else {
                    0
                });
        assert!(tile_doubles * size_of::<f64>() <= MAX_HALFLIFE_TILE_BYTES);
    }
    let lower = [1e-7, 0.0, -1.0, 1.0];
    let draws = [1.0, 2.0, 1.0, 3.0];
    let form = quadratic_form_from_lower(&lower, &draws, 2, 0);
    assert!((form / 1e-14 - 1.0).abs() < 4.0 * f64::EPSILON);
    let log_total = log_sum_exp(&[1000.0, 1000.0]);
    assert!((log_total - 1000.0 - 2.0_f64.ln()).abs() < 1e-13);
    assert!(halflife_cancellation_bound(&[1e308, 1e308], 1).is_finite());
}

#[test]
fn halflife_confidence_set_excludes_exact_monte_carlo_boundaries() {
    for (level, boundary) in [(0.90, 9), (0.95, 4)] {
        assert!(!accepts_halflife_null(boundary - 1, 99, level));
        assert!(!accepts_halflife_null(boundary, 99, level));
        assert!(accepts_halflife_null(boundary + 1, 99, level));
    }
}

fn make_deterministic_emat(q: usize, nrep: usize) -> Vec<f64> {
    // Box-Muller from deterministic LCG for self-contained Rust tests.
    let mut state: u64 = 0x1234_5678_9abc_def0;
    let mut next_u = || -> f64 {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let bits = (state >> 11) as f64 * (1.0 / (1_u64 << 53) as f64);
        bits.clamp(1e-12, 1.0 - 1e-12)
    };
    let mut out = vec![0.0_f64; q * nrep];
    for chunk in out.chunks_mut(2) {
        let u1 = next_u();
        let u2 = next_u();
        let r = (-2.0 * u1.ln()).sqrt();
        let th = 2.0 * std::f64::consts::PI * u2;
        chunk[0] = r * th.cos();
        if chunk.len() > 1 {
            chunk[1] = r * th.sin();
        }
    }
    out
}

#[test]
fn spur_decay_brackets_large_roots_and_checks_the_colocation_floor() {
    let coords: Vec<f64> = (0..10).map(|i| i as f64 * 1e-8).chain([1.0]).collect();
    let mut geom = SpurGeometry::from_coords(&coords, coords.len(), 1, Metric::Euclidean).unwrap();
    let targets = [0.03, 0.001, 0.0001];
    let roots = geom.get_cbar_batch(&targets).unwrap();
    for (&target, root) in targets.iter().zip(roots) {
        assert!(root > 10_240.0);
        assert!((deterministic_mean_exp(&geom.lower_tri, root) - target).abs() < 1e-12);
        assert_eq!(geom.get_cbar(target), Ok(root));
    }
    let mut colocated =
        SpurGeometry::from_coords(&[0.0, 0.0, 1.0], 3, 1, Metric::Euclidean).unwrap();
    let err = colocated.get_cbar(1.0 / 3.0).unwrap_err();
    assert_eq!(err.kind(), ErrorKind::CalibrationFailed);
    assert!(err.message().contains("co-located-pair fraction"));
    assert!(colocated.get_cbar_batch(&[0.5, 0.1]).is_err());
}

#[test]
fn residual_brownian_basis_preserves_orientation_at_paired_maxima() {
    let mut coords: Vec<f64> = (0..60).map(|i| i as f64 * (1.0 / 59.0)).collect();
    coords[59] = 1.0;
    let design: Vec<f64> = coords.iter().flat_map(|&s| [1.0, s - 0.5]).collect();
    let q = 15;
    let design = prepare_resid_design(&design, coords.len(), 2, q).unwrap();
    let geom = SpurGeometry::from_coords(&coords, coords.len(), 1, Metric::Euclidean).unwrap();
    let basis = prepare_resid_basis(&geom, design, q).unwrap();

    // The second direction has opposite signs at its mirrored maxima.
    // The sign convention selects the earlier site when their magnitudes agree.
    let earlier = basis.basis_r[19 * q + 1];
    let later = basis.basis_r[40 * q + 1];
    assert!(earlier > 0.18 && (earlier + later).abs() < 1e-14);
}

#[test]
fn resolved_eigenvalue_gap_survives_projection_to_a_small_scale() {
    let coords = [0.0, 1e-12, 2e-12, 4e-12, 6e-12, 9e-12]
        .into_iter()
        .chain([0.0, 1e-12, 2e-12, 4e-12, 5e-12, 7e-12].map(|d| 1.0 + d))
        .collect::<Vec<_>>();
    let design: Vec<_> = (0..12)
        .flat_map(|i| if i < 6 { [1.0, 0.0] } else { [0.0, 1.0] })
        .collect();
    let geom = SpurGeometry::from_coords(&coords, 12, 1, Metric::Euclidean).unwrap();
    let design = prepare_resid_design(&design, 12, 2, 2).unwrap();
    let basis = prepare_resid_basis(&geom, design, 2).unwrap();
    assert!(!basis.info.tied);
}

#[test]
fn spur_cbar_and_diagnostics_are_deterministic_and_thread_invariant() {
    let (y, x_in, coords) = load_bench_200_spur();
    let n = y.len();
    let q = 8usize;
    let nrep = 1000usize;
    let emat = make_deterministic_emat(2 * q, nrep);

    let mut geom = SpurGeometry::from_coords(&coords, n, 2, Metric::Euclidean).unwrap();
    let c099 = geom.get_cbar(0.999).unwrap();
    let c099_cached = geom.get_cbar(0.999).unwrap();
    assert_eq!(c099.to_bits(), c099_cached.to_bits());
    assert!((c099 - 0.002_662_372_394_202_289_4).abs() < 1e-12);
    let mean_rho = deterministic_mean_exp(&geom.lower_tri, c099);
    assert!((mean_rho - 0.999).abs() < 1e-3, "mean_rho={mean_rho}");

    // Pipeline vs individual tests match to machine precision and are thread invariant.
    let pipe_1 = run_spur_pipeline(&SpurPipelineRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        x_in: &x_in,
        k_x: 3,
        emat: &emat,
        q,
        nrep,
        transform_vars: &x_in,
        n_transform_vars: 3,
        nthreads: Some(1),
    })
    .unwrap();
    let pipe_4 = run_spur_pipeline(&SpurPipelineRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        x_in: &x_in,
        k_x: 3,
        emat: &emat,
        q,
        nrep,
        transform_vars: &x_in,
        n_transform_vars: 3,
        nthreads: Some(4),
    })
    .unwrap();
    let i0_ind = run_spur_i0(&SpurTestRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        n_y: 1,
        x_in: None,
        k_x: 0,
        emat: &emat,
        q,
        nrep,
        nthreads: Some(2),
    })
    .unwrap();
    let i0r_ind = run_spur_i0_resid(&SpurTestRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        n_y: 1,
        x_in: Some(&x_in),
        k_x: 3,
        emat: &emat,
        q,
        nrep,
        nthreads: Some(2),
    })
    .unwrap();

    // Finite Monte Carlo p-values are attainable ranks.
    for p in [
        &pipe_1.i0.pvalue,
        &pipe_1.i1.pvalue,
        &pipe_1.i0resid.pvalue,
        &pipe_1.i1resid.pvalue,
    ] {
        let rank = p[0] * (nrep + 1) as f64;
        assert!(rank >= 1.0 && (rank - rank.round()).abs() < 1e-12);
    }

    assert!((pipe_1.i0.lr[0] - i0_ind.lr[0]).abs() < 1e-12);
    assert!((pipe_1.i0.pvalue[0] - i0_ind.pvalue[0]).abs() < 1e-12);
    assert!((pipe_1.i0.ha_param - i0_ind.ha_param).abs() < 1e-12);
    assert!((pipe_1.i0resid.lr[0] - i0r_ind.lr[0]).abs() < 1e-12);
    assert!((pipe_1.i0resid.pvalue[0] - i0r_ind.pvalue[0]).abs() < 1e-12);
    let transform_ind = run_spur_transform(&SpurTransformRequest {
        coords: Some(&coords),
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &x_in,
        n_y: 3,
        kind: SpurTransformKind::LbmGls,
        cluster_id: None,
        nthreads: Some(1),
    })
    .unwrap();
    assert_eq!(pipe_1.transformed_vars, transform_ind);

    // Check I(0) grid maxima invariant
    assert_eq!(pipe_1.i0.rho_grid.len(), 31);
    let max_p = pipe_1.i0.pvalue_mat.iter().copied().fold(0.0_f64, f64::max);
    assert!((pipe_1.i0.pvalue[0] - max_p).abs() < 1e-15);

    // Thread count invariance
    assert!((pipe_1.i0.lr[0] - pipe_4.i0.lr[0]).abs() < 1e-13);
    assert!((pipe_1.i1.lr[0] - pipe_4.i1.lr[0]).abs() < 1e-13);
    assert!((pipe_1.i0resid.lr[0] - pipe_4.i0resid.lr[0]).abs() < 1e-13);
    assert!((pipe_1.i1resid.lr[0] - pipe_4.i1resid.lr[0]).abs() < 1e-13);
    for (&a, &b) in pipe_1.transformed_vars.iter().zip(&pipe_4.transformed_vars) {
        assert!((a - b).abs() < 1e-13);
    }
}

#[test]
fn spur_halflife_and_lbmgls_invariance_and_colocation() {
    let (y, _, coords) = load_bench_200_spur();
    let n = y.len();
    let q = 8usize;
    let nrep = 300usize;
    let emat = make_deterministic_emat(q, nrep);

    let hl_1 = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        emat: &emat,
        q,
        nrep,
        level: 0.95,
        nthreads: Some(1),
    })
    .unwrap();
    let hl_2 = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        emat: &emat,
        q,
        nrep,
        level: 0.95,
        nthreads: Some(2),
    })
    .unwrap();
    assert!(hl_1.pv_mat.len() >= 129);
    assert_eq!(hl_1.pv_mat.len(), hl_1.hl_grid_ho.len());
    assert!(hl_1.ci_lower > 0.0);
    assert_eq!(hl_1.ci_upper, f64::INFINITY);
    assert_eq!(hl_1.components, hl_2.components);
    for (&a, &b) in hl_1.pv_mat.iter().zip(&hl_2.pv_mat) {
        assert!((a - b).abs() < 1e-14);
    }

    // Affine invariance: y -> -2.5 * y + 17 preserves angular likelihood.
    let y_aff: Vec<f64> = y.iter().map(|&v| -2.5 * v + 17.0).collect();
    let hl_aff = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y_aff,
        emat: &emat,
        q,
        nrep,
        level: 0.95,
        nthreads: Some(2),
    })
    .unwrap();
    for (&a, &b) in hl_1.pv_mat.iter().zip(&hl_aff.pv_mat) {
        assert!((a - b).abs() < 1e-12, "p_orig={a}, p_aff={b}");
    }

    // Co-located LBM-GLS rank < n - 1 and zero constant vector action.
    let coords_col = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    // Pass the 5 x 2 row-major matrix [y, ones].
    let y_mat = [1.0, 1.0, 2.0, 1.0, -1.0, 1.0, 0.5, 1.0, 3.0, 1.0];
    let req = SpurTransformRequest {
        coords: Some(&coords_col),
        n: 5,
        d: 2,
        metric: Metric::Euclidean,
        kind: SpurTransformKind::LbmGls,
        cluster_id: None,
        y: &y_mat,
        n_y: 2,
        nthreads: Some(1),
    };
    let lbm = run_spur_transform(&req).unwrap();

    // Reference H = V diag(lambda^-1/2) V' for the centered covariance -P D P / 2.
    let geometry = SpurGeometry::from_coords(&coords_col, 5, 2, Metric::Euclidean).unwrap();
    let mut sigma: Vec<f64> = geometry.norm_distmat.iter().map(|d| -0.5 * d).collect();
    double_demean(&mut sigma, 5);
    let pairs = basis::dense_eigenpairs(&mut sigma, 5, 5).unwrap();
    let mut h = [0.0; 25];
    let mut rank = 0;
    for (vector, &val) in pairs.vectors.chunks_exact(5).zip(&pairs.values) {
        if val > 1e-10 {
            rank += 1;
            for i in 0..5 {
                for j in 0..5 {
                    h[i * 5 + j] += vector[i] * vector[j] / val.sqrt();
                }
            }
        }
    }
    assert_eq!(rank, 3);
    let mut expected = [0.0; 10];
    linalg::multiply(&h, &y_mat, &mut expected, 5, 5, 2, false);
    center_columns(&mut expected, 5, 2);
    for (&actual, expected) in lbm.iter().zip(expected) {
        assert!((actual - expected).abs() < 1e-12);
    }
    for row in lbm.chunks_exact(2) {
        assert!(row[1].abs() < 1e-12);
    }
    assert!((lbm[0] - lbm[2]).abs() < 1e-12);
    assert!(lbm.chunks_exact(2).map(|r| r[0]).sum::<f64>().abs() < 1e-12);

    let iso = run_spur_transform(&SpurTransformRequest {
        kind: SpurTransformKind::Iso { radius: 0.5 },
        ..req
    })
    .unwrap();
    let delta = y_mat[0] - y_mat[2];
    assert_eq!(
        iso,
        vec![delta, 0.0, -delta, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
}

#[test]
fn spur_validates_residual_and_halflife_dimensions_and_scales_invariantly() {
    let n = 12usize;
    let k_x = 3usize;
    let coords: Vec<f64> = (0..n)
        .flat_map(|i| {
            let fi = i as f64;
            [(fi * 0.7).sin() + 0.1 * fi, (fi * 1.1).cos()]
        })
        .collect();
    let y: Vec<f64> = (0..n).map(|i| ((i as f64) * 0.53).sin() + 0.2).collect();
    let mut x_in = Vec::with_capacity(n * k_x);
    for i in 0..n {
        let fi = i as f64;
        x_in.extend_from_slice(&[1.0, (fi * 0.37).cos(), (fi * 0.89).sin() + 0.1 * fi]);
    }

    // Finding 3a: q > n - k_x is rejected cleanly (q=10, 11 for n=12, k_x=3 -> dof=9).
    for q_bad in [n - k_x + 1, n - 1] {
        let emat_bad = make_deterministic_emat(2 * q_bad, 50);
        let err0 = run_spur_i0_resid(&SpurTestRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            n_y: 1,
            x_in: Some(&x_in),
            k_x,
            emat: &emat_bad,
            q: q_bad,
            nrep: 50,
            nthreads: Some(1),
        })
        .unwrap_err();
        assert_eq!(err0.kind(), ErrorKind::InvalidInput);
        let err1 = run_spur_i1_resid(&SpurTestRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &y,
            n_y: 1,
            x_in: Some(&x_in),
            k_x,
            emat: &emat_bad,
            q: q_bad,
            nrep: 50,
            nthreads: Some(1),
        })
        .unwrap_err();
        assert_eq!(err1.kind(), ErrorKind::InvalidInput);
    }

    // Legal boundary q == n - k_x succeeds cleanly.
    let q_legal = n - k_x;
    let emat_legal = make_deterministic_emat(2 * q_legal, 60);
    let ok_res = run_spur_i0_resid(&SpurTestRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        n_y: 1,
        x_in: Some(&x_in),
        k_x,
        emat: &emat_legal,
        q: q_legal,
        nrep: 60,
        nthreads: Some(1),
    })
    .unwrap();
    assert!(ok_res.lr[0].is_finite() && ok_res.pvalue[0].is_finite());

    // Co-located points reducing covariance rank below q are rejected explicitly.
    let mut coords_col = Vec::with_capacity(n * 2);
    for i in 0..n {
        let site = (i % 6) as f64;
        coords_col.extend_from_slice(&[site, (site * 1.3).sin()]);
    }
    let emat_rank = make_deterministic_emat(2 * 7, 40);
    let err_rank = run_spur_i0_resid(&SpurTestRequest {
        coords: &coords_col,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        n_y: 1,
        x_in: Some(&x_in),
        k_x,
        emat: &emat_rank,
        q: 7,
        nrep: 40,
        nthreads: Some(1),
    })
    .unwrap_err();
    assert_eq!(err_rank.kind(), ErrorKind::InvalidInput);

    // Finding 3b: Half-life rejects q=1, constant y, and zero projection, while accepting n=3, q=2.
    let emat_1 = make_deterministic_emat(1, 50);
    let err_q1 = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        emat: &emat_1,
        q: 1,
        nrep: 50,
        level: 0.95,
        nthreads: Some(1),
    })
    .unwrap_err();
    assert_eq!(err_q1.kind(), ErrorKind::InvalidInput);
    let q1_test = SpurTestRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &y,
        n_y: 1,
        x_in: Some(&x_in),
        k_x,
        emat: &emat_1,
        q: 1,
        nrep: 50,
        nthreads: Some(1),
    };
    assert_eq!(
        run_spur_i0(&q1_test).unwrap_err().kind(),
        ErrorKind::InvalidInput
    );
    assert_eq!(
        run_spur_i1(&q1_test).unwrap_err().kind(),
        ErrorKind::InvalidInput
    );
    assert_eq!(
        run_spur_i0_resid(&q1_test).unwrap_err().kind(),
        ErrorKind::InvalidInput
    );
    assert_eq!(
        run_spur_i1_resid(&q1_test).unwrap_err().kind(),
        ErrorKind::InvalidInput
    );

    let emat_2 = make_deterministic_emat(2, 50);
    // All centered covariance shapes on an equilateral triangle are scalar
    // multiples of the identity, despite retaining two positive directions.
    let equilateral = [0.0, 0.0, 1.0, 0.0, 0.5, 3.0_f64.sqrt() / 2.0];
    let err_proportional = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &equilateral,
        n: 3,
        d: 2,
        metric: Metric::Euclidean,
        y: &[0.0, 1.0, -0.5],
        emat: &emat_2,
        q: 2,
        nrep: 50,
        level: 0.95,
        nthreads: Some(1),
    })
    .unwrap_err();
    assert_eq!(err_proportional.kind(), ErrorKind::DegenerateVariance);
    let equilateral_y = [0.0, 1.0, -0.5];
    let intercept = [1.0, 1.0, 1.0];
    let diagnostic_draws = make_deterministic_emat(4, 50);
    let equilateral_test = SpurTestRequest {
        coords: &equilateral,
        n: 3,
        d: 2,
        metric: Metric::Euclidean,
        y: &equilateral_y,
        n_y: 1,
        x_in: Some(&intercept),
        k_x: 1,
        emat: &diagnostic_draws,
        q: 2,
        nrep: 50,
        nthreads: Some(1),
    };
    assert_eq!(
        run_spur_i0(&equilateral_test).unwrap_err().kind(),
        ErrorKind::DegenerateVariance
    );
    assert_eq!(
        run_spur_i1(&equilateral_test).unwrap_err().kind(),
        ErrorKind::DegenerateVariance
    );
    assert_eq!(
        run_spur_i0_resid(&equilateral_test).unwrap_err().kind(),
        ErrorKind::DegenerateVariance
    );
    assert_eq!(
        run_spur_i1_resid(&equilateral_test).unwrap_err().kind(),
        ErrorKind::DegenerateVariance
    );
    let const_y = vec![5.0; n];
    let err_const = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &const_y,
        emat: &emat_2,
        q: 2,
        nrep: 50,
        level: 0.95,
        nthreads: Some(1),
    })
    .unwrap_err();
    assert_eq!(err_const.kind(), ErrorKind::DegenerateVariance);

    // 3-site legal half-life check (n=3, q=2).
    let coords_3 = [0.0, 0.0, 1.0, 0.0, 0.3, 0.9];
    let y_3 = [0.0, 1.0, -0.5];
    let hl_3 = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords_3,
        n: 3,
        d: 2,
        metric: Metric::Euclidean,
        y: &y_3,
        emat: &emat_2,
        q: 2,
        nrep: 50,
        level: 0.95,
        nthreads: Some(1),
    })
    .unwrap();
    assert!(hl_3.ci_lower.is_finite() && !hl_3.ci_upper.is_nan());

    // Finding 5: Column scaling of X_in (1e-12 and 1e12) leaves projection and inference invariant,
    // while true rank deficiency is still rejected.
    let (_, _, coords_200) = load_bench_200_spur();
    let n200 = 100usize;
    let coords_sub = &coords_200[..n200 * 2];
    let y_sub: Vec<f64> = (0..n200).map(|i| ((i as f64) * 0.31).sin()).collect();
    let mut x_orig = Vec::with_capacity(n200 * 2);
    let mut x_tiny = Vec::with_capacity(n200 * 2);
    let mut x_huge = Vec::with_capacity(n200 * 2);
    for i in 0..n200 {
        let v = ((i as f64) * 0.77).cos() + 0.1 * ((i as f64) * 0.13).sin();
        x_orig.extend_from_slice(&[1.0, v]);
        x_tiny.extend_from_slice(&[1.0, v * 1e-12]);
        x_huge.extend_from_slice(&[1.0, v * 1e12]);
    }
    let emat_s = make_deterministic_emat(2 * 8, 200);
    let run_scaled = |xm: &[f64]| {
        run_spur_i0_resid(&SpurTestRequest {
            coords: coords_sub,
            n: n200,
            d: 2,
            metric: Metric::Euclidean,
            y: &y_sub,
            n_y: 1,
            x_in: Some(xm),
            k_x: 2,
            emat: &emat_s,
            q: 8,
            nrep: 200,
            nthreads: Some(1),
        })
        .unwrap()
    };
    let r_orig = run_scaled(&x_orig);
    let r_tiny = run_scaled(&x_tiny);
    let r_huge = run_scaled(&x_huge);
    assert!((r_orig.lr[0] - r_tiny.lr[0]).abs() < 1e-11);
    assert!((r_orig.lr[0] - r_huge.lr[0]).abs() < 1e-11);
    assert!((r_orig.ha_param - r_tiny.ha_param).abs() < 1e-12);
    assert!((r_orig.pvalue[0] - r_tiny.pvalue[0]).abs() < 1e-12);

    // Collinear column still rejected despite scaling.
    let mut x_dep = Vec::with_capacity(n200 * 3);
    for i in 0..n200 {
        let v = x_orig[i * 2 + 1];
        x_dep.extend_from_slice(&[1.0, v * 1e-12, v * 3e-12]);
    }
    let err_dep = orthonormal_design_basis(&x_dep, n200, 3).unwrap_err();
    assert_eq!(err_dep.kind(), ErrorKind::InvalidInput);

    // Division by the column maximum remains valid even when its reciprocal
    // would overflow, and response location does not determine degeneracy.
    let tiny = [1e-320, 0.0, 0.0, 1e-320, -1e-320, -1e-320];
    assert!(orthonormal_design_basis(&tiny, 3, 2).is_ok());
    let basis = [1.0 / 2.0_f64.sqrt(), -1.0 / 2.0_f64.sqrt(), 0.0];
    let p = prepare_y_projections(&[1e13, 1e13 + 1.0, 1e13 - 1.0], 3, 1, &basis, 1).unwrap();
    assert!(p[0].is_finite() && p[0] != 0.0);

    // Residual degeneracy is invariant to an additive offset when the design
    // spans constants, even with the constant represented by the column sum.
    let x: Vec<f64> = (0..3).flat_map(|i| [i as f64, (1 - i) as f64]).collect();
    let design = prepare_resid_design(&x, 3, 2, 1).unwrap();
    let residual_projection = |y| prepare_residual_y_projections(y, 3, 1, &basis, &design, 1);
    let unshifted = residual_projection(&[-1.0, 1.0, 0.0]).unwrap();
    assert_eq!(
        residual_projection(&[1e13 - 1.0, 1e13 + 1.0, 1e13]).unwrap(),
        unshifted
    );
    assert_eq!(
        residual_projection(&[-f64::MAX, f64::MAX, 0.0]).unwrap(),
        unshifted
    );
}

#[test]
fn optimized_omega_grid_quantiles_halflife_and_tiny_coords_match_reference() {
    // 1. projected_omega_grid parity for n < 256 and n > 256 (serial & parallel)
    let (_, _, coords_200) = load_bench_200_spur();
    let geom_200 = SpurGeometry::from_coords(&coords_200, 200, 2, Metric::Euclidean).unwrap();
    let (basis_200, _, mr_200, _) =
        prepare_lbm_basis(&geom_200.sigma_lbm_dm(), &geom_200, 8).unwrap();
    let c_vals = [0.25, 1.5, 8.0, 42.0, 120.0];
    let grid_serial = run_with_threads(Some(1), || {
        geom_200.projected_omega_grid(&c_vals, &mr_200, 8)
    })
    .unwrap();
    let grid_par = run_with_threads(Some(4), || {
        geom_200.projected_omega_grid(&c_vals, &mr_200, 8)
    })
    .unwrap();
    for (idx, &c) in c_vals.iter().enumerate() {
        let ref_om = geom_200.projected_omega_from_mr(c, &mr_200, 8);
        for k in 0..64 {
            assert!(
                (grid_serial[idx][k] - ref_om[k]).abs() < 1e-12,
                "n=200 serial diff at c={c}: {}",
                (grid_serial[idx][k] - ref_om[k]).abs()
            );
            assert_eq!(grid_serial[idx][k].to_bits(), grid_par[idx][k].to_bits());
        }
    }
    let _ = basis_200;

    let n_big = 513usize;
    let coords_big: Vec<f64> = (0..n_big)
        .flat_map(|i| {
            let f = i as f64;
            [(f * 0.37).sin() + 0.05 * f, (f * 0.83).cos()]
        })
        .collect();
    let geom_big = SpurGeometry::from_coords(&coords_big, n_big, 2, Metric::Euclidean).unwrap();
    let (_, _, mr_big, _) = prepare_lbm_basis(&geom_big.sigma_lbm_dm(), &geom_big, 6).unwrap();
    let big_serial = run_with_threads(Some(1), || {
        geom_big.projected_omega_grid(&c_vals, &mr_big, 6)
    })
    .unwrap();
    let big_par = run_with_threads(Some(4), || {
        geom_big.projected_omega_grid(&c_vals, &mr_big, 6)
    })
    .unwrap();
    for (idx, &c) in c_vals.iter().enumerate() {
        let ref_om = geom_big.projected_omega_from_mr(c, &mr_big, 6);
        for k in 0..36 {
            assert!((big_serial[idx][k] - ref_om[k]).abs() < 1e-12);
            assert_eq!(big_serial[idx][k].to_bits(), big_par[idx][k].to_bits());
        }
    }

    // Rank-test critical values use the corresponding null order statistics.
    let mut draws: Vec<_> = (1..=99).rev().map(f64::from).collect();
    assert_eq!(monte_carlo_cutoffs(&mut draws), vec![99.0, 95.0, 90.0]);
    assert_eq!(monte_carlo_pvalue(&draws, 99.0), 0.02);
    assert_eq!(monte_carlo_pvalue(&draws, 100.0), 0.01);
    assert_eq!(monte_carlo_pvalue(&[5.0; 99], 5.0), 1.0);
    assert_eq!(monte_carlo_cutoffs(&mut [1.0; 8]), vec![f64::INFINITY; 3]);

    // 3. Compare both single- and multiple-tile, fast and log-scale
    // half-life paths to a factor-based calculation.
    let n_hl = 64usize;
    let coords_hl = &coords_200[..n_hl * 2];
    let y_hl: Vec<f64> = (0..n_hl).map(|i| ((i as f64) * 0.41).sin() + 0.1).collect();
    for (q_hl, nrep_hl) in [(4, TILE_COLS + 3), (8, 80), (15, 80), (60, 7), (61, 7)] {
        let emat_hl = make_deterministic_emat(q_hl, nrep_hl);
        let out_fast = run_spur_halflife(&SpurHalfLifeRequest {
            coords: coords_hl,
            n: n_hl,
            d: 2,
            metric: Metric::Euclidean,
            y: &y_hl,
            emat: &emat_hl,
            q: q_hl,
            nrep: nrep_hl,
            level: 0.95,
            nthreads: Some(1),
        })
        .unwrap();

        // Reference log-sum-exp evaluation without the packed fast path.
        let g = SpurGeometry::from_coords(coords_hl, n_hl, 2, Metric::Euclidean).unwrap();
        let sigdm_bm = g.sigma_lbm_dm();
        let (basis_r, _) = select_leading_basis(&sigdm_bm, &g, q_hl).unwrap();
        let mut mr = basis_r.clone();
        center_columns(&mut mr, n_hl, q_hl);
        let p_vec = prepare_y_projections(&y_hl, n_hl, 1, &basis_r, q_hl).unwrap();
        let min_distance = g.min_positive_distance();
        let compact = half_life_search_grid(min_distance);
        let ho_f: Vec<_> = half_life_covariances(&g, &mr, q_hl, &compact)
            .iter()
            .map(|omega| factor_omega(omega, q_hl).unwrap())
            .collect();
        let (half_lives, masses) = half_life_mixture(min_distance);
        let ha_f: Vec<_> = half_lives
            .iter()
            .map(|&h| {
                let c = std::f64::consts::LN_2 / h;
                factor_omega(&g.projected_omega_from_mr(c, &mr, q_hl), q_hl).unwrap()
            })
            .collect();
        let q_f = q_hl as f64;
        for (idx, ho) in ho_f.iter().enumerate() {
            let mut xc_obs = vec![0.0_f64; q_hl];
            linalg::multiply(&ho.inverse_lower, &p_vec, &mut xc_obs, q_hl, q_hl, 1, false);
            let log_ho_obs = ho.log_inverse_sqrt_det
                - 0.5 * q_f * xc_obs.iter().map(|v| v * v).sum::<f64>().ln();
            let lr_obs: f64 = ha_f
                .iter()
                .zip(&masses)
                .map(|(ha, weight)| {
                    let mut xc = vec![0.0_f64; q_hl];
                    linalg::multiply(&ha.inverse_lower, &p_vec, &mut xc, q_hl, q_hl, 1, false);
                    let lha = ha.log_inverse_sqrt_det
                        - 0.5 * q_f * xc.iter().map(|v| v * v).sum::<f64>().ln();
                    weight * (lha - log_ho_obs).exp()
                })
                .sum::<f64>();
            let mut exceed = 0usize;
            for c in 0..nrep_hl {
                let e_col: Vec<f64> = (0..q_hl).map(|r| emat_hl[r * nrep_hl + c]).collect();
                let mut x = vec![0.0_f64; q_hl];
                linalg::multiply(&ho.lower, &e_col, &mut x, q_hl, q_hl, 1, false);
                let mut z0 = vec![0.0_f64; q_hl];
                linalg::multiply(&ho.inverse_lower, &x, &mut z0, q_hl, q_hl, 1, false);
                let log_ho = ho.log_inverse_sqrt_det
                    - 0.5 * q_f * z0.iter().map(|v| v * v).sum::<f64>().ln();
                let avg: f64 = ha_f
                    .iter()
                    .zip(&masses)
                    .map(|(ha, weight)| {
                        let mut z1 = vec![0.0_f64; q_hl];
                        linalg::multiply(&ha.inverse_lower, &x, &mut z1, q_hl, q_hl, 1, false);
                        let lha = ha.log_inverse_sqrt_det
                            - 0.5 * q_f * z1.iter().map(|v| v * v).sum::<f64>().ln();
                        weight * (lha - log_ho).exp()
                    })
                    .sum::<f64>();
                if avg >= lr_obs {
                    exceed += 1;
                }
            }
            let ref_p = (exceed + 1) as f64 / (nrep_hl + 1) as f64;
            assert!(
                (out_fast.pv_mat[idx] - ref_p).abs() < 1e-12,
                "q={q_hl} idx={idx}: fast={} ref={ref_p}",
                out_fast.pv_mat[idx]
            );
        }
    }

    // 4. Tiny positive coordinate scale invariance (coords * 1e-12)
    let coords_tiny: Vec<f64> = coords_hl.iter().map(|&v| v * 1e-12).collect();
    let emat_tiny = make_deterministic_emat(6, 100);
    let hl_orig = run_spur_halflife(&SpurHalfLifeRequest {
        coords: coords_hl,
        n: n_hl,
        d: 2,
        metric: Metric::Euclidean,
        y: &y_hl,
        emat: &emat_tiny,
        q: 6,
        nrep: 100,
        level: 0.95,
        nthreads: Some(1),
    })
    .unwrap();
    let hl_tiny = run_spur_halflife(&SpurHalfLifeRequest {
        coords: &coords_tiny,
        n: n_hl,
        d: 2,
        metric: Metric::Euclidean,
        y: &y_hl,
        emat: &emat_tiny,
        q: 6,
        nrep: 100,
        level: 0.95,
        nthreads: Some(1),
    })
    .unwrap();
    assert!((hl_orig.ci_lower - hl_tiny.ci_lower).abs() < 1e-12);
    assert!((hl_orig.ci_upper - hl_tiny.ci_upper).abs() < 1e-12);
    assert!((hl_tiny.max_dist - hl_orig.max_dist * 1e-12).abs() < 1e-22);

    // A large ratio can underflow its reciprocal power despite a nonzero
    // weighted term. The fast path must mark that draw for log-space repair.
    let mut sums = [0.0];
    add_halflife_likelihood_terms(&mut sums, &[1e11], &[1.0], 1e140, |x| half_power(x, 60));
    assert!(sums[0].is_nan());
    let dominant_log = 140.0_f64 * 10.0_f64.ln() - 30.0 * (1e11_f64).ln();
    let secondary_log = 0.0_f64 - 30.0 * (1e7_f64).ln();
    let true_log_sum = log_sum_exp(&[dominant_log, secondary_log]);
    assert!(true_log_sum > -200.0 * 10.0_f64.ln());
    assert!(true_log_sum - secondary_log > 40.0);

    // Sparse and extreme cluster IDs must map cleanly without OOM or overflow.
    let y_cl = [1.0, 3.0, 10.0, 20.0];
    let tr = run_spur_transform(&SpurTransformRequest {
        kind: SpurTransformKind::Cluster,
        n: 4,
        d: 0,
        metric: Metric::Euclidean,
        coords: None,
        y: &y_cl,
        n_y: 1,
        cluster_id: Some(&[usize::MAX, 1_000_000_000, usize::MAX, 1_000_000_000]),
        nthreads: Some(1),
    })
    .unwrap();
    assert!((tr[0] - (-4.5)).abs() < 1e-12);
    assert!((tr[1] - (-8.5)).abs() < 1e-12);
    assert!((tr[2] - 4.5).abs() < 1e-12);
    assert!((tr[3] - 8.5).abs() < 1e-12);
}

#[test]
fn residual_design_accepts_a_constant_spanned_by_nearly_collinear_columns() {
    let n = 60;
    let x: Vec<_> = (1..=n)
        .flat_map(|i| {
            let scaled = 1e6 * i as f64;
            [scaled + 1.0, -scaled]
        })
        .collect();
    assert!(prepare_resid_design(&x, n, 2, 6).is_ok());
}

#[test]
fn residual_tests_require_constant_span_and_are_invariant_to_x_beta_shifts() {
    let n = 40usize;
    let q = 6usize;
    let nrep = 120usize;
    let coords: Vec<f64> = (0..n)
        .flat_map(|i| [((i * 17) % 41) as f64 / 41.0, ((i * 29) % 43) as f64 / 43.0])
        .collect();
    let x: Vec<f64> = (0..n).map(|i| 1.5 + (i as f64 * 0.23).sin()).collect();
    let x_in: Vec<f64> = x.iter().flat_map(|&xi| [xi, 1.0 - xi]).collect();
    let eps: Vec<f64> = (0..n)
        .map(|i| (i as f64 * 0.37).cos() + 0.1 * (i as f64))
        .collect();
    let shifted_y: Vec<f64> = eps.iter().zip(&x).map(|(&e, &x)| e + 1e6 * x).collect();
    let emat: Vec<f64> = (0..2 * q * nrep)
        .map(|k| (((k * 31 + 7) % 97) as f64 - 48.0) / 25.0)
        .collect();

    let req = SpurTestRequest {
        coords: &coords,
        n,
        d: 2,
        metric: Metric::Euclidean,
        y: &eps,
        n_y: 1,
        x_in: Some(&x),
        k_x: 1,
        emat: &emat,
        q,
        nrep,
        nthreads: Some(1),
    };
    for err in [
        run_spur_i0_resid(&req).unwrap_err(),
        run_spur_i1_resid(&req).unwrap_err(),
        run_spur_pipeline(&SpurPipelineRequest {
            coords: &coords,
            n,
            d: 2,
            metric: Metric::Euclidean,
            y: &eps,
            x_in: &x,
            k_x: 1,
            emat: &emat,
            q,
            nrep,
            transform_vars: &x,
            n_transform_vars: 1,
            nthreads: Some(1),
        })
        .unwrap_err(),
    ] {
        assert_eq!(err.kind(), ErrorKind::InvalidInput);
        assert_eq!(
            err.message(),
            "SPUR residual tests require a design that spans a constant"
        );
    }

    let req = SpurTestRequest {
        x_in: Some(&x_in),
        k_x: 2,
        ..req
    };
    let base_i0 = run_spur_i0_resid(&req).unwrap();
    let shift_i0 = run_spur_i0_resid(&SpurTestRequest {
        y: &shifted_y,
        ..req
    })
    .unwrap();
    assert!(
        (base_i0.lr[0] - shift_i0.lr[0]).abs() <= 1e-8 * base_i0.lr[0].abs().max(1.0),
        "base={}, shift={}",
        base_i0.lr[0],
        shift_i0.lr[0]
    );
    assert!((base_i0.pvalue[0] - shift_i0.pvalue[0]).abs() < 1e-12);

    let base_i1 = run_spur_i1_resid(&req).unwrap();
    let shift_i1 = run_spur_i1_resid(&SpurTestRequest {
        y: &shifted_y,
        ..req
    })
    .unwrap();
    assert!(
        (base_i1.lr[0] - shift_i1.lr[0]).abs() <= 1e-8 * base_i1.lr[0].abs().max(1.0),
        "base={}, shift={}",
        base_i1.lr[0],
        shift_i1.lr[0]
    );
    assert!((base_i1.pvalue[0] - shift_i1.pvalue[0]).abs() < 1e-12);

    // Response strictly inside span(X) must be rejected cleanly.
    let err = run_spur_i0_resid(&SpurTestRequest { y: &x, ..req }).unwrap_err();
    assert!(err.message().contains("zero variation"));
}

#[test]
fn half_life_search_refines_disconnected_components() {
    let grid = [0.0, 0.25, 0.5, 0.75, 1.0];
    let accept = |t: f64| t <= 0.1 || (0.4..=0.6).contains(&t) || t >= 0.9;
    let found = invert_half_life(&grid, &grid.map(accept), |points| {
        Ok(points.iter().map(|&p| accept(p)).collect())
    })
    .unwrap();
    assert_eq!(found.len(), 3);
    assert_eq!(found[0][0], 0.0);
    assert!((found[0][1] - 1.0 / 9.0).abs() < 1e-6);
    assert!((found[1][0] - 2.0 / 3.0).abs() < 1e-6);
    assert!((found[1][1] - 1.5).abs() < 1e-6);
    assert!((found[2][0] - 9.0).abs() < 2e-5);
    assert_eq!(found[2][1], f64::INFINITY);
    let (half_lives, weights) = half_life_mixture(1e-8);
    assert!((weights.iter().sum::<f64>() - 1.0).abs() < 1e-15);
    let mean = half_lives
        .iter()
        .zip(weights)
        .map(|(&h, w)| h * w)
        .sum::<f64>();
    assert!((mean - 0.5).abs() < 1e-15);
}

#[test]
fn half_life_mixture_resolves_the_close_pair_contrast() {
    let sites = [0.0, 1e-8, 0.3, 0.6, 1.0];
    let geom = SpurGeometry::from_coords(&sites, 5, 1, Metric::Euclidean).unwrap();
    // The four Helmert contrasts form an orthonormal basis of the centered space.
    let mut basis = [0.0; 5 * 4];
    for j in 0..4 {
        let span = (j + 1) as f64;
        let scale = (span * (span + 1.0)).sqrt().recip();
        for i in 0..=j {
            basis[i * 4 + j] = scale;
        }
        basis[(j + 1) * 4 + j] = -span * scale;
    }
    let observed = [2.0_f64.sqrt(), 0.0, 0.0, 0.0];
    let (half_lives, weights) = half_life_mixture(geom.min_positive_distance());
    let points: Vec<_> = half_lives.iter().map(|h| h / (1.0 + h)).collect();
    let mixture: f64 = half_life_covariances(&geom, &basis, 4, &points)
        .iter()
        .zip(weights)
        .map(|(omega, weight)| {
            let factor = factor_omega(omega, 4).unwrap();
            let form = quadratic_form_from_lower(&factor.inverse_lower, &observed, 1, 0);
            weight * factor.log_inverse_sqrt_det.exp() / form.powi(2)
        })
        .sum();
    // Independent adaptive integration of the angular Gaussian density in h.
    assert!((mixture / 3.632914971239679e-9 - 1.0).abs() < 1e-5);
}
