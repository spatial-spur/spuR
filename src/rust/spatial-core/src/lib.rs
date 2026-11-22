//! Shared numerical core for spatial SCPC and SPUR.
//!
//! Matrices at this boundary are contiguous and row-major.

pub mod basis;
pub mod calibration;
mod covariance;
pub mod distance;
pub mod linalg;
pub mod nystrom;
pub mod preparation;
mod quadrature;
mod scpc_op;
pub mod spur_op;
mod student_t;
pub mod threads;

pub use distance::Metric;
pub use preparation::{run_model_scpc, ScpcPreparation};
pub use scpc_op::ScpcOutput;
pub use spur_op::{
    run_spur_halflife, run_spur_i0, run_spur_i0_resid, run_spur_i1, run_spur_i1_resid,
    run_spur_pipeline, run_spur_transform, SpurHalfLifeOutput, SpurHalfLifeRequest, SpurI0Output,
    SpurI1Output, SpurPipelineOutput, SpurPipelineRequest, SpurTestRequest, SpurTransformKind,
    SpurTransformRequest, HALF_LIFE_MAX_COMPONENTS,
};

pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const SPEC_VERSION: &str = "0.3.7";
pub const SCPC_API_VERSION: usize = 7;
pub const SPUR_API_VERSION: usize = 5;

/// Returns the active linear algebra provider compiled into `spatial-core`.
pub const fn linalg_provider() -> &'static str {
    #[cfg(all(target_vendor = "apple", not(feature = "portable")))]
    {
        "accelerate"
    }
    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    {
        "faer"
    }
}

use std::fmt;

/// Categorized error kind so host adapters can distinguish recoverable states
/// (e.g. tied eigenspaces, degenerate variance, calibration failure) from bugs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrorKind {
    InvalidInput,
    UnsupportedCase,
    TiedEigenspace,
    InvalidSampledCovariance,
    CalibrationFailed,
    DegenerateVariance,
    InternalFailure,
}

/// A numerical input or calculation error, translated by each host adapter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoreError {
    pub kind: ErrorKind,
    pub message: &'static str,
}

impl CoreError {
    pub const fn new(kind: ErrorKind, message: &'static str) -> Self {
        Self { kind, message }
    }

    pub const fn invalid_input(message: &'static str) -> Self {
        Self::new(ErrorKind::InvalidInput, message)
    }

    pub const fn unsupported(message: &'static str) -> Self {
        Self::new(ErrorKind::UnsupportedCase, message)
    }

    pub const fn tied(message: &'static str) -> Self {
        Self::new(ErrorKind::TiedEigenspace, message)
    }

    pub const fn invalid_sampled_covariance(message: &'static str) -> Self {
        Self::new(ErrorKind::InvalidSampledCovariance, message)
    }

    pub const fn calibration(message: &'static str) -> Self {
        Self::new(ErrorKind::CalibrationFailed, message)
    }

    pub const fn degenerate(message: &'static str) -> Self {
        Self::new(ErrorKind::DegenerateVariance, message)
    }

    pub const fn internal(message: &'static str) -> Self {
        Self::new(ErrorKind::InternalFailure, message)
    }

    pub const fn kind(&self) -> ErrorKind {
        self.kind
    }

    pub const fn message(&self) -> &'static str {
        self.message
    }
}

impl fmt::Display for CoreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.message.fmt(f)
    }
}

impl std::error::Error for CoreError {}

pub(crate) fn checked_len(dimensions: &[usize]) -> Result<usize, CoreError> {
    dimensions
        .iter()
        .try_fold(1usize, |size, dimension| size.checked_mul(*dimension))
        .filter(|&size| size <= isize::MAX as usize / size_of::<f64>())
        .ok_or(CoreError::invalid_input("matrix dimensions are too large"))
}
