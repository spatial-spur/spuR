//! Controlled Rayon thread-pool dispatch for host language adapters.

use crate::CoreError;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::cell::Cell;
use std::sync::{Arc, Mutex};

static CACHED_POOL: Mutex<Option<(usize, Arc<ThreadPool>)>> = Mutex::new(None);

thread_local! {
    static SERIAL_MODE: Cell<bool> = const { Cell::new(false) };
}

/// Returns `true` while executing inside a single-thread `nthreads = Some(1)` scope.
///
/// Parallel sections check this flag and run sequentially on the calling thread,
/// avoiding Rayon pool locks that would deadlock after `fork()`.
pub fn is_serial() -> bool {
    SERIAL_MODE.with(Cell::get)
}

struct SerialGuard(bool);

impl SerialGuard {
    fn enter(serial: bool) -> Self {
        let prev = SERIAL_MODE.with(|flag| {
            let prev = flag.get();
            flag.set(serial);
            prev
        });
        Self(prev)
    }
}

impl Drop for SerialGuard {
    fn drop(&mut self) {
        SERIAL_MODE.with(|flag| flag.set(self.0));
    }
}

/// Execute `f` using the requested thread budget.
///
/// Setting `nthreads = Some(1)` executes the compiled Rust routine sequentially
/// directly on the caller thread without locking or spawning Rayon pools, making
/// single-threaded calls fork-safe even if multi-threaded calls ran previously.
pub fn run_with_threads<R, F>(nthreads: Option<usize>, f: F) -> Result<R, CoreError>
where
    R: Send,
    F: FnOnce() -> R + Send,
{
    let Some(count) = nthreads else {
        let _guard = SerialGuard::enter(false);
        return Ok(f());
    };
    if count == 0 {
        return Err(CoreError::invalid_input("thread count must be at least 1"));
    }
    if count == 1 {
        let _guard = SerialGuard::enter(true);
        return Ok(f());
    }
    let pool = {
        let mut guard = CACHED_POOL
            .lock()
            .map_err(|_| CoreError::internal("thread pool mutex was poisoned"))?;
        if let Some((cached_count, ref pool)) = *guard {
            if cached_count == count {
                Arc::clone(pool)
            } else {
                let built = Arc::new(
                    ThreadPoolBuilder::new()
                        .num_threads(count)
                        .build()
                        .map_err(|_| CoreError::internal("failed to initialize Rayon thread pool"))?,
                );
                *guard = Some((count, Arc::clone(&built)));
                built
            }
        } else {
            let built = Arc::new(
                ThreadPoolBuilder::new()
                    .num_threads(count)
                    .build()
                    .map_err(|_| CoreError::internal("failed to initialize Rayon thread pool"))?,
            );
            *guard = Some((count, Arc::clone(&built)));
            built
        }
    };
    Ok(pool.install(|| {
        let _guard = SerialGuard::enter(false);
        f()
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_thread_budget_runs_serially_on_caller_thread() {
        let caller = std::thread::current().id();
        let (inside_serial, inside_tid) =
            run_with_threads(Some(1), || (is_serial(), std::thread::current().id())).unwrap();
        assert!(inside_serial);
        assert_eq!(inside_tid, caller);
        assert!(!is_serial());

        let (multi_serial, seen_two) =
            run_with_threads(Some(2), || (is_serial(), rayon::current_num_threads())).unwrap();
        assert!(!multi_serial);
        assert_eq!(seen_two, 2);

        // Even after warming a multi-thread pool, Some(1) executes strictly on caller thread.
        let (after_serial, after_tid) =
            run_with_threads(Some(1), || (is_serial(), std::thread::current().id())).unwrap();
        assert!(after_serial);
        assert_eq!(after_tid, caller);
    }
}
