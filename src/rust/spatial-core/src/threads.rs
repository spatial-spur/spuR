//! Controlled Rayon thread-pool dispatch for host language adapters.

use crate::CoreError;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::cell::Cell;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{Arc, Mutex};

static CACHED_POOL: Mutex<Option<(u32, usize, Arc<ThreadPool>)>> = Mutex::new(None);
static POOL_OWNER_PID: AtomicU32 = AtomicU32::new(0);
#[cfg(all(target_vendor = "apple", not(feature = "portable")))]
static ACCELERATE_OWNER_PID: AtomicU32 = AtomicU32::new(0);

pub const MAX_THREADS: usize = 1024;

pub fn validate_nthreads(nthreads: Option<usize>) -> Result<(), CoreError> {
    if nthreads.is_some_and(|n| !(1..=MAX_THREADS).contains(&n)) {
        return Err(CoreError::invalid_input(
            "thread count must be between 1 and 1024",
        ));
    }
    Ok(())
}

pub fn inherited_accelerate_state() -> bool {
    #[cfg(all(target_vendor = "apple", not(feature = "portable")))]
    {
        let owner = ACCELERATE_OWNER_PID.load(Ordering::Acquire);
        owner != 0 && owner != std::process::id()
    }
    #[cfg(any(not(target_vendor = "apple"), feature = "portable"))]
    false
}

pub(crate) fn check_native_process() -> Result<(), CoreError> {
    #[cfg(all(target_vendor = "apple", not(feature = "portable")))]
    {
        let _ = ACCELERATE_OWNER_PID.compare_exchange(
            0,
            std::process::id(),
            Ordering::AcqRel,
            Ordering::Acquire,
        );
        if inherited_accelerate_state() {
            return Err(CoreError::unsupported(
                "Accelerate was used before this process forked. Use a fresh process or the host fallback.",
            ));
        }
    }
    Ok(())
}

thread_local! {
    static SERIAL_MODE: Cell<bool> = const { Cell::new(false) };
}

/// Returns `true` while executing inside a single-thread `nthreads = Some(1)` scope
/// or inside a forked child process after a Rayon pool was initialized in the parent.
///
/// Parallel sections check this flag and run sequentially on the calling thread,
/// avoiding Rayon pool locks and dead worker threads after `fork()`.
pub fn is_serial() -> bool {
    SERIAL_MODE.with(Cell::get) || is_forked_after_pool_init()
}

/// Returns `true` when this process is a forked child of a parent process that
/// already initialized a multi-thread Rayon pool.
pub fn is_forked_after_pool_init() -> bool {
    let owner = POOL_OWNER_PID.load(Ordering::Acquire);
    owner != 0 && owner != std::process::id()
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
/// directly on the caller thread without locking or spawning Rayon pools.
/// In a forked child process whose parent already initialized a multi-thread Rayon
/// pool, calls automatically degrade to serial execution on the calling thread so
/// inherited dead Rayon worker threads are never waited on or dropped.
pub fn run_with_threads<R, F>(nthreads: Option<usize>, f: F) -> Result<R, CoreError>
where
    R: Send,
    F: FnOnce() -> R + Send,
{
    validate_nthreads(nthreads)?;
    check_native_process()?;
    if matches!(nthreads, Some(1)) || is_forked_after_pool_init() {
        let _guard = SerialGuard::enter(true);
        return Ok(f());
    }

    let current_pid = std::process::id();
    if let Err(owner) =
        POOL_OWNER_PID.compare_exchange(0, current_pid, Ordering::AcqRel, Ordering::Acquire)
    {
        if owner != current_pid {
            let _guard = SerialGuard::enter(true);
            return Ok(f());
        }
    }

    let Some(count) = nthreads else {
        let _guard = SerialGuard::enter(false);
        return Ok(f());
    };
    let pool = {
        let mut guard = CACHED_POOL
            .lock()
            .map_err(|_| CoreError::internal("thread pool mutex was poisoned"))?;
        if let Some((owner_pid, cached_count, ref pool)) = *guard {
            if owner_pid != current_pid {
                let _guard = SerialGuard::enter(true);
                return Ok(f());
            }
            if cached_count == count {
                Arc::clone(pool)
            } else {
                let built = Arc::new(
                    ThreadPoolBuilder::new()
                        .num_threads(count)
                        .build()
                        .map_err(|_| {
                            CoreError::internal("failed to initialize Rayon thread pool")
                        })?,
                );
                *guard = Some((current_pid, count, Arc::clone(&built)));
                built
            }
        } else {
            let built = Arc::new(
                ThreadPoolBuilder::new()
                    .num_threads(count)
                    .build()
                    .map_err(|_| CoreError::internal("failed to initialize Rayon thread pool"))?,
            );
            *guard = Some((current_pid, count, Arc::clone(&built)));
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

    #[cfg(unix)]
    #[test]
    fn child_of_a_warmed_process_obeys_the_provider_fork_policy() {
        extern "C" {
            fn fork() -> i32;
            fn _exit(status: i32) -> !;
            fn waitpid(pid: i32, status: *mut i32, options: i32) -> i32;
        }

        // Warm both the cached multi-thread pool and verify owner PID tracking in parent.
        let parent_threads = run_with_threads(Some(2), || rayon::current_num_threads()).unwrap();
        assert_eq!(parent_threads, 2);
        assert!(!is_forked_after_pool_init());

        let pid = unsafe { fork() };
        assert!(pid >= 0, "fork() failed");
        if pid == 0 {
            let ok = is_forked_after_pool_init()
                && [None, Some(1), Some(2)].into_iter().all(|nthreads| {
                    let out = run_with_threads(nthreads, is_serial);
                    if inherited_accelerate_state() {
                        matches!(out, Err(err) if err.kind() == crate::ErrorKind::UnsupportedCase)
                    } else {
                        matches!(out, Ok(true))
                    }
                });
            unsafe { _exit(if ok { 0 } else { 1 }) };
        }
        let mut status: i32 = -1;
        let waited = unsafe { waitpid(pid, &mut status, 0) };
        assert_eq!(waited, pid);
        assert_eq!(
            status, 0,
            "forked child exited with non-zero status: {status}"
        );
    }
}
