//! Phase 11.4 test: `ComputeContext::from_shared` produces a context that
//! can dispatch a kernel using a borrowed `Arc<Device>`/`Arc<Queue>`.
//!
//! Skipped if no GPU adapter is available.

use cell_simulator_x::compute::{vec_add_gpu, ComputeContext};

#[test]
fn from_shared_dispatches_correctly() {
    let owner = match ComputeContext::new_headless_blocking() {
        Ok(c) => c,
        Err(e) => {
            eprintln!("skipping GPU test: no adapter: {e}");
            return;
        }
    };

    // Borrow the device + queue and rebuild a context.
    let shared = ComputeContext::from_shared(owner.device.clone(), owner.queue.clone());
    assert!(shared.adapter.is_none(), "shared context must not own an adapter");

    // Dispatch a vec_add through the shared context.
    let n = 1024;
    let a: Vec<f32> = (0..n).map(|i| i as f32).collect();
    let b: Vec<f32> = (0..n).map(|i| (i as f32) * 2.0).collect();
    let out = vec_add_gpu(&shared, &a, &b).expect("vec_add_gpu");

    assert_eq!(out.len(), n);
    for i in 0..n {
        assert!((out[i] - (a[i] + b[i])).abs() < 1e-6);
    }

    // Adapter summary should be the borrowed-context message.
    let summary = shared.adapter_summary();
    assert!(
        summary.contains("borrowed"),
        "expected 'borrowed' in adapter_summary, got: {summary}"
    );
}
