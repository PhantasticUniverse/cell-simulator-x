//! Phase 11.5 — runtime backend selection via the `WORLD_BACKEND`
//! environment variable.
//!
//! The simulator can run physics on either the CPU (`PhysicsSolver`) or
//! the GPU (`compute::PhysicsBackend`). For the moment, multi-cell `World`
//! still steps on CPU; the env-var infrastructure here lets future
//! callers switch backends at runtime without recompiling.
//!
//! ## Usage
//!
//! ```rust,ignore
//! use cell_simulator_x::world::backend::{Backend, env_default};
//!
//! match env_default() {
//!     Backend::Cpu => { /* CPU path */ }
//!     Backend::Gpu => { /* GPU path */ }
//! }
//! ```

/// Backend choice for physics + biochemistry stepping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Backend {
    /// CPU-side `PhysicsSolver` + `FullyIntegratedSolver` per cell.
    Cpu,
    /// GPU-side `PhysicsBackend` (Phase 11.3.E) + `compute::run_full_biochem_batch_with_hb`
    /// (Phase 11.2.E) per dispatch.
    Gpu,
}

impl Default for Backend {
    fn default() -> Self {
        Backend::Cpu
    }
}

impl Backend {
    /// Parse a backend name. Accepts "cpu" or "gpu" case-insensitively.
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "cpu" => Some(Backend::Cpu),
            "gpu" => Some(Backend::Gpu),
            _ => None,
        }
    }
}

/// The environment variable name. Constants like this are exported so
/// tests and tools can document/validate the contract.
pub const ENV_VAR: &str = "WORLD_BACKEND";

/// Read `WORLD_BACKEND` from the environment, returning `Backend::Cpu`
/// when unset or unparseable.
pub fn env_default() -> Backend {
    match std::env::var(ENV_VAR) {
        Ok(s) => Backend::parse(&s).unwrap_or_default(),
        Err(_) => Backend::default(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_known_values() {
        assert_eq!(Backend::parse("cpu"), Some(Backend::Cpu));
        assert_eq!(Backend::parse("CPU"), Some(Backend::Cpu));
        assert_eq!(Backend::parse(" gpu "), Some(Backend::Gpu));
        assert_eq!(Backend::parse("Gpu"), Some(Backend::Gpu));
    }

    #[test]
    fn parse_unknown_returns_none() {
        assert_eq!(Backend::parse(""), None);
        assert_eq!(Backend::parse("whatever"), None);
    }

    #[test]
    fn default_is_cpu() {
        assert_eq!(Backend::default(), Backend::Cpu);
    }
}
