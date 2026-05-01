//! Per-experiment validation routines.
//!
//! Each submodule reproduces a published reference dataset and compares
//! the simulator's prediction against it. The functions here are the
//! entry points referenced by [`crate::validation::run_full_suite`].

pub mod dao_2003;
pub mod imai_1981;
pub mod mulquiney_1999;
pub mod rief_1999;
pub mod waugh_evans_1979;
