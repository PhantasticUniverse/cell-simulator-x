//! Multi-cell simulation world (Phase 10.5).
//!
//! `World` owns N independent [`Cell`]s and steps them in parallel via rayon.
//! It is the new canonical multi-cell API. The legacy single-cell
//! [`crate::coupling::CoupledSolver`] continues to work unchanged for
//! existing callers; new code targeting Phase 11+ GPU compute should prefer
//! `World`.
//!
//! ## Why a separate type
//!
//! `CoupledSolver` exposes an external-state API (mesh, spectrin,
//! physics_state, metabolites are passed in by reference each step). That
//! shape is incompatible with a buffer-resident multi-cell layout, where
//! the simulator needs to own all state in contiguous arrays for kernel
//! launches.
//!
//! `World` owns its `Cell`s. Each `Cell` owns its mesh, spectrin network,
//! physics state, metabolite pool, biochemistry solver, and tension
//! history. `World::step` parallelizes per-cell work via [`rayon::prelude`]
//! — Phase 10.5 lands the API; Phase 11 will retarget the per-cell hot
//! loops onto wgpu compute shaders without changing this surface.
//!
//! ## Linear-scaling exit criterion
//!
//! Phase 10.5's exit criterion is that `World::step` scales linearly in
//! the number of cells on CPU. The cells are independent (no cell-cell
//! mechanical interactions in Phase 10.5; that is Phase 13's territory),
//! so parallelization is embarrassingly parallel. The
//! `--diagnose-multi-cell` CLI flag measures this directly.

use rayon::prelude::*;

use crate::config::Parameters;
use crate::coupling::{CoupledConfig, SpectrinModulator};

mod cell;

pub use cell::Cell;

/// Opaque handle to a cell inside a [`World`].
///
/// Indices are stable across `add_cell`. `World` does not currently support
/// removing cells; if removal is added, indices will become non-monotonic
/// and consumers should treat `CellHandle` as opaque.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CellHandle(usize);

impl CellHandle {
    /// Index of this cell within the world's storage. Useful for
    /// per-cell logging; do not rely on it for cross-world stability.
    #[inline]
    pub fn index(self) -> usize { self.0 }
}

/// Multi-cell simulation world.
pub struct World {
    cells: Vec<Cell>,
    config: CoupledConfig,
    spectrin_modulator: SpectrinModulator,
    /// Latest synchronized simulation time (the cells step in lockstep).
    pub time_sec: f64,
}

impl World {
    /// Create an empty world. Cells must be added with [`World::add_cell`].
    pub fn new(config: CoupledConfig) -> Self {
        Self {
            cells: Vec::new(),
            config,
            spectrin_modulator: SpectrinModulator::default(),
            time_sec: 0.0,
        }
    }

    /// Add a fresh cell to the world. Returns a handle for later access.
    pub fn add_cell(&mut self, params: &Parameters) -> CellHandle {
        let cell = Cell::new(params, &self.config);
        let idx = self.cells.len();
        self.cells.push(cell);
        CellHandle(idx)
    }

    /// Step every cell forward by one biochemistry timestep
    /// (== `config.physics_substeps` physics substeps each).
    ///
    /// Cells are stepped in parallel via rayon. Per-cell state is fully
    /// owned by the cell, so there is no aliasing across the parallel
    /// iteration.
    pub fn step(&mut self) {
        // Pull config + modulator out of `self` so the closure does not
        // need a mutable borrow that conflicts with `cells.par_iter_mut`.
        let config = &self.config;
        let modulator = &self.spectrin_modulator;

        self.cells.par_iter_mut().for_each(|cell| {
            cell.step(config, modulator);
        });

        if let Some(first) = self.cells.first() {
            self.time_sec = first.time_sec;
        }
    }

    /// Run for `duration_sec` seconds of simulated time.
    pub fn run(&mut self, duration_sec: f64) {
        let n_steps = (duration_sec / self.config.biochem_dt_sec).ceil() as usize;
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// Borrow a cell.
    #[inline]
    pub fn cell(&self, h: CellHandle) -> &Cell {
        &self.cells[h.0]
    }

    /// Borrow a cell mutably.
    #[inline]
    pub fn cell_mut(&mut self, h: CellHandle) -> &mut Cell {
        &mut self.cells[h.0]
    }

    /// All cells, read-only.
    #[inline]
    pub fn cells(&self) -> &[Cell] { &self.cells }

    /// All cells, mutable.
    #[inline]
    pub fn cells_mut(&mut self) -> &mut [Cell] { &mut self.cells }

    /// Number of cells in the world.
    #[inline]
    pub fn len(&self) -> usize { self.cells.len() }

    /// True iff no cells have been added.
    #[inline]
    pub fn is_empty(&self) -> bool { self.cells.is_empty() }

    /// Read-only access to the configuration.
    #[inline]
    pub fn config(&self) -> &CoupledConfig { &self.config }

    /// Read-only access to the spectrin modulator (for diagnostics).
    #[inline]
    pub fn spectrin_modulator(&self) -> &SpectrinModulator { &self.spectrin_modulator }

    /// Reset all cells. Time is reset to 0; cell positions and
    /// metabolite pools are *not* reset (use re-create for that).
    pub fn reset(&mut self) {
        for cell in &mut self.cells {
            cell.reset();
        }
        self.time_sec = 0.0;
    }

    /// Apply an override tension to every cell (debugging / testing).
    pub fn set_tension_override(&mut self, tension_pN_per_nm: f64) {
        for cell in &mut self.cells {
            cell.set_tension_override(tension_pN_per_nm);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::{GeometryParameters, Parameters};

    fn small_params() -> Parameters {
        Parameters {
            geometry: GeometryParameters {
                cell_radius_um: 3.91,
                fung_tong_c0_um: 0.81,
                fung_tong_c2_um: 7.83,
                fung_tong_c4_um: -4.39,
                mesh_resolution: 8,
                spectrin_target_count: 50,
            },
            ..Parameters::default()
        }
    }

    fn small_config() -> CoupledConfig {
        CoupledConfig {
            physics_substeps: 5,
            ..CoupledConfig::default()
        }
    }

    #[test]
    fn world_starts_empty() {
        let world = World::new(small_config());
        assert!(world.is_empty());
        assert_eq!(world.len(), 0);
    }

    #[test]
    fn add_cell_returns_handle() {
        let mut world = World::new(small_config());
        let params = small_params();
        let h = world.add_cell(&params);
        assert_eq!(h.index(), 0);
        assert_eq!(world.len(), 1);
    }

    #[test]
    fn world_steps_n_equals_one() {
        let mut world = World::new(small_config());
        let params = small_params();
        let h = world.add_cell(&params);
        world.step();
        assert!(world.cell(h).time_sec > 0.0);
        assert!(world.time_sec > 0.0);
    }

    #[test]
    fn world_steps_n_equals_ten() {
        let mut world = World::new(small_config());
        let params = small_params();
        for _ in 0..10 {
            world.add_cell(&params);
        }
        world.step();
        assert_eq!(world.len(), 10);
        for cell in world.cells() {
            assert!(cell.time_sec > 0.0);
        }
    }

    #[test]
    fn world_handles_independent() {
        // Each cell maintains its own state; perturbing one must not
        // affect another.
        let mut world = World::new(small_config());
        let params = small_params();
        let h0 = world.add_cell(&params);
        let h1 = world.add_cell(&params);

        // Apply a tension override only to cell 0.
        let atp_idx = world.cell(h0).biochemistry.indices.glycolysis.atp;
        world.cell_mut(h0).metabolites.set(atp_idx, 0.5);

        let atp0 = world.cell(h0).metabolites.get(atp_idx);
        let atp1 = world.cell(h1).metabolites.get(atp_idx);
        assert!((atp0 - 0.5).abs() < 1e-9);
        assert!((atp1 - 2.0).abs() < 1e-9);
    }
}
