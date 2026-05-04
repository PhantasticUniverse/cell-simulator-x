# Cell Simulator X — figure-reproduction Makefile (Phase 17.2 + 17.4)
#
# Each `figure-*` target runs the test that emits the corresponding CSV
# under `target/`, then renders the matching PNG into `figures/` via the
# Python plotting script in `scripts/figures/`.
#
# Plotting requires `pandas>=2.0` and `matplotlib>=3.7` — install via
# `pip install -r scripts/figures/requirements.txt`.

CARGO ?= cargo
RELEASE_FLAGS ?= --release
TARGET_DIR ?= target
FIGURES_DIR ?= figures
PYTHON ?= python3
SCRIPTS_DIR ?= scripts/figures

.PHONY: build test figures figures-png help \
	figure-storage figure-additives figure-sensitivity \
	figure-tank-treading figure-parachute \
	figure-splenic-transit \
	figure-gpu-parity \
	figures-dir

help:
	@echo "Cell Simulator X figure-reproduction Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  figures               Reproduce all headline figure CSVs + PNGs."
	@echo "  figures-png           Re-render PNGs from existing CSVs (no cargo)."
	@echo ""
	@echo "  figure-storage        Figure 1: 42-day storage lesion (Hess 2010 ions/ATP/DPG)."
	@echo "  figure-additives      Figure 2A: AS-3 / SAGM / PAGGSM additive comparator."
	@echo "  figure-sensitivity    Figure 2B: ±20% storage envelope sensitivity sweep."
	@echo "  figure-tank-treading  Figure 3: Fischer 2007 tank-treading."
	@echo "  figure-parachute      Figure 4: Skalak 1973 parachute aspect ratio."
	@echo "  figure-splenic-transit Figure 5: Splenic transit at storage day N (the demo)."
	@echo "  figure-gpu-parity     Technical fig: CPU/GPU parity tolerances."
	@echo ""
	@echo "  build                 cargo build --release."
	@echo "  test                  cargo test --release (full regression)."
	@echo ""
	@echo "Plotting prereq: pip install -r $(SCRIPTS_DIR)/requirements.txt"

build:
	$(CARGO) build $(RELEASE_FLAGS)

test:
	$(CARGO) test $(RELEASE_FLAGS)

figures-dir:
	@mkdir -p $(FIGURES_DIR)

figures: figure-storage figure-additives figure-sensitivity \
         figure-tank-treading figure-parachute figure-splenic-transit
	@echo ""
	@echo "All figure CSVs reproduced under $(TARGET_DIR)/ and PNGs under $(FIGURES_DIR)/"

# Aggregate target: re-render every PNG from existing CSVs without
# re-running cargo. Useful when only the plotting code changed.
figures-png: figures-dir
	$(PYTHON) $(SCRIPTS_DIR)/figure1_storage.py
	$(PYTHON) $(SCRIPTS_DIR)/figure2a_additive.py
	$(PYTHON) $(SCRIPTS_DIR)/figure2b_sensitivity.py
	$(PYTHON) $(SCRIPTS_DIR)/figure3_tank_treading.py
	$(PYTHON) $(SCRIPTS_DIR)/figure4_parachute.py
	$(PYTHON) $(SCRIPTS_DIR)/figure5_splenic_transit.py

figure-storage: figures-dir
	@echo "=== Figure 1: 42-day storage lesion ==="
	$(CARGO) test $(RELEASE_FLAGS) --test storage_curve \
		ion_gradients_match_hess_2010_quantitatively -- --nocapture
	$(CARGO) test $(RELEASE_FLAGS) --test storage_curve \
		write_csv_succeeds -- --nocapture
	$(PYTHON) $(SCRIPTS_DIR)/figure1_storage.py

figure-additives: figures-dir
	@echo "=== Figure 2A: Additive comparator (AS-3 / SAGM / PAGGSM) ==="
	$(CARGO) test $(RELEASE_FLAGS) --test storage_additive_comparator \
		additive_comparator_emits_csv -- --nocapture
	$(PYTHON) $(SCRIPTS_DIR)/figure2a_additive.py

figure-sensitivity: figures-dir
	@echo "=== Figure 2B: ±20% sensitivity sweep ==="
	$(CARGO) test $(RELEASE_FLAGS) --test storage_sensitivity \
		baseline_oat_sweep_emits_csv -- --nocapture
	@echo "  → $(TARGET_DIR)/storage_sensitivity.csv"
	$(PYTHON) $(SCRIPTS_DIR)/figure2b_sensitivity.py

figure-tank-treading: figures-dir
	@echo "=== Figure 3: Fischer 2007 tank-treading ==="
	$(CARGO) test $(RELEASE_FLAGS) --test headline_validation_csvs \
		tank_treading_emits_csv -- --nocapture
	$(PYTHON) $(SCRIPTS_DIR)/figure3_tank_treading.py

figure-parachute: figures-dir
	@echo "=== Figure 4: Skalak 1973 parachute ==="
	$(CARGO) test $(RELEASE_FLAGS) --test headline_validation_csvs \
		parachute_emits_csv -- --nocapture
	$(PYTHON) $(SCRIPTS_DIR)/figure4_parachute.py

figure-splenic-transit: figures-dir
	@echo "=== Figure 5: Splenic transit at storage day N ==="
	$(CARGO) test $(RELEASE_FLAGS) --test splenic_transit_sweep \
		transit_sweep_emits_csv -- --nocapture
	@echo "  → $(TARGET_DIR)/splenic_transit_storage_curve.csv"
	$(PYTHON) $(SCRIPTS_DIR)/figure5_splenic_transit.py

figure-gpu-parity:
	@echo "=== Technical: CPU/GPU parity validation ==="
	$(CARGO) test $(RELEASE_FLAGS) --test biochem_gpu_full_solver_parity -- --nocapture
	$(CARGO) test $(RELEASE_FLAGS) --test physics_backend_parity -- --nocapture
