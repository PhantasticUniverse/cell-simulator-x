# Cell Simulator X — figure-reproduction Makefile (Phase 17.2)
#
# Each `figure-*` target runs the test that emits the corresponding CSV
# under `target/`. The `figures` aggregate target runs every figure.
#
# Plotting is left to the consumer; the CSVs are the canonical
# scientific artifact.

CARGO ?= cargo
RELEASE_FLAGS ?= --release
TARGET_DIR ?= target

.PHONY: build test figures help \
	figure-storage figure-additives figure-sensitivity \
	figure-tank-treading figure-parachute \
	figure-splenic-transit \
	figure-gpu-parity

help:
	@echo "Cell Simulator X figure-reproduction Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  figures               Reproduce all headline figure CSVs."
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

build:
	$(CARGO) build $(RELEASE_FLAGS)

test:
	$(CARGO) test $(RELEASE_FLAGS)

figures: figure-storage figure-additives figure-sensitivity \
         figure-tank-treading figure-parachute figure-splenic-transit
	@echo ""
	@echo "All figure CSVs reproduced under $(TARGET_DIR)/"

figure-storage:
	@echo "=== Figure 1: 42-day storage lesion ==="
	$(CARGO) test $(RELEASE_FLAGS) --test storage_curve \
		ion_gradients_match_hess_2010_quantitatively -- --nocapture
	$(CARGO) test $(RELEASE_FLAGS) --test storage_curve \
		write_csv_succeeds -- --nocapture

figure-additives:
	@echo "=== Figure 2A: Additive comparator (AS-3 / SAGM / PAGGSM) ==="
	$(CARGO) test $(RELEASE_FLAGS) --test storage_additive_comparator -- --nocapture

figure-sensitivity:
	@echo "=== Figure 2B: ±20% sensitivity sweep ==="
	$(CARGO) test $(RELEASE_FLAGS) --test storage_sensitivity \
		baseline_oat_sweep_emits_csv -- --nocapture
	@echo "  → $(TARGET_DIR)/storage_sensitivity.csv"

figure-tank-treading:
	@echo "=== Figure 3: Fischer 2007 tank-treading ==="
	$(CARGO) test $(RELEASE_FLAGS) --features validation \
		--lib validation::experiments::fischer_2007 -- --nocapture

figure-parachute:
	@echo "=== Figure 4: Skalak 1973 parachute ==="
	$(CARGO) test $(RELEASE_FLAGS) --features validation \
		--lib validation::experiments::skalak_1973 -- --nocapture

figure-splenic-transit:
	@echo "=== Figure 5: Splenic transit at storage day N ==="
	$(CARGO) test $(RELEASE_FLAGS) --test splenic_transit_sweep \
		transit_sweep_emits_csv -- --nocapture
	@echo "  → $(TARGET_DIR)/splenic_transit_storage_curve.csv"

figure-gpu-parity:
	@echo "=== Technical: CPU/GPU parity validation ==="
	$(CARGO) test $(RELEASE_FLAGS) --test biochem_gpu_full_solver_parity -- --nocapture
	$(CARGO) test $(RELEASE_FLAGS) --test physics_backend_parity -- --nocapture
