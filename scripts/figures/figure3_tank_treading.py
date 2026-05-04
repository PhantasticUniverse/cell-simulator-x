"""Figure 3: Tank-treading frequency vs shear rate.

Reads target/tank_treading_fischer_2007.csv. Plots Keller-Skalak 1982
analytic prediction vs shear rate, with the Fischer 2007 K(λ) ∈
[0.04, 0.15] band shaded as the published reference.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
CSV = ROOT / "target" / "tank_treading_fischer_2007.csv"
OUT = ROOT / "figures" / "figure3_tank_treading.png"


def main() -> None:
    df = pd.read_csv(CSV).sort_values("shear_rate_per_sec")

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.fill_between(
        df["shear_rate_per_sec"],
        df["fischer_low_hz"],
        df["fischer_high_hz"],
        color="#bdc3c7", alpha=0.55,
        label="Fischer 2007 band (K(λ) ∈ [0.04, 0.15])",
    )
    ax.plot(
        df["shear_rate_per_sec"], df["predicted_freq_hz"],
        color="#c0392b", marker="o", lw=2.2, ms=7,
        label="Keller-Skalak 1982 (this work)",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Shear rate γ̇ (1/s)")
    ax.set_ylabel("Tank-treading frequency f_TT (Hz)")
    ax.set_title(
        "Figure 3. Tank-treading frequency vs shear rate "
        "(Keller-Skalak 1982 / Fischer 2007)",
        fontsize=11, fontweight="bold",
    )
    ax.legend(loc="lower right", framealpha=0.9)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
