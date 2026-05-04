"""Figure 2B: ±20% storage envelope sensitivity tornado plot.

Reads target/storage_sensitivity.csv (10 rows: 5 parameters × ±20%) and
shows the relative change in day-42 ATP from baseline as horizontal bars,
sorted by absolute magnitude. The dominant parameter (ATP half-life)
sits at the top.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
CSV = ROOT / "target" / "storage_sensitivity.csv"
OUT = ROOT / "figures" / "figure2b_sensitivity.png"


def main() -> None:
    df = pd.read_csv(CSV)

    # Pivot so each parameter has a -20% and +20% bar.
    pivot = df.pivot_table(
        index="parameter",
        columns="perturbation_pct",
        values="day42_atp_relative_change",
    )
    # Sort by max-abs swing across the two perturbations.
    pivot["abs_max"] = pivot.abs().max(axis=1)
    pivot = pivot.sort_values("abs_max", ascending=True)
    pivot = pivot.drop(columns=["abs_max"])

    params = list(pivot.index)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    y = range(len(params))
    minus = pivot[-20.0].values * 100.0
    plus = pivot[20.0].values * 100.0

    ax.barh(y, minus, color="#c0392b", alpha=0.85, label="-20%")
    ax.barh(y, plus, color="#2980b9", alpha=0.85, label="+20%")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(params, fontsize=10)
    ax.set_xlabel("Day-42 ATP relative change (%)")
    ax.set_title(
        "Figure 2B. ±20% storage envelope sensitivity — tornado plot",
        fontsize=12, fontweight="bold",
    )
    ax.legend(loc="lower right", framealpha=0.9)
    ax.grid(alpha=0.3, axis="x")

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
