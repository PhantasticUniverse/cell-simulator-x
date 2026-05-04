"""Figure 5: Splenic transit displacement at storage day × slit width.

Reads target/splenic_transit_storage_curve.csv (21 rows: 7 days × 3
widths) and plots a grouped bar chart of centroid displacement. The
finding from Stream A (drag-sensitivity sweep): the system is
stiffness-saturated at the current SpectrinModulator
(max_stiffening_factor=0.5); slit-width effect dominates.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
CSV = ROOT / "target" / "splenic_transit_storage_curve.csv"
OUT = ROOT / "figures" / "figure5_splenic_transit.png"

WIDTH_COLORS = {
    0.5: "#c0392b",
    0.7: "#e67e22",
    1.0: "#2980b9",
}


def main() -> None:
    df = pd.read_csv(CSV)

    days = sorted(df["storage_day"].unique())
    widths = sorted(df["slit_width_um"].unique())
    n_days = len(days)
    n_widths = len(widths)
    bar_w = 0.8 / n_widths

    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(n_days)
    for i, w in enumerate(widths):
        sub = df[df["slit_width_um"] == w].sort_values("storage_day")
        ax.bar(
            x + (i - (n_widths - 1) / 2) * bar_w,
            sub["centroid_displacement_um"],
            width=bar_w,
            color=WIDTH_COLORS.get(w, "grey"),
            label=f"slit = {w:.1f} μm",
        )

    ax.set_xticks(list(x))
    ax.set_xticklabels([f"{int(d)}" for d in days])
    ax.set_xlabel("Storage day")
    ax.set_ylabel("Centroid displacement (μm)")
    ax.set_title(
        "Figure 5. Splenic-slit transit is stiffness-saturated at current\n"
        "SpectrinModulator (max_stiffening_factor=0.5); slit-width effect dominates",
        fontsize=11, fontweight="bold",
    )
    ax.legend(title="Slit width", loc="upper right", framealpha=0.9)
    ax.grid(alpha=0.3, axis="y")

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
