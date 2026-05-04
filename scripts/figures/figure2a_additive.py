"""Figure 2A: ATP retention by additive solution over 42 days.

Reads target/storage_additive_comparator.csv (4 additives × 43 days) and
overlays ATP_mM trajectories. Annotates day-42 retention values.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
CSV = ROOT / "target" / "storage_additive_comparator.csv"
OUT = ROOT / "figures" / "figure2a_additive.png"

# Phase 14.D.1 calibrated colors per additive (warm → cool = worst → best).
ADDITIVE_COLORS = {
    "CPD": "#c0392b",
    "SAGM": "#e67e22",
    "AS-3": "#2980b9",
    "PAGGSM": "#16a085",
}


def main() -> None:
    df = pd.read_csv(CSV)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for additive, group in df.groupby("additive"):
        color = ADDITIVE_COLORS.get(additive, "grey")
        ax.plot(
            group["day"], group["atp_mM"],
            color=color, lw=2.2, label=additive,
        )
        # Annotate day-42 retention.
        d42 = group[group["day"] >= 41.5]
        if not d42.empty:
            atp42 = d42["atp_mM"].iloc[-1]
            day_x = d42["day"].iloc[-1]
            ax.annotate(
                f"{atp42:.2f} mM",
                xy=(day_x, atp42),
                xytext=(2, 0),
                textcoords="offset points",
                fontsize=9, color=color, fontweight="bold",
                va="center",
            )

    ax.set_xlabel("Storage day")
    ax.set_ylabel("ATP (mM)")
    ax.set_title(
        "Figure 2A. Day-42 ATP retention by additive (Hess 2010 / Burger 2008)",
        fontsize=12, fontweight="bold",
    )
    ax.set_xlim(0, 45)
    ax.set_ylim(0, 2.3)
    ax.legend(title="Additive", loc="upper right", framealpha=0.9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
