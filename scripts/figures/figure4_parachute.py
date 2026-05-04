"""Figure 4: Parachute aspect ratio vs shear rate / capillary number.

Reads target/parachute_skalak_1973.csv. Plots predicted aspect ratio vs
shear rate (with capillary number annotation) and shades the Skalak 1973
reported band [1.5, 2.0] as the published reference.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
CSV = ROOT / "target" / "parachute_skalak_1973.csv"
OUT = ROOT / "figures" / "figure4_parachute.png"


def main() -> None:
    df = pd.read_csv(CSV).sort_values("shear_rate_per_sec")

    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.fill_between(
        df["shear_rate_per_sec"],
        df["skalak_low"],
        df["skalak_high"],
        color="#bdc3c7", alpha=0.55,
        label="Skalak 1973 band ([1.5, 2.0])",
    )
    ax.plot(
        df["shear_rate_per_sec"], df["predicted_aspect_ratio"],
        color="#2980b9", marker="o", lw=2.2, ms=7,
        label="Empirical fit AR = 1 + 1.5·tanh(3·Ca) (this work)",
    )
    # Annotate Ca on each marker.
    for _, row in df.iterrows():
        ax.annotate(
            f"Ca={row['capillary_number']:.2f}",
            xy=(row["shear_rate_per_sec"], row["predicted_aspect_ratio"]),
            xytext=(6, -10), textcoords="offset points",
            fontsize=8, color="#34495e",
        )

    ax.set_xlabel("Wall shear rate γ̇_wall (1/s)")
    ax.set_ylabel("Parachute aspect ratio (axial / transverse)")
    ax.set_title(
        "Figure 4. Parachute aspect ratio under capillary Poiseuille flow "
        "(Skalak 1973)",
        fontsize=11, fontweight="bold",
    )
    ax.set_ylim(1.0, 2.6)
    ax.legend(loc="lower right", framealpha=0.9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
