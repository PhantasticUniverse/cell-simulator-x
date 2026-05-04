"""Figure 1: 42-day storage lesion in CPD baseline.

Reads target/storage_curve.csv and produces a four-panel time-series of
ATP, 2,3-DPG, ion gradients (Na+, K+), and relative deformability over
the 42-day storage window. Reproduces Hess 2010 / Pivkin 2011 envelopes.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
CSV = ROOT / "target" / "storage_curve.csv"
OUT = ROOT / "figures" / "figure1_storage.png"


def main() -> None:
    df = pd.read_csv(CSV)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    fig.suptitle(
        "42-day storage lesion in CPD baseline (Hess 2010 / Pivkin 2011)",
        fontsize=13, fontweight="bold",
    )

    # Panel A: ATP
    ax = axes[0, 0]
    ax.plot(df["day"], df["atp_mM"], color="#c0392b", lw=2.0, label="ATP")
    ax.axhline(2.0, color="grey", ls=":", lw=0.8, alpha=0.6)
    ax.set_ylabel("ATP (mM)")
    ax.set_title("A. ATP — Hess 2010 21-day half-life")
    ax.set_ylim(0, 2.3)
    ax.grid(alpha=0.3)

    # Panel B: 2,3-DPG
    ax = axes[0, 1]
    ax.plot(df["day"], df["dpg23_mM"], color="#2980b9", lw=2.0, label="2,3-DPG")
    ax.set_ylabel("2,3-DPG (mM)")
    ax.set_title("B. 2,3-DPG — Zimrin 2009")
    ax.set_ylim(0, 5.5)
    ax.grid(alpha=0.3)

    # Panel C: Ion gradients
    ax = axes[1, 0]
    ax.plot(df["day"], df["na_cyt_mM"], color="#e67e22", lw=2.0, label="Na+")
    ax.plot(df["day"], df["k_cyt_mM"], color="#16a085", lw=2.0, label="K+")
    ax.set_xlabel("Storage day")
    ax.set_ylabel("Cytosolic concentration (mM)")
    ax.set_title("C. Ion gradients — Hess 2010 / Luten 2008")
    ax.legend(loc="center right", framealpha=0.9)
    ax.grid(alpha=0.3)

    # Panel D: Deformability
    ax = axes[1, 1]
    ax.plot(df["day"], df["deformability_relative"], color="#8e44ad", lw=2.0)
    ax.set_xlabel("Storage day")
    ax.set_ylabel("Deformability (relative to day 0)")
    ax.set_title("D. Deformability — Pivkin 2011 ~30% decline")
    ax.set_ylim(0.65, 1.05)
    ax.grid(alpha=0.3)

    for ax in axes.flat:
        ax.set_xlim(0, 42)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
