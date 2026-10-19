#!/usr/bin/env python3
"""Plot the saved M6Tw025 Li closure comparison for the thesis.

Run with Python 3, numpy, pandas and matplotlib. No solver is executed and no
analysis output is changed. The original study directory is the default;
--source-dir can point to an unchanged copy of its saved CSVs.

The two panels reproduce the middle row of run_comparison.py's
wall_flux_comparison.png. Both input-profile inversions use the DNS wall fluxes
as denominators, rather than each profile's own resolved wall values. The
horizontal coordinate is the common physical first-point height expressed in
DNS wall units. The y_1,DNS^+=200 point is retained from the original plot even
though its second reference point, y_2=3y_1, is outside both boundary layers.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "thesis-tm-closure-mpl")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_SOURCE = Path(
    "/home/mgnadlinger/OpenFOAM/mgnadlinger-v2412/src/LiWallFunctions/"
    "tests/ReverseSolver/M6Tw025_TM_closure_isolation"
)
DEFAULT_OUTPUT = Path(__file__).with_name("TM_prof_err_Li-aprio.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    data = pd.read_csv(args.source_dir / "comparison.csv")
    references = pd.read_csv(args.source_dir / "reference_states.csv", index_col="run")
    dns = references.loc["A_DNS"]
    heights = data["common_y1_DNS_plus"]
    np.testing.assert_allclose(heights, data["y1_m"] * dns["yplus_per_m"])
    np.testing.assert_allclose(data["y2_m"], 3.0 * data["y1_m"])

    # The shared denominator is essential to this comparison. Verify the
    # reconstruction against independently saved comparison columns.
    errors = {}
    for profile in ("A", "B"):
        for quantity, recovered, reference in (
            ("tau", "tau_w_recovered_Pa", "tau_Pa"),
            ("q", "q_w_recovered_W_per_m2", "q_W_m2"),
        ):
            error = 100.0 * (data[f"{profile}_{recovered}"] / dns[reference] - 1.0)
            saved = (
                f"A_{quantity}_w_error_percent"
                if profile == "A"
                else f"B_{quantity}_error_vs_DNS_percent"
            )
            np.testing.assert_allclose(error, data[saved], rtol=1e-10, atol=1e-10)
            errors[profile, quantity] = error

    plt.rcParams.update(
        {"font.size": 11, "axes.spines.top": False, "axes.spines.right": False}
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), dpi=180)
    fig.subplots_adjust(left=0.07, right=0.995, bottom=0.19, top=0.89, wspace=0.15)

    for profile, color, label in (
        ("A", "#222222", "DNS profile"),
        ("B", "#dc6803", r"Wall-resolved SA, $Pr_t=0.9$"),
    ):
        for ax, quantity in zip(axes, ("tau", "q")):
            ax.plot(
                heights, errors[profile, quantity], color=color,
                label=label, linewidth=2, marker="o", markersize=np.sqrt(20),
            )

    for ax, ylabel in zip(
        axes,
        ("Shear error vs DNS reference [%]", "Heat-flux error vs DNS reference [%]"),
    ):
        ax.set_xscale("log")
        ax.set_xlim(1, 200)
        ax.set_xlabel(r"Common sampling height $y_{1,\mathrm{DNS}}^+$")
        ax.set_ylabel(ylabel)
        ax.grid(True, which="both", alpha=0.18)
        ax.axhline(0, color="#888888", linewidth=0.8)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.savefig(args.output, dpi=180)
    plt.close(fig)
    print(f"Saved {args.output} (2160 x 774 pixels); all saved error checks passed.")


if __name__ == "__main__":
    main()
