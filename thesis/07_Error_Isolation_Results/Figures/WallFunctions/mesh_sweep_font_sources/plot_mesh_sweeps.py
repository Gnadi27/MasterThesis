#!/usr/bin/env python3
"""Replot the two 3% mesh sweeps at 170 mm wide, with 10 pt figure text.

Run normally to reproduce the PDF/PNG figures from retained plotting_data.json.
Use --refresh-data only when the original WSL report inputs are available.
No simulation, mesh, interpolation, envelope or normalization is recomputed.
"""

from pathlib import Path
import argparse
import csv
import hashlib
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/thesis-mesh-sweep-font-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import LogLocator
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
FONT_ROOT = Path("/mnt/c/texlive/2026/texmf-dist/fonts/opentype/public/tex-gyre")
MM = 1 / 25.4
WIDTH, HEIGHT = 170 * MM, 120 * MM


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def collect_data():
    """Freeze the original report inputs and exact plotted envelope samples."""
    import pandas as pd

    roots = {
        "LowHM3": Path("/home/mgnadlinger/WF_validation/LowH_M3/opt_meshes/marching_design"),
        "M6Tw025": Path("/home/mgnadlinger/Xue_et_al/flat_plate/M6Tw025/opt_mesh"),
    }
    result = {}
    provenance = {
        "purpose": "Typography-only replots of the existing cropped mesh-sweep figures.",
        "physical_size_mm": [170, 120],
        "font_size_pt": 10,
        "source_files_sha256": {},
        "scientific_operations": "The original complete faces_plus arrays and target_yplus values are retained. The left-panel centre separation remains (diff(f)[:-1]+diff(f)[1:])/2. The same delta99 mask, axis limits, full case names, and original envelope samples are used. No mesh, profile, coefficient, wall normalization or acceptance criterion is changed.",
        "layout_changes": "A shared full-size legend replaces the small left-panel inset legend; axis labels wrap to two lines. Report titles and footnotes remain omitted, as in the pre-existing manuscript crops.",
    }
    for name, root in roots.items():
        lowh = name == "LowHM3"
        search = root if lowh else root / "searches"
        meta_path = root / "reference/normalization.json"
        meta = json.loads(meta_path.read_text())
        csv_path = root / ("mesh_report_cubic_spacing_limits.csv" if lowh else "spacing_limits/symmetric_finite_stencil_limits.csv")
        script = root / ("plot_mesh_report.py" if lowh else "scripts/plot_analysis.py")
        sources = [meta_path, csv_path, script]
        records = []
        for path in (search / "3pct").glob("*_grid.json"):
            source = json.loads(path.read_text())
            records.append({key: source[key] for key in ("case", "target_yplus", "faces_plus")})
            sources.append(path)
        records.sort(key=lambda item: item["target_yplus"])
        assert len(records) == 12
        if lowh:
            rows = [r for r in csv.DictReader(csv_path.open()) if int(r["tolerance_percent"]) == 3]
            envelope_x = [float(r["y_face_plus"]) for r in rows]
            envelope_y = [float(r["h_plus"]) for r in rows]
        else:
            # This is exactly the CSV load/filter/grouping in plot_analysis.py.
            limits = pd.read_csv(csv_path)
            cubic = limits[(limits.tolerance_percent == 3) & (limits.interpolant == "cubic")].groupby("y_face_plus").h_plus.min()
            envelope_x, envelope_y = cubic.index.tolist(), cubic.values.tolist()
        edge = meta["delta99_plus"]
        result[name] = {
            "case": "LowH_M3" if lowh else name,
            "wall_units": "WR RANS" if lowh else "DNS",
            "delta99_plus": edge,
            "records": records,
            "envelope_y_face_plus": envelope_x,
            "envelope_h_plus": envelope_y,
            "geometry_y_face_plus": np.geomspace(.1, edge, 600 if lowh else 300).tolist(),
        }
        for path in sources:
            provenance["source_files_sha256"][str(path)] = sha256(path)
        print(f"Retained {name}: {len(records)} meshes and {len(envelope_x)} envelope points.")
    (HERE / "plotting_data.json").write_text(json.dumps(result, indent=2) + "\n")
    provenance["retained_plotting_data_sha256"] = sha256(HERE / "plotting_data.json")
    (HERE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def configure_font():
    candidates = [FONT_ROOT, Path("C:/texlive/2026/texmf-dist/fonts/opentype/public/tex-gyre")]
    for directory in candidates:
        if (directory / "texgyreheros-regular.otf").exists():
            for path in directory.glob("texgyreheros-*.otf"):
                font_manager.fontManager.addfont(str(path))
            break
    # findfont raises if neither the retained installation nor a system font exists.
    font_manager.findfont("TeX Gyre Heros", fallback_to_default=False)
    plt.rcParams.update({
        "font.family": "TeX Gyre Heros", "font.size": 10,
        "axes.labelsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.fontsize": 10, "mathtext.fontset": "cm",
        "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6,
        "xtick.minor.width": .45, "ytick.minor.width": .45,
    })


def plot_one(name, data):
    fig = plt.figure(figsize=(WIDTH, HEIGHT))
    # Figure coordinates in mm: the inter-panel gap accommodates complete case names.
    ax = fig.add_axes([16 / 170, 20 / 120, 59 / 170, 74 / 120])
    mesh = fig.add_axes([108 / 170, 20 / 120, 57 / 170, 74 / 120])
    edge = data["delta99_plus"]
    records = data["records"]
    for i, record in enumerate(records):
        f = np.asarray(record["faces_plus"])
        w = np.diff(f)
        mask = f[1:-1] <= edge + 1e-8
        ax.plot(f[1:-1][mask], ((w[:-1] + w[1:]) / 2)[mask],
                ".-", ms=1.8, lw=.6, label=record["case"])
        b = f[f <= edge + 1e-8]
        mesh.vlines(b[1:], i - .3, i + .3, color="#75889C", lw=.5)
        mesh.plot(record["target_yplus"], i, "o", color="#BA5A30", ms=3.2)
        # Validate all displayed mesh coordinates against retained source arrays.
        np.testing.assert_array_equal(ax.lines[-1].get_xdata(), f[1:-1][mask])
        np.testing.assert_array_equal(ax.lines[-1].get_ydata(), ((w[:-1] + w[1:]) / 2)[mask])
        np.testing.assert_array_equal([segment[0, 0] for segment in mesh.collections[-1].get_segments()], b[1:])
        assert mesh.lines[-1].get_xdata()[0] == record["target_yplus"]
    ax.plot(data["envelope_y_face_plus"], data["envelope_h_plus"],
            "--", color="#202C37", lw=.85, label="Cubic + geometry envelope")
    q = data["geometry_y_face_plus"]
    ax.plot(q, q, ":", color="#707D87", lw=.8, label=r"Geometry: $h^+=y_f^+$")
    units = data["wall_units"]
    ax.set(xscale="log", yscale="log", xlim=(.1, edge), ylim=(.005, 2000),
           xlabel=rf"Face location $y_f^+$" + f"\n({units} wall units)",
           ylabel=r"Full centre separation $h^+$")
    ax.grid(True, which="both", alpha=.23, linewidth=.4)
    mesh.set(xscale="log", xlim=(.05, edge * 1.15),
             yticks=range(len(records)), yticklabels=[r["case"] for r in records],
             xlabel=r"Face or sampling location $y^+$" + f"\n({units} wall units)")
    mesh.grid(axis="x", alpha=.25, linewidth=.4)
    # Match the source autoscale exactly for [-.3, 11.3] face bars.
    mesh.set_ylim(-.88, 11.88)
    for axis in (ax, mesh):
        axis.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
        axis.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
        axis.tick_params(axis="both", which="major", pad=3, length=3)
        axis.tick_params(axis="both", which="minor", length=1.7)
        axis.xaxis.labelpad = 4
    ax.yaxis.labelpad = 4
    lines = ax.lines
    legend = fig.legend(lines[:12], [r["case"] for r in records],
                        loc="upper center", bbox_to_anchor=(.5, .995), ncol=4,
                        frameon=False, handlelength=1.6, columnspacing=1.5,
                        handletextpad=.5, labelspacing=.25, borderaxespad=0)
    envelope_legend = fig.legend(lines[12:], [line.get_label() for line in lines[12:]],
                                 loc="upper center", bbox_to_anchor=(.52, .86), ncol=2,
                                 frameon=False, handlelength=2.3, columnspacing=1.4,
                                 handletextpad=.5, borderaxespad=0)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    frame = fig.bbox
    # get_tightbbox includes rendered tick labels, unlike the unfiltered Text
    # children, which also include out-of-view logarithmic locator ticks.
    bb = fig.get_tightbbox(renderer).transformed(fig.dpi_scale_trans)
    if bb.x0 < -.5 or bb.y0 < -.5 or bb.x1 > frame.x1 + .5 or bb.y1 > frame.y1 + .5:
        raise ValueError(f"Rendered content outside figure: {bb}")
    assert legend.get_window_extent(renderer).y0 > envelope_legend.get_window_extent(renderer).y1
    assert envelope_legend.get_window_extent(renderer).y0 > ax.bbox.y1
    output = OUT / f"{name}_mesh_sweep_3pct"
    fig.savefig(output.with_suffix(".pdf"))
    fig.savefig(output.with_suffix(".png"), dpi=300)
    print(f"Wrote {output.name}.pdf/.png: 170 × 120 mm; all figure text 10 pt.")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-data", action="store_true")
    args = parser.parse_args()
    if args.refresh_data:
        collect_data()
    configure_font()
    data = json.loads((HERE / "plotting_data.json").read_text())
    provenance = json.loads((HERE / "provenance.json").read_text())
    assert sha256(HERE / "plotting_data.json") == provenance["retained_plotting_data_sha256"]
    for name, record in data.items():
        plot_one(name, record)


if __name__ == "__main__":
    main()
