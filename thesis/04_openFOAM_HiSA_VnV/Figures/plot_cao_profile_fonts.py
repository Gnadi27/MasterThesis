"""Reproduce Figure 4.29's two Cao profiles with enlarged DejaVu Sans text.

Run with Python, NumPy and Matplotlib. Source CSV snapshots and their original
WSL paths are recorded in cao_profile_sources/source_provenance.json.
The snapshots preserve the currently included legacy dy1=0.001 mm, dx=1.5 mm
cases at saved pseudo-iteration 50000 for BOTH pressures. The older manuscript
comment about a finest-mesh 1550 Pa export does not describe the retained PNG:
that PNG exactly matched rows 157:1215 of the selected_legacy_1550 export.

Only typography/layout change. All plotted coordinate hashes are compared
against captures of the original generating scripts before exporting. The
1550 Pa plot retains no title; the 760 Pa plot retains its original title.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text
import numpy as np


HERE = Path(__file__).resolve().parent
SOURCES = HERE / "cao_profile_sources"
SCALE = 1.25
DPI = 200


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(record, key):
    source = record["snapshots"][key]
    path = SOURCES / source["file"]
    if sha(path) != source["sha256"]:
        raise ValueError(f"Source CSV changed: {path}")
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def line_hashes(axes):
    return [hashlib.sha256(np.column_stack((line.get_xdata(), line.get_ydata()))
                          .astype("<f8").tobytes()).hexdigest()
            for ax in axes for line in ax.lines]


def verify_text(fig):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ignored = set()
    for ax in fig.axes:
        for axis, bounds in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if not bounds[0] <= tick.get_loc() <= bounds[1]:
                    ignored.update((id(tick.label1), id(tick.label2)))
    for text in fig.findobj(Text):
        if id(text) in ignored or not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < -.5 or box.y0 < -.5 or box.x1 > fig.bbox.width + .5 or box.y1 > fig.bbox.height + .5:
            raise ValueError(f"Text outside figure: {text.get_text()}")
    for ax in fig.axes:
        legend_box = ax.get_legend().get_window_extent(renderer).padded(3)
        for line in ax.lines:
            if line.get_label().startswith("_"):
                continue  # The dotted zero-velocity guide is not a profile.
            if line.get_path().transformed(line.get_transform()).intersects_bbox(legend_box, filled=False):
                raise ValueError(f"Legend overlaps plotted line: {line.get_label()}")


def plot(record, output):
    pressure = record["pressure_Pa"]
    d = read_csv(record, "source_profile_csv")
    u = read_csv(record, "reference_velocity_processed_csv")
    t = read_csv(record, "reference_temperature_processed_csv")
    mask = d["eta"] <= 18
    # Fixed 1823-pixel output width equals the preceding PNGs. Thus the 1.25
    # point-size factor remains exactly 1.25 at the same manuscript linewidth.
    height = 1058 if pressure == 1550 else 1215
    fig, axes = plt.subplots(1, 2, figsize=(1823 / DPI, height / DPI),
                             sharey=True, dpi=DPI)
    cfd_label = ("OpenFOAM/HiSA, " + r"$x/L=0.4$" + "\n"
                 + r"$\Delta y_1=0.001\,\mathrm{mm},\;\Delta x=1.5\,\mathrm{mm}$")
    axes[0].plot(d["U_over_Uinf"][mask], d["eta"][mask], color="tab:blue", label=cfd_label)
    axes[1].plot(d["T_over_Tinf"][mask], d["eta"][mask], color="tab:blue", label=cfd_label)
    style = dict(color="black", linestyle="-", linewidth=1.15, marker="o",
                 markersize=3.4, markerfacecolor="white", markeredgewidth=.75,
                 label="Cao et al. (digitized)", zorder=6)
    axes[0].plot(u["U/Uinf"], u["eta"], **style)
    axes[1].plot(t["T/Tinf"], t["eta"], **style)
    axes[0].axvline(0, color="0.45", linestyle=":", linewidth=1)
    for ax, original, legend_location in zip(axes, record["original_axes"], ("upper left", "upper right")):
        ax.set_xlabel(original["xlabel"])
        ax.set_ylabel(original["ylabel"])
        ax.set_xticks(original["xticks"])
        ax.set_yticks(original["yticks"])
        ax.set_xlim(original["xlim"])
        ax.set_ylim(original["ylim"])
        ax.minorticks_on()
        ax.legend(loc=legend_location)
    if record["previous_title_visible"]:
        fig.suptitle("Flat-plate profile comparison\n"
                     + rf"$p_\infty={pressure}\,\mathrm{{Pa}},\quad Re_\infty\approx 4.2\times10^6\,\mathrm{{m^{{-1}}}}$")
    fig.tight_layout(pad=.55, w_pad=.55)
    hashes = line_hashes(axes)
    if hashes != [line["sha256"] for line in record["original_lines"]]:
        raise ValueError(f"Original plotted values changed: {pressure} Pa")
    verify_text(fig)
    outputs = []
    for extension in ("png", "pdf"):
        path = output / f"05_{pressure}_pinf_boundary_layer_profiles.{extension}"
        # Replace atomically rather than writing through a historical hardlink.
        temporary = path.with_name(path.name + ".new")
        fig.savefig(temporary, format=extension, dpi=DPI, bbox_inches=None)
        temporary.replace(path)
        outputs.append({"path": str(path), "sha256": sha(path)})
    plt.close(fig)
    return {"pressure_Pa": pressure, "source_provenance": str(SOURCES / "source_provenance.json"),
            "profile_snapshot_pseudo_iteration": record["profile_snapshot_pseudo_iteration"],
            "line_sha256": hashes, "matches_original_plotted_values": True,
            "axes": record["original_axes"], "title_visible": record["previous_title_visible"],
            "font_scale_at_manuscript_width": SCALE, "output_pixels": [1823, height],
            "text_within_canvas": True, "legends_clear_of_plotted_lines": True, "outputs": outputs}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcdefaults()
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans", "text.usetex": False,
        "font.size": 12 * SCALE, "axes.labelsize": 13 * SCALE,
        "legend.fontsize": 9.5 * SCALE, "figure.titlesize": 14 * SCALE,
        "xtick.labelsize": 11 * SCALE, "ytick.labelsize": 11 * SCALE,
        "lines.linewidth": 1.8, "axes.grid": True, "grid.alpha": .24, "grid.linewidth": .7,
        "savefig.facecolor": "white", "figure.facecolor": "white", "axes.linewidth": .8,
        "pdf.fonttype": 42,
    })
    records = json.loads((SOURCES / "source_provenance.json").read_text())["records"]
    report = {"font_family": "DejaVu Sans", "mathtext_fontset": "dejavusans",
              "font_scale": SCALE, "matplotlib_version": matplotlib.__version__,
              "records": [plot(record, args.output_dir.resolve()) for record in records]}
    target = args.output_dir / "cao_profile_sources/font_update_provenance.json"
    target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(report, indent=2) + "\n")
    print("Saved both Cao profile figures; all original plotted coordinates are unchanged.")


if __name__ == "__main__":
    main()
