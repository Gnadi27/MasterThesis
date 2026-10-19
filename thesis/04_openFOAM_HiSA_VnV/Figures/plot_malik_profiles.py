"""Update the typography of the three Malik--Anderson manuscript figures.

Run from any directory with Python, numpy, pandas, matplotlib and Pillow.
The default data root is the adjacent Files/ValData archive; override it with
--data-root. Mach-15 curves are read through the original normalization
function and reference arrays. Original source scripts/CSVs are not modified.

Mach-10 source CSVs are unavailable. Its retained original PNG is therefore
used ONLY as a raster layer: no curves are digitized, fitted or substituted.
Only the two empty-background legend regions are cleared; all remaining
plot-area pixels are retained. Enlarged labels and legends are newly typeset.
The user explicitly selected this typography-only treatment on 2026-10-19.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
from PIL import Image


HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[2]
SCALE = 1.25
NAMES = {
    "ma10": "04_MalikAnderson_Ma10_perfGas_lamBL_profiles",
    "ma15_perfect": "04_MalikAnderson_Ma15_perfGas_lamBL_profiles",
    "ma15_imperfect": "04_MalikAnderson_Ma15_imperfGas_lamBL_profiles",
}
SOURCES = {
    "ma15_perfect": "Tunel/00 LamZPGFP/OLDMalikAnderson/Ma15/01 perfectGas",
    "ma15_imperfect": "FIN/MalikAnderson/Ma15_impG/03 imperfectGas(newSuth,LEshock)",
}
COLORS = plt.get_cmap("viridis")(np.linspace(0.15, 0.85, 5))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fonts(base, title, label, tick, legend):
    plt.rcdefaults()
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
        "text.usetex": False, "font.size": base,
        "axes.titlesize": title, "axes.labelsize": label,
        "xtick.labelsize": tick, "ytick.labelsize": tick,
        "legend.fontsize": legend, "axes.linewidth": 0.8,
        "savefig.facecolor": "white", "figure.facecolor": "white",
        "pdf.fonttype": 42,
    })


def decorate(axes, ymax):
    for i, ax in enumerate(axes):
        ax.set_xlim(0, 1.05 if i == 0 else 19)
        ax.set_ylim(0, ymax)
        ax.set_xticks(np.arange(0, 1.01, .2) if i == 0 else np.arange(0, 19, 2.5))
        ax.set_yticks(np.arange(0, ymax + 1, 5))
        ax.set_xlabel(r"$U / u_e$" if i == 0 else r"$T / T_e$")
        ax.set_ylabel(r"$Y = y / L_r(x),\quad L_r = \sqrt{\nu_e\, x / u_e}$"
                      if i == 0 else r"$Y = y / L_r(x)$")
        ax.set_title("(a) Velocity profile" if i == 0 else "(b) Temperature profile")


def save(fig, key, output, dpi):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for txt in fig.findobj(Text):
        if txt.get_visible() and txt.get_text():
            b = txt.get_window_extent(renderer)
            if b.x0 < -.5 or b.y0 < -.5 or b.x1 > fig.bbox.width + .5 or b.y1 > fig.bbox.height + .5:
                raise ValueError(f"Text outside figure: {key}: {txt.get_text()}")
    for extension in ("png", "pdf"):
        fig.savefig(output / (NAMES[key] + "." + extension), dpi=dpi)


def original_reader(script):
    """Load only scientific constants, reference arrays and the original reader."""
    tree = ast.parse(script.read_text(encoding="utf-8"), filename=str(script))
    names = {"T_E", "U_E", "RE_UNIT", "M_E", "ref_U", "ref_T", "ref_Y_u", "ref_Y_t"}
    nodes = []
    for node in tree.body:
        if isinstance(node, ast.Import) and all(a.name in {"numpy", "pandas"} for a in node.names):
            nodes.append(node)
        elif isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id in names for t in node.targets):
            nodes.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name == "load_and_normalize":
            nodes.append(node)
    namespace = {}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(script), "exec"), namespace)
    return namespace


def line_hashes(axes):
    return [hashlib.sha256(np.column_stack((line.get_xdata(), line.get_ydata()))
                          .astype("<f8").tobytes()).hexdigest()
            for ax in axes for line in ax.lines]


def ma15(key, data_root, output):
    script = data_root / SOURCES[key] / "malik_anderson_plot.py"
    source = original_reader(script)
    csvs = [script.parent / f"x0p{s}.csv" for s in (1, 2, 4, 6, 8)]
    stations = [source["load_and_normalize"](p) for p in csvs]
    # Keep the source sizes and multiply every text category, including the
    # explicit 8.5 pt legend override, by 1.25. Keep the original 12-inch width.
    fonts(12*SCALE, 14*SCALE, 13*SCALE, 11*SCALE, 8.5*SCALE)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=180)
    for s, color in zip(stations, COLORS):
        label = f"x = {s['x']:.3f} m  ($\\bar\\chi$ = {s['chi_bar']:.2f})"
        for ax, quantity in zip(axes, ("U_star", "T_star")):
            ax.plot(s[quantity], s["Y_star"], lw=2, color=color, label=label)
    for ax, x, y in zip(axes, ("ref_U", "ref_T"), ("ref_Y_u", "ref_Y_t")):
        ax.plot(source[x], source[y], color="black", lw=1, marker="o", ms=5,
                mfc="white", label="Malik & Anderson 1991")
        ax.grid(True, alpha=.3)
    decorate(axes, 30)
    for ax, loc in zip(axes, ("upper left", "upper right")):
        ax.legend(loc=loc, framealpha=.95)
    fig.tight_layout(pad=.7, w_pad=.8)
    hashes = line_hashes(axes)
    historical = WORKSPACE / "tmp/laminar_font_review/render_audit.json"
    checked = False
    if historical.exists():
        old = json.loads(historical.read_text())[key]["exports"][0]["lines"]
        assert hashes == [line["sha256"] for line in old], "Original plotted values changed"
        checked = True
    save(fig, key, output, 180)
    plt.close(fig)
    return {"method": "replot from original CSVs and digitized reference arrays",
            "script": str(script), "script_sha256": sha(script),
            "csvs": [{"path": str(p), "sha256": sha(p)} for p in csvs],
            "line_sha256": hashes, "matches_previous_render_audit": checked,
            "font_scale": SCALE}


def clear_legend(rgb, box, quiet_columns, quiet_rows):
    """Restore only grid/background under an old legend, never curve pixels.

    The selected rectangles contain legends and empty plot background in the
    retained original. Median samples from quiet rows/columns in the SAME
    panel recover its horizontal/vertical grid without estimating any data.
    """
    x0, y0, x1, y1 = box
    horizontal = np.median(rgb[y0:y1, quiet_columns[0]:quiet_columns[1]], axis=1)
    vertical = np.median(rgb[quiet_rows[0]:quiet_rows[1], x0:x1], axis=0)
    rgb[y0:y1, x0:x1] = np.rint(horizontal[:, None] * vertical[None, :] / 255).astype(np.uint8)


def ma10(output):
    original = HERE / "malik_plot_sources/ma10_original.png"
    rgb = np.asarray(Image.open(original).convert("RGB")).copy()
    untouched = rgb.copy()
    # Original legend rectangles, confirmed visually to contain no curves.
    boxes = [(426, 290, 630, 421), (1062, 32, 1265, 164)]
    clear_legend(rgb, boxes[0], (80, 380), (200, 280))
    clear_legend(rgb, boxes[1], (900, 1050), (170, 220))
    allowed = np.zeros(rgb.shape[:2], dtype=bool)
    for x0, y0, x1, y1 in boxes:
        allowed[y0:y1, x0:x1] = True
        patch = rgb[y0:y1, x0:x1]
        assert np.array_equal(patch[..., 0], patch[..., 1])
        assert np.array_equal(patch[..., 1], patch[..., 2])
    assert np.array_equal(rgb[~allowed], untouched[~allowed])

    # Slightly wider margins accommodate larger labels while retaining both
    # original 566 x 397 pixel panels. Compensate typography for the wider
    # canvas when included at the same manuscript linewidth.
    width, height, dpi = 1324, 520, 100
    compensation = width / 1280
    scale = SCALE * compensation
    # Original raster's effective font sizes at 100 dpi: 13 pt titles,
    # 10.8 pt labels/ticks, 9.18 pt legends (matching DejaVu Sans glyphs).
    fonts(10.8*scale, 13*scale, 10.8*scale, 10.8*scale, 9.18*scale)
    fig = plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
    axes = [fig.add_axes([x/width, 67/height, 566/width, 397/height]) for x in (90, 744)]
    decorate(axes, 40)
    for ax, x0 in zip(axes, (68, 704)):
        panel = rgb[28:425, x0:x0+566]
        ax.imshow(panel, extent=(*ax.get_xlim(), 0, 40), origin="upper",
                  interpolation="none", aspect="auto", zorder=0)
        for spine in ax.spines.values():
            spine.set_visible(False)  # Original raster retains its borders.
    handles = [Line2D([], [], color=c, lw=2,
                    label=f"x = {x:.3f} m  ($\\bar\\chi$ = {chi:.2f})")
               for c, x, chi in zip(COLORS, (.05, .1, .2, .4, .606), (1.74, 1.23, .87, .62, .50))]
    handles.append(Line2D([], [], color="black", lw=1, marker="o", ms=5,
                          mfc="white", label="Malik & Anderson 1991"))
    for ax, loc in zip(axes, ("upper left", "upper right")):
        ax.legend(handles=handles, loc=loc, framealpha=.95)
    save(fig, "ma10", output, dpi)
    plt.close(fig)
    return {"method": "typography only; retained original curve raster, no recovered numerical data",
            "original": str(original), "original_sha256": sha(original),
            "original_data_status": "missing: x0p05m.csv, x0p1m.csv, x0p2m.csv, x0p4m.csv, x0_606m.csv",
            "legend_background_rectangles_xyxy": boxes,
            "pixels_outside_legend_rectangles_unchanged": True,
            "original_panel_pixels": [566, 397], "font_scale_at_manuscript_width": SCALE}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=WORKSPACE.parents[1] / "Files/ValData")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {"font_family": "DejaVu Sans", "math_font": "DejaVu Sans", "figures": {}}
    for key in SOURCES:
        report["figures"][key] = ma15(key, args.data_root.resolve(), args.output_dir)
    report["figures"]["ma10"] = ma10(args.output_dir)
    (args.output_dir / "malik_plot_sources/font_update_provenance.json").parent.mkdir(exist_ok=True)
    (args.output_dir / "malik_plot_sources/font_update_provenance.json").write_text(json.dumps(report, indent=2) + "\n")
    print("Saved all three figures; Mach-15 line values verified against the previous render audit.")


if __name__ == "__main__":
    main()
