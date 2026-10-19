"""Replot the Cao low-Re Stanton figure with capital Delta spacing labels.

Uses validate_cao.py and its retained full_sweep_output CSVs from
/home/mgnadlinger/Cao_lamZPGFP/validation. The preceding manuscript PNG was
pixel-identical to rows 772:1454 of that original 1408x1454, 200-dpi export.
Keep this crop, styling, and numerical data; change only the legend notation.
Requires NumPy, Matplotlib, and LaTeX (lmodern/dvipng) on PATH.
"""

from pathlib import Path
from types import SimpleNamespace
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.transforms import Bbox
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def geometry_hash(figure):
    arrays = []
    for axis in figure.axes:
        arrays.extend(np.column_stack(line.get_data()) for line in axis.lines)
        arrays.extend(path.vertices for collection in axis.collections
                      for path in collection.get_paths())
    return hashlib.sha256(b"".join(array.astype("<f8").tobytes()
                                   for array in arrays)).hexdigest()


def main():
    default_source = (Path(r"\\wsl.localhost\Ubuntu\home\mgnadlinger\Cao_lamZPGFP\validation")
                      if os.name == "nt" else Path("/home/mgnadlinger/Cao_lamZPGFP/validation"))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", type=Path, default=default_source)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = args.validation_dir / "validate_cao.py"
    spec = importlib.util.spec_from_file_location("cao_validation", source)
    base = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = base
    spec.loader.exec_module(base)
    data = args.validation_dir / "full_sweep_output"
    wall_path = data / "760_pinf_dy1_1e-3_dx_1p5_wall.csv"
    reference_path = data / "760_pinf_reference_stanton.csv"
    wall = read_csv(wall_path)
    reference = {"stanton": read_csv(reference_path)}
    result = SimpleNamespace(spec=SimpleNamespace(label="dy1_1e-3_dx_1p5"), wall=wall)
    audit = {"original_script": str(source), "script_sha256": sha(source),
             "matplotlib_version": matplotlib.__version__,
             "inputs": [{"path": str(path), "sha256": sha(path)}
                        for path in (wall_path, reference_path)],
             "original_crop_pixels": [0, 772, 1408, 1454]}
    original_savefig = Figure.savefig

    def export(figure, destination, *unused, **kwargs):
        figure.set_dpi(200)
        figure.canvas.draw()
        box = figure.get_tightbbox(figure.canvas.get_renderer()).padded(
            plt.rcParams["savefig.pad_inches"])
        # Anchor at the original export's bottom: its raster height is rounded
        # down. This preserves the pixel alignment of its final 682 rows.
        # The epsilon avoids truncating 681.9999999999999 pixels to 681.
        crop = Bbox.from_bounds(box.x0, box.y0, 1408 / 200 + 1e-8,
                               682 / 200 + 1e-8)
        before = geometry_hash(figure)
        legend = figure.axes[1].get_legend()
        labels = []
        for text in legend.get_texts():
            old = text.get_text()
            new = old.replace(r"\delta y_1", r"\Delta y_1").replace(r"\delta x", r"\Delta x")
            text.set_text(new)
            labels.append(new)
        assert r"\Delta y_1" in labels[0] and r"\Delta x" in labels[0]
        assert geometry_hash(figure) == before
        audit.update(legend=labels, curve_and_band_sha256=before,
                     numerical_geometry_unchanged=True)
        for extension in ("png", "pdf"):
            path = Path(destination).with_suffix("." + extension)
            temporary = path.with_name(path.name + ".new")
            original_savefig(figure, temporary, format=extension, dpi=200,
                             bbox_inches=crop, pad_inches=0)
            temporary.replace(path)

    plt.rcdefaults()
    base.configure_plot_style()
    Figure.savefig = export
    try:
        base.plot_wall_coefficients(base.reference_specs(base.REPO_ROOT)["760_pinf"],
                                    [result], reference,
                                    args.output_dir / "05_760_pinf_wall_coefficients.png", 200)
    finally:
        Figure.savefig = original_savefig
    (args.output_dir / "cao_low_re_stanton_provenance.json").write_text(
        json.dumps(audit, indent=2) + "\n")
    print("Saved Cao low-Re Stanton figure with capital Delta y_1 and Delta x.")


if __name__ == "__main__":
    main()
