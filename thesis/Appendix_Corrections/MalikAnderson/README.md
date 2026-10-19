# Malik--Anderson Sutherland-coefficient appendix

The thesis includes `Response_to_Cosimo_d1_cancellation.tex`. The exact original
report, dated 31 August 2026, is preserved in `source_as_supplied/`.

Source directory:
`../../Files/ValData/FIN/MalikAnderson/` (relative to the `Thesis_tex` workspace).

Source SHA-256: `88407bf69e9bc65e25a7ca2cd4ba8cf989692a6795a057003412f7b630f2f3b0`

The adaptation removes the standalone document setup and bibliography, adds
the thesis appendix chapter, scopes and renames report-specific commands and
colours, namespaces internal labels with `app:ma:`, replaces the report's
internal appendix reset with continuing section numbering, and maps citations
to the thesis bibliography. Layout-only changes wrap a long equation and the
opening DOI, and make the original boxes fit the thesis text width. The coefficient tests,
derivations, numerical results, and caveats are retained. A missing semicolon
in the marching-diagnostic paragraph is corrected.

The tested coefficient 1.5869e-6 remains a candidate, not an established unique
correction to the printed 1.0869e-6. The failed printed-coefficient march remains
explicitly identified as a non-converged diagnostic.

Figure copies (all unmodified):
- `01 Ma10_prfG_x0p1/MA1991_Fig1_Ma10_perfect_profile_validation.png` -> `Figures/Ma10_perfect.png`
- `02 Ma15_prfG_x0p1/MA1991_Fig3_Ma15_perfect_profile_validation.png` -> `Figures/Ma15_perfect.png`
- `03 Ma15_impG_x0p1/MA1991_Fig3_Ma15_imperfect_profile_validation_with_CN_diagnostic.png` -> `Figures/Ma15_imperfect_CN.png`

Main-document dependencies:
- `\include{Appendix_Corrections/MalikAnderson/Response_to_Cosimo_d1_cancellation}`
  after `\appendix`;
- bibliography database `XXX_Bibliography/Appendix_MalikAnderson`;
- existing `MalikAnderson1991RealGas` entry;
- packages already supplied by `tumphthesis`: `amsmath`, `graphicx`, `booktabs`,
  `xcolor`, `hyperref`, and `natbib`; `xurl` is already loaded by `main.tex`.

The standalone source is archival and is not included in the thesis build.
