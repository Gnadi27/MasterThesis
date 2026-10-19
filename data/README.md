# Thesis data archives

The archives in this directory use Git LFS. Install Git LFS before cloning,
or run `git lfs install` and `git lfs pull` in an existing clone.

`01_VnV.7z`, `02_WallResolvedTMs.7z` and `04_WallFunctions.7z` are stored as
numbered volumes of at most 1,900,000,000 bytes each. Download every volume for
an archive into the same directory, then open its `.7z.001` file with 7-Zip to
extract it. The volumes
concatenate byte for byte to the original archive; no recompression was used.
On Linux or WSL, for example:

```sh
cat 01_VnV.7z.??? > 01_VnV.7z
cat 02_WallResolvedTMs.7z.??? > 02_WallResolvedTMs.7z
cat 04_WallFunctions.7z.??? > 04_WallFunctions.7z
sha256sum -c SHA256SUMS
```

Run these commands from `data/`. `SHA256SUMS` contains hashes for the stored
files and all three reconstructed complete archives. Complete local copies of
the split archives are ignored by Git so that they are not uploaded twice.

# Thesis support archive

Simulation cases, saved results, analysis material and custom code supporting Michael Roman Gnadlinger's thesis, **Assessment of RANS Wall Functions for Predicting Aerothermal Heating in Hypersonic Cold-Wall Flows**.

The manuscript is maintained separately in `Docs_and_Presentations/Thesis_tex/thesis/`, with `main.tex` defining the chapter order. The guide below follows that structure.

## Where to find things

| Folder | Contents | Connection to the thesis |
| --- | --- | --- |
| [OpenFOAM_customLibs/](OpenFOAM_customLibs/) | Modified HiSA solver; Li, Griffin and Wang wall functions; thermophysical models; density-corrected Spalart–Allmaras model. | **Chapter 3: Numerical framework**, **Chapter 5: Wall-function formulations and implementation**, and the SA density-treatment appendix. |
| [VnV/](VnV/) | Shock, laminar-flow and adiabatic turbulent-flow validation cases, with reference data and supporting material. | **Chapter 4: Verification and validation of the numerical framework**; also the Malik–Anderson Sutherland-coefficient appendix. |
| [WallResolvedTMs/](WallResolvedTMs/) | Wall-resolved SA/SST benchmark cases, mesh studies and comparison data. | **Chapter 6: Benchmark cases and evaluation methodology**; **Chapter 7, wall-resolved turbulence-model results**; mesh-convergence and additional-profile appendices. |
| [WallFunctionClosures/A-Priori-Solver/](WallFunctionClosures/A-Priori-Solver/) | Li and Griffin closure tests using supplied boundary-layer profiles. | **Chapter 7, wall-function closure assessment**, and the **a-priori wall-function appendix**. |
| [WallFunctions/](WallFunctions/) | Coupled CFD cases on conventional and boundary-layer-informed meshes. | **Chapter 7, spatial-discretization error and sampling sensitivity**, using the benchmark and mesh definitions in Chapter 6. |

Within **WallResolvedTMs**, `TM_12-ZPGs_wall_resolved_baselines/` contains the benchmark baselines; `TM_12-ZPGs_mesh_convergence/` and `M6Tw025_2D_mesh_convergence/` contain mesh studies.

Within **WallFunctions**, start with [README.txt](WallFunctions/README.txt), then choose `LowH_M3/` (the nominal Mach-3.5 case) or `M6Tw025/`. Each contains `Conventional/` and `BoundaryLayerInformed/` collections, with run indexes for locating individual cases.

## Chapter 4 validation cases

The numbered folders inside **VnV** follow the order of the chapter's sections:

1. `01_Shock_Angle_Post_Shock/` — wedge shock angles and post-shock states.
2. `02_Shock_Stand_Off_Carbuncle/` — Kim/Hornung shock stand-off and carbuncle studies.
3. `03_Laminar_Boundary_Layer_Profiles/` — Malik–Anderson and Marxen profiles; the separate [Sutherland-coefficient study](VnV/03_Laminar_Boundary_Layer_Profiles/MalikAnderson/Appendix_Sutherland_Coefficient_Study/) supports the appendix.
4. `04_Laminar_Heat_Transfer_Compression_Corner/` — compression-corner heating comparisons.
5. `05_Adiabatic_Turbulent_Boundary_Layer_Profiles/` — Coles SA/SST cases, including their thermophysical support package.

## Reading and reusing the archive

Prefer the extracted folders for browsing: `.7z` packages may predate recovered files and later additions. Local READMEs, run indexes and copy manifests provide more detail and original source locations.

In a typical CFD case, `system/` contains run settings, `constant/` contains model settings and the mesh, `0/` or `0.org/` holds initial fields or templates, numbered folders hold saved solution states, and `postProcessing/` contains extracted results.

Some runs are incomplete, and selected VnV cases still lack original setup files. See the [VnV recovery summary](VnV/Recovery_summary_2026-10-19.txt) for the remaining gaps. Historical launch paths and software dependencies may need adaptation; this guide is a navigation aid, not a universal run recipe.
