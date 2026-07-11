# Gao Wall Functions - Design And Deviations Log

## Scope And Packaging

This implementation is a self-contained OpenFOAM v2412 user library:

- Source directory: `$WM_PROJECT_USER_DIR/src/GaoWallFunctions`
- Built library: `$FOAM_USER_LIBBIN/libGaoWallFunctions.so`
- Build command used: `source /opt/of-v2412/etc/bashrc && cd $WM_PROJECT_USER_DIR/src/GaoWallFunctions && wmake libso`
- No OpenFOAM core library or existing turbulence-library source was edited or rebuilt.

Runtime load entry:

```foam
libs
(
    "libGaoWallFunctions.so"
);
```

Boundary condition type names:

- `nut`: `nutGaoWallFunction`
- `alphat`: `compressible::alphatGaoWallFunction`
- `omega`: `omegaGaoWallFunction`
- `k`: `kGaoWallFunction`

The compiled environment was:

- `WM_PROJECT_DIR=/opt/of-v2412`
- `WM_PROJECT_USER_DIR=/home/mgnadlinger/OpenFOAM/mgnadlinger-v2412`
- `WM_OPTIONS=linux64GccDPInt32Opt`
- `FOAM_USER_LIBBIN=/home/mgnadlinger/OpenFOAM/mgnadlinger-v2412/platforms/linux64GccDPInt32Opt/lib`

## Make Files

`Make/files`:

```make
GaoWallFunction/GaoWallFunction.C
nutGaoWallFunction/nutGaoWallFunctionFvPatchScalarField.C
alphatGaoWallFunction/alphatGaoWallFunctionFvPatchScalarField.C
omegaGaoWallFunction/omegaGaoWallFunctionFvPatchScalarField.C
kGaoWallFunction/kGaoWallFunctionFvPatchScalarField.C

LIB = $(FOAM_USER_LIBBIN)/libGaoWallFunctions
```

`Make/options` links only the OpenFOAM libraries needed for finite-volume fields, thermophysical properties, and compressible turbulence:

```make
EXE_INC = \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/transportModels/compressible/lnInclude \
    -I$(LIB_SRC)/thermophysicalModels/basic/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
    -I$(LIB_SRC)/TurbulenceModels/compressible/lnInclude \
    -I$(LIB_SRC)/thermoTools/lnInclude \
    -IGaoWallFunction

LIB_LIBS = \
    -lfiniteVolume \
    -lmeshTools \
    -lcompressibleTransportModels \
    -lfluidThermophysicalModels \
    -lturbulenceModels \
    -lcompressibleTurbulenceModels \
    -lthermoTools
```

## Implemented Equations

The shared helper `GaoWallFunction::calcPatch` evaluates the wall-face coupled solve used by all four BCs.

Gao et al. is treated as the authority where it corrects Nichols and Nelson:

- `kappa=0.41`, `B=5.0` from Gao section 2.2.1, not Nichols' `0.4/5.5`.
- `beta = 1.2*qw*mu_w/(rho_w*Tw*k_w*uTau)` from Gao Eq. 16.
- `dyWhitePlus/dyPlus` uses the `-0.5` exponent from Gao Eq. 9.
- Modified Temperature Law-of-the-wall, specifically Gao Eq. 6 using the modified beta expression from Eq. 16.
- The energy wall flux uses Gao Eq. 18, `qHalf = qw - r*uHalf*tauw/Prt`.

Nichols and Nelson is used only where Gao omits the k-omega first-cell values:

- `omega_i = 6*mu_w/(0.075*rho_w*y^2)` from Nichols Eq. 19.
- `omega_o = uTau/(sqrt(Cmu)*kappa*y)` from Nichols Eq. 20.
- `omega = sqrt(omega_i^2 + omega_o^2)` from Nichols Eq. 21.
- `k = omega*mu_t/rho` from Nichols Eq. 22.

## Coupled Solve

For each wall face, the nonlinear velocity law is solved in `uPlus`. For a trial `uPlus`, the helper computes:

- `uTau = |U_cell - U_wall|/uPlus`
- `yPlus = rho_w*uTau*y/mu_w`
- `Gamma = r*uTau^2/(2*cp*Tw)`
- `beta` from Gao Eq. 6 for isothermal walls, or `beta=0` for adiabatic walls
- the Gao/White velocity law residual from Gao Eq. 5

Initial guess:

- `uTau = sqrt((nut + nu)*|snGrad(U)|)`, using the current OpenFOAM wall-field values.
- If that is too small but the slip-relative cell velocity is nonzero, the fallback is `sqrt(nu*|U|/y)`.
- This is converted to the initial `uPlus = |U_cell - U_wall|/uTau`.

Iteration:

- Bracketed scalar bisection on `uPlus`.
- The residual is `lawYPlus(uPlus) - rho_w*|U|*y/(mu_w*uPlus)`.
- Defaults: `maxIter=50`, `tolerance=1e-7`.
- Converged if relative bracket width or absolute velocity-law residual is below tolerance.

The original implementation used damped Newton iteration in `uTau`. That was too fragile for coarse-wall-cell starts: on a yPlus-about-100 mesh it could miss the physical `uPlus` root, land on a very large `uPlus`, and feed an enormous Eq. 8 eddy viscosity into Nichols Eq. 22 for `k`. The bracketed `uPlus` solve was added to remove that implementation failure mode.

## Thermal Mode

The dictionary entry `thermalMode` can be:

- `auto` (default): fixed-value `T` patches are treated as isothermal; non-fixed-value `T` patches are treated as adiabatic.
- `isothermal`: use the current wall `T` patch value as `Tw`.
- `adiabatic`: set `qw=0` and compute `Tw = T_cell + r*u^2/(2*cp)` from Gao Eq. 7.

For wall density, pressure is taken from the wall pressure patch value and `rho_w` is evaluated through the thermophysical model EoS using `rhoEoS(p_wall, Tw, faceCells)`.

## Finite-Volume Mapping Choices

`nutGaoWallFunction`:

- Solves Gao Eq. 5/Eq. 6/Eq. 16 for `uTau` and `tauw`.
- Sets OpenFOAM kinematic `nut` so the near-wall discrete shear stress reproduces `tauw` through OpenFOAM's actual wall-stress path: `nut = (tauw - mu_w*|snGrad(U)|)/(rho_patch*|snGrad(U)|)`, clamped non-negative.
- `rho_patch` is `rho.boundaryField()[patchi]`, not Gao's thermodynamic wall density `rho_w`. Gao's `rho_w` is still used inside the wall-law solve, but the final kinematic `nut` must be inverted with the same density that OpenFOAM uses when it computes viscous wall stress. This distinction matters for cold-wall hypersonic cases where `rho_patch` and `rho_w` can differ strongly.
- Gao Eq. 8 is still evaluated as the first-cell dynamic eddy viscosity `mu_t` and is used for Nichols Eq. 22 in `kGaoWallFunction`.

`compressible::alphatGaoWallFunction`:

- Computes `qw` from Gao Eq. 16 after the coupled solve.
- Computes `uHalf` by inverting Gao Eq. 5 at `0.5*yPlus`, then applies Gao Eq. 18.
- By default, sets `alphat` so `alphaEff*he.snGrad()` matches Gao's physical wall heat flux `qw`, after converting from Gao's wall-normal sign convention to OpenFOAM's outward-patch-normal sign convention. The OpenFOAM target flux is `-qw`.
- For comparison/testing, `thermalFluxTarget qHalf;` can instead target Gao Eq. 18's half-cell flux `qHalf`, with OpenFOAM target flux `-qHalf`.
- For OpenFOAM's `heThermo` sensible-internal-energy path, `alphaEff = CpByCpv*(alpha + alphat)`. The target `alphaEff = -qTarget/he.snGrad()` is therefore inverted by `CpByCpv` before assigning `alphat`.
- If the local enthalpy gradient is zero, no finite scalar diffusivity can impose a nonzero flux, so it leaves `alphaEff=alpha`.
- Negative turbulent diffusivity is clamped to zero.

`omegaGaoWallFunction`:

- Derives from OpenFOAM's `omegaWallFunction` to preserve its master-patch and corner-averaging machinery.
- Replaces the omega value with Nichols Eqs. 19-21, using Gao's `uTau`, `kappa`, and `Cmu`.
- The turbulence production `G` insertion is not specified by Gao or Nichols for OpenFOAM. It keeps OpenFOAM's standard wall-function production form, but uses the Gao/Nichols wall `k` and `nut` values.

`kGaoWallFunction`:

- Sets the patch and adjacent cell value from Nichols Eq. 22 using the Gao Eq. 8 first-cell `mu_t`.
- It constrains the scalar equation matrix similarly to a fixed-value wall function.
- Corner averaging across different wall patches is not implemented for `k`; `omega` retains OpenFOAM's existing averaging. This is a remaining OpenFOAM translation choice, not a paper equation.

## Validation-Driven Corrections

The `WF_yPlus_47_91` validation case exposed two finite-volume mapping errors that were not paper-equation errors:

- Thermal mapping: diagnostics showed `qOpenFOAM/qGaoHalf` was exactly about `1.4` on faces where `alphat` was active. This came from assuming `alphaEff = alpha + alphat`, while OpenFOAM's internal-energy thermo computes `alphaEff = CpByCpv*(alpha + alphat)`. The `alphat` wall function now inverts the `CpByCpv` factor before assigning `alphat`.
- Thermal sign convention: Gao's wall-normal coordinate is positive from the wall into the fluid, while OpenFOAM's `he.snGrad()` uses the outward patch normal. For the bottom wall in the flat-plate validation case, the same physical wall heating has positive Gao `qw/qHalf` but negative OpenFOAM `wallHeatFlux`. The `alphat` wall function therefore targets `qOpenFOAM = -qTarget`.
- Thermal target selection: subsequent testing showed that imposing `qHalf` exactly did not improve the coarse-cell temperature profile. The default target was therefore changed to Gao's physical wall heat flux `qw`, while `qHalf` remains available through `thermalFluxTarget qHalf;` for diagnostic comparison.
- Momentum mapping: near the sampled station around `x = 2.5`, Gao's solved `uTau` implied a wall shear close to the pure-SST value, but OpenFOAM's reported `wallShearStress` was much lower. The cause was converting Gao's dynamic target `tauw` to kinematic `nut` with Gao's wall EoS density `rho_w`. OpenFOAM's stress calculation uses the patch `rho` field instead, so `nutShear` is now inverted using `rho.boundaryField()[patchi]`.

These fixes make the OpenFOAM finite-volume operators reproduce the wall-function targets where the targets can be represented by non-negative `nut` and `alphat`. They do not remove the remaining thermal limitation: if the selected OpenFOAM-signed target has the opposite sign from `he.snGrad()`, or requires an effective diffusivity below the laminar value, a non-negative `alphat` wall condition cannot impose it and the value is clamped to zero.

## Numerical Guards

The papers do not specify numerical guards. The implementation adds the following:

- `asin` arguments are clamped to `[-1, 1]` before calling `asin`; the principal branch is used.
- `Gamma`, `Q`, wall temperature, viscosity, conductivity, heat capacity, density, and wall distance are bounded away from zero with OpenFOAM `VSMALL`.
- The White/Gao exponential is clamped to the range `[-50, 50]` before `exp`.
- The Eq. 9 derivative denominator uses `max(1e-12, 1 - arg^2)` after inverse-sine argument clamping. This prevents an out-of-domain or branch-limit value from producing an effectively infinite first-cell eddy viscosity.
- `alphat`, `nut`, and Eq. 8 `mu_t` are clamped non-negative.

These guards prevent floating-point failures near zero velocity, very small heat-transfer parameters, or roundoff at the inverse-sine domain boundary. They are not additional physics.

## Dictionary Entries

All BCs accept the same optional entries:

```foam
kappa       0.41;      // Gao section 2.2.1
B           5.0;       // Gao section 2.2.1
Prt         0.9;       // Gao Eq. 18 discussion
Cmu         0.09;      // Nichols k-omega values
maxIter     50;
tolerance   1e-7;
thermalMode auto;      // auto | isothermal | adiabatic
```

For faithful Gao use, leave `kappa`, `B`, `Prt`, and `Cmu` at their defaults.

`compressible::alphatGaoWallFunction` also accepts diagnostic entries:

```foam
thermalFluxTarget        qw;     // qw | qHalf
diagnostics              false;  // enable Gao/OpenFOAM heat-flux diagnostics
diagnosticLog            true;   // write min/max/mean summary to solver log
diagnosticWrite          true;   // write face-wise .dat output
diagnosticWriteInterval  1;      // time-index interval for diagnostic output
diagnosticWritePrecision 12;
```

When enabled, the face-wise files are written to:

```text
postProcessing/GaoWallFunctionDiagnostics/<patch>/<time>.dat
```

For parallel/decomposed runs, each rank writes to:

```text
postProcessing/GaoWallFunctionDiagnostics/processor<N>/<patch>/<time>.dat
```

The columns are:

```text
face x y z qGaoHalf qOpenFOAM qError qRatio alphat alphaEff heSnGrad qw uTau beta yPlus Gamma
```

where `qOpenFOAM = alphaEff*he.snGrad()`. `qGaoHalf` and `qw` are written with Gao's paper sign convention. `qError` and `qRatio` compare against the selected OpenFOAM-signed target: `-qw` for `thermalFluxTarget qw;`, or `-qGaoHalf` for `thermalFluxTarget qHalf;`.

Example field snippets:

```foam
// constant/turbulenceProperties or system/controlDict must load:
libs ("libGaoWallFunctions.so");

// 0/nut
wall
{
    type        nutGaoWallFunction;
    value       uniform 0;
}

// 0/alphat
wall
{
    type        compressible::alphatGaoWallFunction;
    value       uniform 0;
}

// 0/omega
wall
{
    type        omegaGaoWallFunction;
    value       uniform 1;
}

// 0/k
wall
{
    type        kGaoWallFunction;
    value       uniform 1e-12;
}
```

## Build Status

`wmake libso` completed successfully and wrote:

`/home/mgnadlinger/OpenFOAM/mgnadlinger-v2412/platforms/linux64GccDPInt32Opt/lib/libGaoWallFunctions.so`

No validation case was created or run.
