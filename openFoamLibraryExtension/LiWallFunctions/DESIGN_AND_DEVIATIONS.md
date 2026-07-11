# Li Wall Functions - Design, Simplifications, and Deviations

## Scope

The library implements the Li et al. wall-function closure for:

- one translational/rotational temperature;
- frozen composition;
- an ideal-gas equation of state;
- constant or temperature-dependent heat capacity and transport properties;
- smooth-wall Spalart-Allmaras RANS calculations.

The implementation is a user library and does not modify OpenFOAM under
`/opt/of-v2412`.

## Frozen-composition reduction

The species equations, reaction source terms, species-diffusion enthalpy flux,
reaction mechanism, and catalytic-wall model in paper Eq. (1) are removed.
Composition is an input of the selected OpenFOAM thermophysical model and must
remain fixed. Consequently, `dYs/dy=0` in Eq. (24) is exact and

```text
dh/dy = Cp(T,Yfixed) dT/dy.
```

OpenFOAM supplies density, `Cp`, molecular viscosity, and thermal conductivity.
The library does not contain Gupta transport fits, Wilke mixing, or chemistry.
A frozen multi-species mixture can be used if its OpenFOAM thermo presents the
mixture through `fluidThermo`; no species boundary conditions are added.

The generalized enthalpy-velocity law is retained. In this reduced model its
generalized recovery factor describes cold-wall and variable-property effects,
not dissociation or recombination.

## Enthalpy definition

Li et al. use static enthalpy ratios but do not define a datum that can be
mapped unambiguously to OpenFOAM's sensible-energy formulations. The library
uses one wall-local frozen sensible-enthalpy datum:

```text
hw   = Cp,w Tw + enthalpyOffset
h1   = hw + (hs,1   - hs,w)
href = hw + (hs,ref - hs,w)
```

For an internal-energy thermo, `hs=e+p/rho`; for an enthalpy thermo, `hs=h`.
This exactly recovers `h=Cp*T` for a calorically perfect ideal gas and keeps the
thermo model's exact sensible-enthalpy differences when `Cp(T)` varies. It
removes the arbitrary OpenFOAM sensible-energy reference and the constant
formation enthalpy, which cannot affect a frozen-composition energy equation.

This datum is a documented modelling choice because Li's ratios are not
invariant to an additive enthalpy constant. `enthalpyOffset` permits controlled
comparison with a different convention.

## Equation mapping

For isothermal walls, the first-cell state and automatic second-cell reference
are inserted into Eq. (40). Solving Eq. (40) gives `hrg`; Eq. (41) then gives
`rg`. For adiabatic walls, `qw=0`, `hrg=hw`, and Eq. (41) gives `rg` directly.

The shared evaluator then performs:

1. Solve Eq. (42) for `uTau` by bracketed bisection.
2. Evaluate `tauw=rho_w*uTau^2`, `beta`, `Gamma`, and `yPlus`.
3. Evaluate `Prw=mu_w*Cp,w/lambda_w` from wall properties, not a constant.
4. Evaluate `qw` from Eq. (33), except that adiabatic mode enforces `qw=0`.
5. Invert Eq. (42) at half of the first-cell `yPlus`, then evaluate Eq. (34).
6. Evaluate Eqs. (37)-(38) with the actual first-cell `mu1/muw`.
7. Invert Eq. (39) for the SA first-cell `nuTilda`, using `cv1=7.1`.

The paper constants are the defaults: `kappa=0.41`, `B=6.0`, `Prt=0.9`, and
`omegaH=0.72`.

## OpenFOAM finite-volume mapping

### Momentum

Li's `tauw` replaces the wall-adjacent shear flux. The boundary `nut` is
obtained by inverting OpenFOAM's actual discrete stress relation:

```text
nut = (tauw - muw*|snGrad(Ut)|)/(rhoPatch*|snGrad(Ut)|).
```

`rhoPatch` is used in this last conversion because it is the density used by
OpenFOAM's stress operator. The thermodynamic wall density remains the density
used in `uTau` and `yPlus`.

### Energy

Li's sign convention points from the wall into the fluid. OpenFOAM's patch
normal points out of the fluid, so the target OpenFOAM flux is `-qw` or
`-qHalf`. For sensible-internal-energy thermodynamics, the inversion includes
OpenFOAM's `CpByCpv` factor in `alphaEff`.

Li's finite-difference diagram places `q1/2` between the wall node and the first
off-wall point. An OpenFOAM boundary face is the physical wall. The default
uses `qHalf` for isothermal cases to follow paper Fig. 16, but uses physical
`qw=0` for adiabatic cases so that adiabaticity is not violated. Both values are
reported.

### Spalart-Allmaras

The `nuTilda` patch value remains zero. The adjacent cell equation is
intended to be constrained to Eq. (39) through `manipulateMatrix`, matching the
paper's first-grid-point prescription without assigning a nonzero value at the
wall. Stock OpenFOAM v2412 `SpalartAllmarasBase::correct()` does not call
`boundaryManipulate()` before solving, so the callback is not executed and the
constraint is inactive unless the SA solve is patched to invoke it after
relaxation/`fvOptions.constrain()` and before `solve()`.
The local `/opt/of-v2412` source and compressible turbulence library contain
this patch as of 2026-07-11; externally installed or cluster OpenFOAM builds
must be checked independently.

## Reference-point selection

For every wall face, the boundary owner is point 1. Candidate cell neighbours
are ranked by the projection of their centre displacement on the inward wall
normal; the largest positive projection is point 2. Processor-patch neighbours
and their field values are exchanged through `patchNeighbourField` and are
eligible candidates.

The candidate must be farther inward and have tangential speed greater than
point 1. This is the cell-centred OpenFOAM interpretation of Sec. IV.C and is
not specified by the paper for skewed or general polyhedral meshes.

## Numerical guards and fallback

- Inverse-sine arguments are clamped to `[-1,1]`.
- Exponential arguments are limited to `[-50,50]`.
- Positive thermodynamic denominators are protected from zero.
- Eq. (42) and Eq. (39) use bracketed scalar solves.
- Negative equivalent `nut`, `alphat`, and `mut1` are clamped to zero.
- Invalid reference states use zero turbulent transport unless `strict=true`.

These guards are implementation choices, not additional paper physics. A
non-negative scalar diffusivity cannot reproduce a target below the laminar
flux or a target with the wrong gradient sign; diagnostics expose the mismatch.

## Remaining paper assumptions and ambiguities

- `omegaH=0.72` was fitted to the reacting-air cases in the paper. Its accuracy
  is not established for arbitrary frozen gases or transport models.
- Eq. (29) is retained. Frozen composition makes `R/Rw=1`, but the neglected
  `Cp,w/Cp` variation remains for calorically imperfect gases.
- The constant-shear, attached, smooth-wall, approximately zero-pressure-
  gradient assumptions are retained.
- The method is not validated here for stagnation, separation, roughness,
  strong pressure gradients, catalytic walls, or reacting composition.
- Only Spalart-Allmaras is supported because Eq. (39) is model-specific.
- The reacting Wedge15H25 and Wedge20H25 results cannot be reproduced by this
  deliberately frozen-composition implementation and are not claimed as
  validation data.
