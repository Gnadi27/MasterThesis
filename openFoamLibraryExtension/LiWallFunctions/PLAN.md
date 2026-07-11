# Li Frozen-Composition Wall-Function Library

## Summary

Implement a standalone OpenFOAM v2412 user library in `LiWallFunctions` without modifying `/opt/of-v2412`. The library will target single-temperature, ideal-gas, frozen-composition thermodynamics with constant or temperature-dependent properties and Spalart–Allmaras turbulence.

Runtime types:

- `nutLiWallFunction`
- `compressible::alphatLiWallFunction`
- `nuTildaLiWallFunction`

Build output: `$FOAM_USER_LIBBIN/libLiWallFunctions.so`.

## Implementation

- Create a shared Li wall-law evaluator using the paper’s reference-point formulation:
  - Automatically select the second off-wall cell by maximum inward wall-normal displacement, including processor-neighbor access.
  - Use Eq. (40) to infer the generalized recovery enthalpy for isothermal walls.
  - For adiabatic walls, impose \(q_w=0\), obtain \(r_g^*\) from Eq. (41), and use \(h_{rg}=h_w\).
  - Solve Eq. (42) for \(u_\tau\) with bracketed iteration using \(\kappa=0.41\), \(B=6.0\), and \(\omega^*=0.72\).
  - Compute \(\tau_w\), \(q_w\), \(q_{1/2}\), and the half-cell velocity from Eqs. (33)–(35).
  - Evaluate first-cell eddy viscosity from Eqs. (37)–(38), then solve Eq. (39) for SA `nuTilda` using \(c_{v1}=7.1\).

- Map the results into OpenFOAM:
  - `nutLiWallFunction` converts \(\tau_w\) to the wall `nut` needed by OpenFOAM’s discrete stress operator.
  - `nuTildaLiWallFunction` keeps the physical wall value zero and constrains the adjacent cell to Eq. (39).
  - `alphatLiWallFunction` converts the Li energy-flux target to OpenFOAM `alphaEff`, including `CpByCpv` for internal-energy thermodynamics and the opposite wall-normal sign convention.
  - Isothermal mode targets paper Eq. (34), \(q_{1/2}\), by default; `qw` remains an optional diagnostic/compatibility target.
  - Adiabatic mode imposes physical \(q_w=0\); the generally nonzero \(q_{1/2}\) is diagnostic only because applying it at OpenFOAM’s physical boundary face would violate adiabaticity.
  - Clamp negative equivalent `nut`/`alphat` to zero and report when the requested flux cannot be represented by non-negative diffusivity.

- Reconstruct frozen-gas enthalpy without OpenFOAM’s arbitrary reference datum:
  - Set \(h_w=C_{p,w}T_w+h_{\mathrm{offset}}\).
  - Obtain first/reference enthalpies from exact OpenFOAM sensible-enthalpy differences relative to the wall.
  - This gives \(h=C_pT\) exactly for a calorically perfect gas and retains temperature-dependent enthalpy differences for calorically imperfect gases.
  - Default `enthalpyOffset` is zero and remains configurable.

- Expose shared options on all boundary conditions:
  - `kappa 0.41`, `B 6.0`, `Prt 0.9`, `omegaH 0.72`, `cv1 7.1`
  - `thermalMode auto|isothermal|adiabatic`
  - nonlinear iteration tolerance and limit
  - optional strict failure handling and diagnostics
  - `thermalFluxTarget qHalf|qw` for isothermal walls

- Add diagnostics containing reference-cell identity, \(u_\tau\), \(y^+\), \(r_g^*\), \(\beta^*\), \(\Gamma^*\), \(\tau_w\), \(q_w\), \(q_{1/2}\), \(\mu_{t1}\), `nuTilda1`, residuals, and any fallback/clamping status.

## Simplifications and Documented Deviations

- Remove species transport, diffusion enthalpy fluxes, chemical source terms, reaction kinetics, and catalytic-wall behavior from Eq. (1).
- Treat composition as fixed, making \(dY_s/dy=0\) in Eq. (24) exact. OpenFOAM supplies mixture \(C_p(T)\), \(\mu(T)\), conductivity, density, and sensible enthalpy at the configured frozen composition.
- Retain Eqs. (25), (33)–(42); chemical-reaction effects disappear from their interpretation, while cold-wall and variable-property effects remain.
- Do not implement Gupta transport fits or Wilke mixing internally; use the active OpenFOAM thermophysical and transport model.
- Retain the paper’s \(\omega^*=0.72\), but flag that it was fitted for reacting high-temperature air and is not established as universal for arbitrary frozen gases.
- Retain Eq. (29). Frozen composition makes \(R=R_w\), but \(C_p/C_{p,w}\neq1\) for calorically imperfect gases, so the density–enthalpy relation remains approximate.
- Document the enthalpy datum reconstruction because the paper does not define a datum sufficiently for direct use with OpenFOAM sensible-energy models.
- Document second-cell selection as a cell-centered OpenFOAM interpretation of the paper’s “second grid point.”
- Document the \(q_{1/2}\)-versus-\(q_w\) boundary-face ambiguity and the adiabatic choice above.
- Limit v1 to smooth-wall, attached, approximately zero-pressure-gradient equilibrium boundary layers and SA turbulence. Separation, stagnation regions, rough walls, chemistry, SST, and species wall functions are out of scope.
- On zero tangential velocity, invalid/non-monotonic reference states, or an unavailable second cell, use a reported laminar fallback by default; strict mode terminates with a detailed error.

## Tests and Acceptance

- Build successfully with `source /opt/of-v2412/etc/bashrc && wmake libso`.
- Add equation-level regression tests for Eqs. (33)–(42), finite-difference verification of Eq. (38), nonlinear-root recovery, and Eq. (39).
- Verify the calorically perfect limit produces \(h=C_pT\), and a JANAF/polynomial-transport case preserves \(dh/dT=C_p(T)\).
- Add small compressible SA integration cases for:
  - isothermal wall and \(q_{1/2}\) coupling;
  - adiabatic wall with zero physical heat flux;
  - first cells in viscous and logarithmic-layer ranges;
  - serial and two-rank decomposed execution.
- Check that reported OpenFOAM wall stress and energy flux match the Li targets wherever non-negative equivalent diffusivities exist.
- Provide `README.md`, example field dictionaries, and `DESIGN_AND_DEVIATIONS.md` with an equation-to-code map and all limitations.
- Do not claim reproduction of the paper’s reacting Wedge15H25/Wedge20H25 results, since chemistry is deliberately excluded.
