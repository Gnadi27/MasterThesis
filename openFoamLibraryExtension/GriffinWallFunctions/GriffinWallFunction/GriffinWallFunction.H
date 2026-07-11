/*---------------------------------------------------------------------------*\
  Griffin, Fu & Moin wall-model support for OpenFOAM v2412.

  RANS (Spalart-Allmaras) deployment of the near-wall model of
  Griffin, Fu & Moin, "Near-wall model for compressible turbulent boundary
  layers based on an inverse velocity transformation", CTR Annual Research
  Briefs 2022 (paper Eqs. 2.1-2.8), built on
   - the total-stress-based velocity transformation of Griffin, Fu & Moin,
     PNAS 118 (2021), Eqs. 4-5 (constant-stress-layer form), inverted as
     paper Eq. 2.5;
   - the algebraic temperature-velocity relation of Zhang et al. (2014)
     with matching-point correction, paper Eq. 2.8;
   - Reynolds-number dependent incompressible coefficients kappa(Re_tau*),
     A+(Re_tau*) calibrated on Lee & Moser (2015) channel DNS;
   - boundary-layer edge quantities from the local-reconstruction method of
     Griffin, Fu & Moin, PRF 6 (2021), compressible (stagnation-enthalpy)
     form of its Appendix D (see GriffinBoundaryLayerEdge).

  The model returns the wall shear stress tau_w by shooting on the inverse
  transformation ODE and the wall heat flux algebraically,
      q_w = s tau_w Cp (T_w - T_r)/U_e,
  which is positive for heat entering the fluid (Griffin's B_q sign
  convention; this coincides with OpenFOAM's outward-normal boundary-flux
  sign, so no sign flip is needed when imposing it).

  Energy coupling is Neumann-only: the heat flux is imposed through the
  prescribed-gradient temperature boundary condition
  GriffinWallHeatFluxTemperature, never through an alphat inversion (see
  LiWallFunctions/FINDINGS_M6Tw076_energy_coupling.md for why). The patch
  T value therefore drifts and the physical wall temperature must be
  supplied through the coefficients (Tw), and viscosity through an
  explicit temperature law (the inner integration needs mu(T) at arbitrary
  T in any case).
\*---------------------------------------------------------------------------*/

#ifndef GriffinWallFunction_H
#define GriffinWallFunction_H

#include "fvPatch.H"
#include "labelField.H"
#include "scalarField.H"
#include "turbulentFluidThermoModel.H"
#include "word.H"

namespace Foam
{
class dictionary;
class Ostream;

namespace GriffinWallFunction
{

struct Coeffs
{
    // Generalized Reynolds-analogy factor s of paper Eq. 2.8 (Zhang et al.
    // 2014 convention).
    scalar s;

    // Molecular Prandtl number: enters Eq. 2.8 as s*Pr and the recovery
    // factor r = Pr^(1/3).
    scalar Pr;

    // Spalart-Allmaras fv1 constant for the nuTilda first-cell inversion.
    scalar cv1;

    // Shooting/root-finding controls.
    label maxIter;
    scalar tolerance;

    // Inner-profile integration: number of log-spaced intervals between
    // the viscous seed point and the matching location, and the number of
    // per-node property sub-iterations (paper Sec. 2.3).
    label integrationNodes;
    label subIterations;

    word thermalMode;   // auto | isothermal | adiabatic
    bool strict;

    // Physical wall temperature; the patch T value drifts under the
    // Neumann coupling and can no longer carry it. Negative means unset.
    scalar Tw;

    // Viscosity law mu(T) used at the wall and along the inner profile:
    //   sutherland - mu = As*sqrt(T)/(1 + Ts/T)
    //   powerLaw   - mu = muRef*(T/TRef)^n
    //   constant   - fixed mu
    word viscosityMode;
    scalar viscosityAs;
    scalar viscosityTs;
    scalar viscosityMuRef;
    scalar viscosityTRef;
    scalar viscosityN;
    scalar viscosityMu;

    // Boundary-layer edge:
    //   reconstruct - PRF-2021 local reconstruction along wall-normal cell
    //                 walks (compressible/stagnation-enthalpy form)
    //   prescribed  - uniform Ue, Te, delta from this dictionary
    word edgeMode;
    scalar Ue;
    scalar Te;
    scalar delta;
    label maxEdgeCells;
    scalar maxEdgeDistance;
    label edgeUpdateInterval;
    scalar edgeThreshold;   // U/U_I at the edge; paper n/100 = 0.99

    Coeffs();
    explicit Coeffs(const dictionary& dict);

    void writeEntries(Ostream& os) const;

    // Viscosity law evaluation.
    scalar mu(const scalar T) const;
    scalar dmudT(const scalar T) const;
};


struct PatchData
{
    scalarField uTau;
    scalarField tauw;
    scalarField yPlus;
    scalarField yStar;
    scalarField uPlus;
    scalarField ReTauStar;
    scalarField kappaCal;
    scalarField APlusCal;
    scalarField Tw;
    scalarField Tr;
    scalarField T1;
    scalarField U1;
    scalarField Ue;
    scalarField Te;
    scalarField delta;
    scalarField qw;         // Griffin/OpenFOAM sign: > 0 into the fluid
    scalarField nutShear;
    scalarField mut1;
    scalarField nuTilda1;
    scalarField residual;
    labelField iterations;
    labelField converged;
    labelField valid;
    labelField edgeStatus;
    labelField status;

    explicit PatchData(const label size);
};


// Interpolate the calibrated incompressible coefficients (paper Sec. 2):
// linear in Re_tau*, constant extrapolation outside [182, 5190].
void calibration(const scalar ReTauStar, scalar& kappa, scalar& APlus);


// Van Driest-damped mixing length of paper Eq. 2.2 in semi-local units:
// mut+ = kappa y* (1 - exp(-y*/A+))^2. The incompressible constant-stress
// shear follows algebraically as S_t+ = 1/(1 + mut+).
scalar mutPlusIncompressible
(
    const scalar yStar,
    const scalar kappa,
    const scalar APlus
);


// Algebraic temperature-velocity relation, paper Eq. 2.8.
scalar temperatureFromU
(
    const scalar U,
    const scalar UM,
    const scalar TM,
    const scalar Ue,
    const scalar Tw,
    const scalar Tr,
    const scalar s,
    const scalar Pr
);


// dT/dU of paper Eq. 2.8 (for the analytic property-gradient terms of the
// inverse transformation, paper Eq. 2.5).
scalar dTemperaturedU
(
    const scalar U,
    const scalar UM,
    const scalar TM,
    const scalar Ue,
    const scalar Tw,
    const scalar Tr,
    const scalar s,
    const scalar Pr
);


// Recovery temperature T_r = T_e + r U_e^2/(2 Cp), r = Pr^(1/3).
scalar recoveryTemperature
(
    const scalar Te,
    const scalar Ue,
    const scalar Cp,
    const scalar Pr
);


// Spalart-Allmaras inversion nut = nuTilda*fv1(nuTilda/nu).
scalar solveNuTilda
(
    const scalar nutTarget,
    const scalar nu,
    const scalar cv1,
    const label maxIter,
    const scalar tolerance
);


// Per-face inner problem: everything the 1D model needs, mesh-free.
struct InnerInput
{
    scalar yM;          // matching distance (first-cell centre)
    scalar UM;          // wall-tangential velocity magnitude at yM
    scalar TM;          // temperature at yM
    scalar rhoM;        // density at yM (fixes p = rhoM R TM, constant)
    scalar Tw;          // physical wall temperature (ignored if adiabatic)
    bool adiabatic;
    scalar Ue;          // boundary-layer edge velocity
    scalar Te;          // boundary-layer edge temperature
    scalar delta;       // boundary-layer thickness
    scalar Cp;          // specific heat (calorically-perfect-gas form)
    scalar tauwGuess;   // initial wall shear stress guess
};


struct InnerResult
{
    scalar uTau;
    scalar tauw;
    scalar qw;          // > 0 into the fluid (Griffin B_q sign convention)
    scalar Tr;
    scalar Tw;          // effective wall temperature (= Tr if adiabatic)
    scalar yPlus;
    scalar yStar;
    scalar ReTauStar;
    scalar kappa;
    scalar APlus;
    scalar mut1;        // dynamic eddy viscosity at yM from the profile
    scalar residual;
    label iterations;
    label converged;
};


// Solve the inner problem: integrate the inverse-transformation ODE
// (paper Eq. 2.5) with per-node algebraic T(U) property sub-iteration and
// shoot on tau_w such that U(yM) = UM (paper Sec. 2.3).
InnerResult solveInner(const InnerInput& in, const Coeffs& coeffs);


// Resolve the thermal mode against the actual T boundary condition,
// recognising the GriffinWallHeatFluxTemperature boundary condition.
bool isAdiabatic
(
    const Coeffs& coeffs,
    const fvPatchScalarField& Tpatch
);


// Calculate all Griffin wall-function quantities for one wall patch.
PatchData calcPatch
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const Coeffs& coeffs
);

} // End namespace GriffinWallFunction
} // End namespace Foam

#endif

// ************************************************************************* //
