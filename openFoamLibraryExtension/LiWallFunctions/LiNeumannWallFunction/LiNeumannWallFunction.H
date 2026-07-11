/*---------------------------------------------------------------------------*\
  Li et al. frozen-composition wall-function support for OpenFOAM v2412.

  Neumann-coupling variant: identical Li et al. closure evaluation, but the
  wall heat flux is imposed through a prescribed-gradient temperature
  boundary condition (LiNeumannWallHeatFluxTemperature) instead of an
  alphat inversion. Because the patch temperature value no longer equals
  the physical wall temperature under that boundary condition, the wall
  state (Tw, rho_w, mu_w, lambda_w, Cp_w, h_w) is evaluated from the
  dictionary Tw and a wallTransport sub-dictionary instead of patch values.
\*---------------------------------------------------------------------------*/

#ifndef LiNeumannWallFunction_H
#define LiNeumannWallFunction_H

#include "fvPatch.H"
#include "labelField.H"
#include "scalarField.H"
#include "turbulentFluidThermoModel.H"
#include "word.H"

namespace Foam
{
class dictionary;
class Ostream;

namespace LiNeumannWallFunction
{

struct Coeffs
{
    scalar kappa;
    scalar B;
    scalar Prt;
    scalar omegaH;
    scalar cv1;
    scalar enthalpyOffset;
    label maxIter;
    scalar tolerance;
    word thermalMode;
    bool strict;

    // Physical wall temperature; the patch T value drifts under the
    // Neumann coupling and can no longer carry it. Negative means unset
    // (fall back to patch values; only valid without the Neumann T BC).
    scalar Tw;

    // Wall transport evaluation at (p_wall, Tw):
    //   fromPatch  - patch mu/kappa (only valid when patch T is pinned
    //                at Tw, i.e. without the Neumann T BC)
    //   sutherland - mu = As*sqrt(Tw)/(1 + Ts/Tw); kappa from Pr if
    //                given, otherwise modified Eucken
    //   constant   - fixed mu and kappa
    word wallTransportMode;
    scalar wallTransportAs;
    scalar wallTransportTs;
    scalar wallTransportPr;
    scalar wallTransportMu;
    scalar wallTransportKappa;

    Coeffs();
    explicit Coeffs(const dictionary& dict);

    void writeEntries(Ostream& os) const;
};


struct PatchData
{
    scalarField uTau;
    scalarField tauw;
    scalarField yPlus;
    scalarField uPlus;
    scalarField recovery;
    scalarField rg;
    scalarField beta;
    scalarField Gamma;
    scalarField Pr;
    scalarField Tw;
    scalarField hw;
    scalarField h1;
    scalarField hRef;
    scalarField hrg;
    scalarField qw;
    scalarField qHalf;
    scalarField uHalf;
    scalarField nutShear;
    scalarField magGradU;
    scalarField tauDeficit;
    scalarField mut1;
    scalarField nuTilda1;
    scalarField residual;
    labelField iterations;
    labelField converged;
    labelField valid;
    labelField referenceCell;
    labelField status;

    explicit PatchData(const label size);
};


// Equation (42), evaluated with fixed beta* and Gamma*.
scalar lawYPlus
(
    scalar uPlus,
    scalar beta,
    scalar Gamma,
    const Coeffs& coeffs
);


// Equation (38), derivative of the generalized White term.
scalar dYWhiteDuPlus
(
    scalar uPlus,
    scalar beta,
    scalar Gamma,
    const Coeffs& coeffs
);


// Equation (39): invert nut = nuTilda*fv1.
scalar solveNuTilda
(
    scalar nutTarget,
    scalar nu,
    scalar cv1,
    label maxIter,
    scalar tolerance
);


// Resolve the thermal mode against the actual T boundary condition,
// recognising the LiNeumannWallHeatFluxTemperature boundary condition.
bool isAdiabatic
(
    const Coeffs& coeffs,
    const fvPatchScalarField& Tpatch
);


// Calculate all Li wall-function quantities for one wall patch.
PatchData calcPatch
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const Coeffs& coeffs
);

} // End namespace LiNeumannWallFunction
} // End namespace Foam

#endif

// ************************************************************************* //
