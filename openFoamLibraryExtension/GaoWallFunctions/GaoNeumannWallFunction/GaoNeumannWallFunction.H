/*---------------------------------------------------------------------------*\
  Gao et al. compressible wall-function support for OpenFOAM v2412.

  Neumann-coupling variant: identical Gao/Nichols closure evaluation, but
  the wall heat flux is imposed through a prescribed-gradient temperature
  boundary condition (GaoNeumannWallHeatFluxTemperature) instead of an
  alphat inversion. Because the patch temperature value no longer equals
  the physical wall temperature under that boundary condition, the wall
  state (Tw, rho_w, mu_w, k_w, cp_w) is evaluated from the dictionary Tw
  and a wallTransport sub-dictionary instead of patch values.
\*---------------------------------------------------------------------------*/

#ifndef GaoNeumannWallFunction_H
#define GaoNeumannWallFunction_H

#include "fvPatch.H"
#include "labelField.H"
#include "scalarField.H"
#include "turbulentFluidThermoModel.H"
#include "word.H"

namespace Foam
{
class dictionary;
class Ostream;

namespace GaoNeumannWallFunction
{

struct Coeffs
{
    scalar kappa;
    scalar B;
    scalar Prt;
    scalar Cmu;
    label maxIter;
    scalar tolerance;
    word thermalMode;

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
    scalarField Gamma;
    scalarField beta;
    scalarField Pr;
    scalarField recovery;
    scalarField Tw;
    scalarField cp;
    scalarField rhoW;
    scalarField qw;
    scalarField qHalf;
    scalarField uHalf;
    scalarField nutShear;
    scalarField mutEq8;
    scalarField omega;
    scalarField k;
    scalarField residual;
    labelField iterations;
    labelField converged;

    explicit PatchData(const label size);
};


//- Resolve the thermal mode against the actual T boundary condition,
//  recognising the GaoNeumannWallHeatFluxTemperature boundary condition.
bool isAdiabatic
(
    const Coeffs& coeffs,
    const fvPatchScalarField& Tpatch
);


//- Calculate all Gao/Nichols wall-function quantities for one wall patch
PatchData calcPatch
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const Coeffs& coeffs
);

} // End namespace GaoNeumannWallFunction
} // End namespace Foam

#endif

// ************************************************************************* //
