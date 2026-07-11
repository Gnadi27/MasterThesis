/*---------------------------------------------------------------------------*\
  Li et al. frozen-composition wall-function support for OpenFOAM v2412.
\*---------------------------------------------------------------------------*/

#ifndef LiWallFunction_H
#define LiWallFunction_H

#include "fvPatch.H"
#include "labelField.H"
#include "scalarField.H"
#include "turbulentFluidThermoModel.H"
#include "word.H"

namespace Foam
{
class dictionary;
class Ostream;

namespace LiWallFunction
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


// Calculate all Li wall-function quantities for one wall patch.
PatchData calcPatch
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const Coeffs& coeffs
);

} // End namespace LiWallFunction
} // End namespace Foam

#endif

// ************************************************************************* //
