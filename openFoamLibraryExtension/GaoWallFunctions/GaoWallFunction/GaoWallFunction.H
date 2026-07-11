/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef GaoWallFunction_H
#define GaoWallFunction_H

#include "fvPatch.H"
#include "labelField.H"
#include "scalarField.H"
#include "turbulentFluidThermoModel.H"
#include "word.H"

namespace Foam
{
class dictionary;
class Ostream;

namespace GaoWallFunction
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


//- Calculate all Gao/Nichols wall-function quantities for one wall patch
PatchData calcPatch
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const Coeffs& coeffs
);

} // End namespace GaoWallFunction
} // End namespace Foam

#endif

// ************************************************************************* //
