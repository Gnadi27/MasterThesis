/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef alphatGaoWallFunctionFvPatchScalarField_H
#define alphatGaoWallFunctionFvPatchScalarField_H

#include "GaoWallFunction.H"
#include "fixedValueFvPatchFields.H"
#include "turbulentFluidThermoModel.H"

namespace Foam
{
namespace compressible
{

class alphatGaoWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    GaoWallFunction::Coeffs coeffs_;

    word thermalFluxTarget_;

    bool diagnostics_;

    bool diagnosticLog_;

    bool diagnosticWrite_;

    label diagnosticWriteInterval_;

    label diagnosticWritePrecision_;

    label diagnosticTimeIndex_;

    void checkType() const;

    void checkThermalFluxTarget() const;

    void writeDiagnostics
    (
        const scalarField& alphaEff,
        const scalarField& heSnGrad,
        const GaoWallFunction::PatchData& data
    );


public:

    TypeName("compressible::alphatGaoWallFunction");

    alphatGaoWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    alphatGaoWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    alphatGaoWallFunctionFvPatchScalarField
    (
        const alphatGaoWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    alphatGaoWallFunctionFvPatchScalarField
    (
        const alphatGaoWallFunctionFvPatchScalarField&
    );

    alphatGaoWallFunctionFvPatchScalarField
    (
        const alphatGaoWallFunctionFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchField<scalar>> clone() const
    {
        return fvPatchField<scalar>::Clone(*this);
    }

    virtual tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return fvPatchField<scalar>::Clone(*this, iF);
    }

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;
};

} // End namespace compressible
} // End namespace Foam

#endif

// ************************************************************************* //
