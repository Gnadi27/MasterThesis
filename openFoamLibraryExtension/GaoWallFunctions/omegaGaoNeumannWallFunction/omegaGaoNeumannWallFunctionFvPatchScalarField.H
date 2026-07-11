/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef omegaGaoNeumannWallFunctionFvPatchScalarField_H
#define omegaGaoNeumannWallFunctionFvPatchScalarField_H

#include "GaoNeumannWallFunction.H"
#include "omegaWallFunctionFvPatchScalarField.H"

namespace Foam
{

class omegaGaoNeumannWallFunctionFvPatchScalarField
:
    public omegaWallFunctionFvPatchScalarField
{
protected:

    GaoNeumannWallFunction::Coeffs coeffs_;

    virtual void calculate
    (
        const turbulenceModel& turbulence,
        const List<scalar>& cornerWeights,
        const fvPatch& patch,
        scalarField& G,
        scalarField& omega
    );


public:

    TypeName("omegaGaoNeumannWallFunction");

    omegaGaoNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    omegaGaoNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    omegaGaoNeumannWallFunctionFvPatchScalarField
    (
        const omegaGaoNeumannWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    omegaGaoNeumannWallFunctionFvPatchScalarField
    (
        const omegaGaoNeumannWallFunctionFvPatchScalarField&
    );

    omegaGaoNeumannWallFunctionFvPatchScalarField
    (
        const omegaGaoNeumannWallFunctionFvPatchScalarField&,
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

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
