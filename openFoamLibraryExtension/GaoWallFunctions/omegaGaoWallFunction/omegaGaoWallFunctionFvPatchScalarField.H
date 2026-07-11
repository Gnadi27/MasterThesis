/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef omegaGaoWallFunctionFvPatchScalarField_H
#define omegaGaoWallFunctionFvPatchScalarField_H

#include "GaoWallFunction.H"
#include "omegaWallFunctionFvPatchScalarField.H"

namespace Foam
{

class omegaGaoWallFunctionFvPatchScalarField
:
    public omegaWallFunctionFvPatchScalarField
{
protected:

    GaoWallFunction::Coeffs coeffs_;

    virtual void calculate
    (
        const turbulenceModel& turbulence,
        const List<scalar>& cornerWeights,
        const fvPatch& patch,
        scalarField& G,
        scalarField& omega
    );


public:

    TypeName("omegaGaoWallFunction");

    omegaGaoWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    omegaGaoWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    omegaGaoWallFunctionFvPatchScalarField
    (
        const omegaGaoWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    omegaGaoWallFunctionFvPatchScalarField
    (
        const omegaGaoWallFunctionFvPatchScalarField&
    );

    omegaGaoWallFunctionFvPatchScalarField
    (
        const omegaGaoWallFunctionFvPatchScalarField&,
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
