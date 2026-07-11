/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef kGaoNeumannWallFunctionFvPatchScalarField_H
#define kGaoNeumannWallFunctionFvPatchScalarField_H

#include "GaoNeumannWallFunction.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

class kGaoNeumannWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    GaoNeumannWallFunction::Coeffs coeffs_;

    void checkType() const;


public:

    TypeName("kGaoNeumannWallFunction");

    kGaoNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    kGaoNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    kGaoNeumannWallFunctionFvPatchScalarField
    (
        const kGaoNeumannWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    kGaoNeumannWallFunctionFvPatchScalarField
    (
        const kGaoNeumannWallFunctionFvPatchScalarField&
    );

    kGaoNeumannWallFunctionFvPatchScalarField
    (
        const kGaoNeumannWallFunctionFvPatchScalarField&,
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

    virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
