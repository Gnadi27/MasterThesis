/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef kGaoWallFunctionFvPatchScalarField_H
#define kGaoWallFunctionFvPatchScalarField_H

#include "GaoWallFunction.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

class kGaoWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    GaoWallFunction::Coeffs coeffs_;

    void checkType() const;


public:

    TypeName("kGaoWallFunction");

    kGaoWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    kGaoWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    kGaoWallFunctionFvPatchScalarField
    (
        const kGaoWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    kGaoWallFunctionFvPatchScalarField
    (
        const kGaoWallFunctionFvPatchScalarField&
    );

    kGaoWallFunctionFvPatchScalarField
    (
        const kGaoWallFunctionFvPatchScalarField&,
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
