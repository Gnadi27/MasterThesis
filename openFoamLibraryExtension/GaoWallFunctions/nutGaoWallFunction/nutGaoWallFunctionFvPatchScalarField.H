/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef nutGaoWallFunctionFvPatchScalarField_H
#define nutGaoWallFunctionFvPatchScalarField_H

#include "GaoWallFunction.H"
#include "nutWallFunctionFvPatchScalarField.H"

namespace Foam
{

class nutGaoWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    GaoWallFunction::Coeffs coeffs_;

    virtual tmp<scalarField> calcNut() const;


public:

    TypeName("nutGaoWallFunction");

    nutGaoWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nutGaoWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nutGaoWallFunctionFvPatchScalarField
    (
        const nutGaoWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nutGaoWallFunctionFvPatchScalarField
    (
        const nutGaoWallFunctionFvPatchScalarField&
    );

    nutGaoWallFunctionFvPatchScalarField
    (
        const nutGaoWallFunctionFvPatchScalarField&,
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

    virtual tmp<scalarField> yPlus() const;

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
