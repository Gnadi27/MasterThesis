#ifndef nutLiWallFunctionFvPatchScalarField_H
#define nutLiWallFunctionFvPatchScalarField_H

#include "LiWallFunction.H"
#include "nutWallFunctionFvPatchScalarField.H"

namespace Foam
{

class nutLiWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:
    LiWallFunction::Coeffs coeffs_;

    virtual tmp<scalarField> calcNut() const;

public:
    TypeName("nutLiWallFunction");

    nutLiWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nutLiWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nutLiWallFunctionFvPatchScalarField
    (
        const nutLiWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nutLiWallFunctionFvPatchScalarField
    (
        const nutLiWallFunctionFvPatchScalarField&
    );

    nutLiWallFunctionFvPatchScalarField
    (
        const nutLiWallFunctionFvPatchScalarField&,
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
