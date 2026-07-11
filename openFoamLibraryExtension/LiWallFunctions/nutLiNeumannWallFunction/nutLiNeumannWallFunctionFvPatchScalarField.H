#ifndef nutLiNeumannWallFunctionFvPatchScalarField_H
#define nutLiNeumannWallFunctionFvPatchScalarField_H

#include "LiNeumannWallFunction.H"
#include "nutWallFunctionFvPatchScalarField.H"

namespace Foam
{

class nutLiNeumannWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:
    LiNeumannWallFunction::Coeffs coeffs_;

    virtual tmp<scalarField> calcNut() const;

public:
    TypeName("nutLiNeumannWallFunction");

    nutLiNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nutLiNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nutLiNeumannWallFunctionFvPatchScalarField
    (
        const nutLiNeumannWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nutLiNeumannWallFunctionFvPatchScalarField
    (
        const nutLiNeumannWallFunctionFvPatchScalarField&
    );

    nutLiNeumannWallFunctionFvPatchScalarField
    (
        const nutLiNeumannWallFunctionFvPatchScalarField&,
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
