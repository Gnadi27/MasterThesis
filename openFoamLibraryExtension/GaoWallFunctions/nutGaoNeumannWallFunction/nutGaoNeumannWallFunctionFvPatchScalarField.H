#ifndef nutGaoNeumannWallFunctionFvPatchScalarField_H
#define nutGaoNeumannWallFunctionFvPatchScalarField_H

#include "GaoNeumannWallFunction.H"
#include "nutWallFunctionFvPatchScalarField.H"

namespace Foam
{

class nutGaoNeumannWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    GaoNeumannWallFunction::Coeffs coeffs_;

    virtual tmp<scalarField> calcNut() const;


public:

    TypeName("nutGaoNeumannWallFunction");

    nutGaoNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nutGaoNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nutGaoNeumannWallFunctionFvPatchScalarField
    (
        const nutGaoNeumannWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nutGaoNeumannWallFunctionFvPatchScalarField
    (
        const nutGaoNeumannWallFunctionFvPatchScalarField&
    );

    nutGaoNeumannWallFunctionFvPatchScalarField
    (
        const nutGaoNeumannWallFunctionFvPatchScalarField&,
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
