#ifndef nuTildaLiNeumannWallFunctionFvPatchScalarField_H
#define nuTildaLiNeumannWallFunctionFvPatchScalarField_H

#include "LiNeumannWallFunction.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

class nuTildaLiNeumannWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    LiNeumannWallFunction::Coeffs coeffs_;

    void checkType() const;

public:
    TypeName("nuTildaLiNeumannWallFunction");

    nuTildaLiNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nuTildaLiNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nuTildaLiNeumannWallFunctionFvPatchScalarField
    (
        const nuTildaLiNeumannWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nuTildaLiNeumannWallFunctionFvPatchScalarField
    (
        const nuTildaLiNeumannWallFunctionFvPatchScalarField&
    );

    nuTildaLiNeumannWallFunctionFvPatchScalarField
    (
        const nuTildaLiNeumannWallFunctionFvPatchScalarField&,
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
