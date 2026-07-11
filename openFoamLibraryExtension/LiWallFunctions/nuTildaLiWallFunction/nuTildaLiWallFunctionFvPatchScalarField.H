#ifndef nuTildaLiWallFunctionFvPatchScalarField_H
#define nuTildaLiWallFunctionFvPatchScalarField_H

#include "LiWallFunction.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

class nuTildaLiWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    LiWallFunction::Coeffs coeffs_;

    void checkType() const;

public:
    TypeName("nuTildaLiWallFunction");

    nuTildaLiWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nuTildaLiWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nuTildaLiWallFunctionFvPatchScalarField
    (
        const nuTildaLiWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nuTildaLiWallFunctionFvPatchScalarField
    (
        const nuTildaLiWallFunctionFvPatchScalarField&
    );

    nuTildaLiWallFunctionFvPatchScalarField
    (
        const nuTildaLiWallFunctionFvPatchScalarField&,
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
