/*---------------------------------------------------------------------------*\
  Momentum coupling for the Griffin wall function: nut is set so that
  OpenFOAM's discrete wall stress (mu + rho*nut)*snGrad(U) reproduces the
  model tau_w exactly.
\*---------------------------------------------------------------------------*/

#ifndef nutGriffinWallFunctionFvPatchScalarField_H
#define nutGriffinWallFunctionFvPatchScalarField_H

#include "GriffinWallFunction.H"
#include "nutWallFunctionFvPatchScalarField.H"

namespace Foam
{

class nutGriffinWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:
    GriffinWallFunction::Coeffs coeffs_;

    virtual tmp<scalarField> calcNut() const;

public:
    TypeName("nutGriffinWallFunction");

    nutGriffinWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nutGriffinWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nutGriffinWallFunctionFvPatchScalarField
    (
        const nutGriffinWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nutGriffinWallFunctionFvPatchScalarField
    (
        const nutGriffinWallFunctionFvPatchScalarField&
    );

    nutGriffinWallFunctionFvPatchScalarField
    (
        const nutGriffinWallFunctionFvPatchScalarField&,
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
