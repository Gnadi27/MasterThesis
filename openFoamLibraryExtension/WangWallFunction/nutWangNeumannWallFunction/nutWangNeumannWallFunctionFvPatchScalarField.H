/*---------------------------------------------------------------------------*\
  nut boundary condition of the Wang et al. (2023) AWF^Cmp Neumann family.

  Sets the kinematic turbulent viscosity so that OpenFOAM's discrete
  wall-stress path reproduces the analytical wall shear stress
  tau_w = -(rho_v sqrt(kP)/mu_v) A1 (thesis Eq. T8.12, paper Eqs. WF30-36):

      nut = max(0, (tau_w - mu_patch*|snGrad(U)|)
                  /(rho_patch*|snGrad(U)|))

  The inversion deliberately uses the *patch-field* mu and rho (the values
  the stress operator actually uses), while the closure itself is evaluated
  at the physical wall/sublayer state (Tw, Tv, wallTransport).

  All WangNeumann* conditions on a patch must carry identical coefficients.
\*---------------------------------------------------------------------------*/

#ifndef nutWangNeumannWallFunctionFvPatchScalarField_H
#define nutWangNeumannWallFunctionFvPatchScalarField_H

#include "WangNeumannWallFunction.H"
#include "nutWallFunctionFvPatchScalarField.H"

namespace Foam
{

class nutWangNeumannWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    WangNeumannWallFunction::Coeffs coeffs_;

    virtual tmp<scalarField> calcNut() const;


public:

    TypeName("nutWangNeumannWallFunction");

    nutWangNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nutWangNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nutWangNeumannWallFunctionFvPatchScalarField
    (
        const nutWangNeumannWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nutWangNeumannWallFunctionFvPatchScalarField
    (
        const nutWangNeumannWallFunctionFvPatchScalarField&
    );

    nutWangNeumannWallFunctionFvPatchScalarField
    (
        const nutWangNeumannWallFunctionFvPatchScalarField&,
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

    //- Reports the AWF non-dimensional wall distance y* (paper Eq. WF26)
    virtual tmp<scalarField> yPlus() const;

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
