/*---------------------------------------------------------------------------*\
  omega boundary condition of the Wang et al. (2023) AWF^Cmp Neumann
  family, for k-omega models (k-omega SST in particular).

  NOT part of the Wang/Craft papers: this is a translation of their
  k-epsilon near-wall treatment (see epsilonWangNeumannWallFunction) to
  omega, flagged as a deviation in DESIGN_AND_DEVIATIONS.md:

      omega_P = epsilon_P/(Cmu kP) = sqrt(kP)/(Cmu c_l y_P)

  With the k-equation sink beta* omega k and beta* = Cmu, the implicit
  near-wall cell sink equals epsilon_P, so the generation field carries
  the same combination as the epsilon variant:

      G = Pk/rho + epsilon_P - epsBar

  making the net k source Pk/rho - epsBar at convergence. Derives from
  omegaWallFunction to reuse its master-patch/corner-averaging and
  matrix-constraint machinery. No Nichols-style inner/outer blending is
  applied. Use with kqRWallFunction on k.
\*---------------------------------------------------------------------------*/

#ifndef omegaWangNeumannWallFunctionFvPatchScalarField_H
#define omegaWangNeumannWallFunctionFvPatchScalarField_H

#include "WangNeumannWallFunction.H"
#include "omegaWallFunctionFvPatchScalarField.H"

namespace Foam
{

class omegaWangNeumannWallFunctionFvPatchScalarField
:
    public omegaWallFunctionFvPatchScalarField
{
protected:

    WangNeumannWallFunction::Coeffs coeffs_;

    //- Calculate the omega and G contributions (AWF values)
    virtual void calculate
    (
        const turbulenceModel& turbulence,
        const List<scalar>& cornerWeights,
        const fvPatch& patch,
        scalarField& G,
        scalarField& omega
    );


public:

    TypeName("omegaWangNeumannWallFunction");

    omegaWangNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    omegaWangNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    omegaWangNeumannWallFunctionFvPatchScalarField
    (
        const omegaWangNeumannWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    omegaWangNeumannWallFunctionFvPatchScalarField
    (
        const omegaWangNeumannWallFunctionFvPatchScalarField&
    );

    omegaWangNeumannWallFunctionFvPatchScalarField
    (
        const omegaWangNeumannWallFunctionFvPatchScalarField&,
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

    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
