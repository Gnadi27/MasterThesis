/*---------------------------------------------------------------------------*\
  epsilon boundary condition of the Wang et al. (2023) AWF^Cmp Neumann
  family, for the standard high-Re k-epsilon model (the model the paper
  couples the AWF to).

  Derives from epsilonWallFunction to reuse its master-patch/corner-
  averaging and matrix-constraint machinery, and replaces the near-wall
  cell treatment by the AWF one:

    - the epsilon transport equation is not solved in near-wall cells;
      the nodal value is prescribed as
          epsilon_P = kP^{3/2}/(c_l y_P)          (Craft 2002 Eq. 4.9)
    - the k equation in the near-wall cell uses the cell-averaged
      production Pk (paper Eq. WF13) and dissipation epsBar (WF16).
      OpenFOAM's k equation sinks the *nodal* cell epsilon implicitly, so
      the cell-averaged sink is delivered through the generation field:
          G = Pk/rho + epsilon_P - epsBar
      which makes the net k source Pk/rho - epsBar at convergence.

  Use with kqRWallFunction on k (zero-gradient; Craft: molecular diffusion
  of k to the wall is neglected). All WangNeumann* conditions on a patch
  must carry identical coefficients.
\*---------------------------------------------------------------------------*/

#ifndef epsilonWangNeumannWallFunctionFvPatchScalarField_H
#define epsilonWangNeumannWallFunctionFvPatchScalarField_H

#include "WangNeumannWallFunction.H"
#include "epsilonWallFunctionFvPatchScalarField.H"

namespace Foam
{

class epsilonWangNeumannWallFunctionFvPatchScalarField
:
    public epsilonWallFunctionFvPatchScalarField
{
protected:

    WangNeumannWallFunction::Coeffs coeffs_;

    //- Calculate the epsilon and G contributions (AWF values)
    virtual void calculate
    (
        const turbulenceModel& turbulence,
        const List<scalar>& cornerWeights,
        const fvPatch& patch,
        scalarField& G,
        scalarField& epsilon
    );


public:

    TypeName("epsilonWangNeumannWallFunction");

    epsilonWangNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    epsilonWangNeumannWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    epsilonWangNeumannWallFunctionFvPatchScalarField
    (
        const epsilonWangNeumannWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    epsilonWangNeumannWallFunctionFvPatchScalarField
    (
        const epsilonWangNeumannWallFunctionFvPatchScalarField&
    );

    epsilonWangNeumannWallFunctionFvPatchScalarField
    (
        const epsilonWangNeumannWallFunctionFvPatchScalarField&,
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
