/*---------------------------------------------------------------------------*\
  Spalart-Allmaras coupling for the Griffin wall function: nuTilda = 0 at
  the wall, and the first off-wall cell is pinned (through the matrix) to
  the value whose nut equals the Griffin eddy viscosity at the first-cell
  centre (PNAS 2021 Eq. 6 in dimensional, constant-stress form). The paper
  is silent on RANS turbulence-variable coupling; this follows the
  validated Li wall-function pattern (flagged deviation).
\*---------------------------------------------------------------------------*/

#ifndef nuTildaGriffinWallFunctionFvPatchScalarField_H
#define nuTildaGriffinWallFunctionFvPatchScalarField_H

#include "GriffinWallFunction.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

class nuTildaGriffinWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    GriffinWallFunction::Coeffs coeffs_;

    void checkType() const;

public:
    TypeName("nuTildaGriffinWallFunction");

    nuTildaGriffinWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    nuTildaGriffinWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    nuTildaGriffinWallFunctionFvPatchScalarField
    (
        const nuTildaGriffinWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    nuTildaGriffinWallFunctionFvPatchScalarField
    (
        const nuTildaGriffinWallFunctionFvPatchScalarField&
    );

    nuTildaGriffinWallFunctionFvPatchScalarField
    (
        const nuTildaGriffinWallFunctionFvPatchScalarField&,
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
