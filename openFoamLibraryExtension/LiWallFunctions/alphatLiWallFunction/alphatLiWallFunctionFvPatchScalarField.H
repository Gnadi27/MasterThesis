#ifndef alphatLiWallFunctionFvPatchScalarField_H
#define alphatLiWallFunctionFvPatchScalarField_H

#include "LiWallFunction.H"
#include "fixedValueFvPatchFields.H"
#include "turbulentFluidThermoModel.H"

namespace Foam
{
namespace compressible
{

class alphatLiWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    LiWallFunction::Coeffs coeffs_;
    word thermalFluxTarget_;
    bool diagnostics_;
    bool diagnosticLog_;
    bool diagnosticWrite_;
    label diagnosticWriteInterval_;
    label diagnosticWritePrecision_;
    label diagnosticTimeIndex_;

    void checkType() const;
    void checkThermalFluxTarget() const;
    bool adiabatic(const fvPatchScalarField& Tpatch) const;

    void writeDiagnostics
    (
        const word& selectedTarget,
        const scalarField& alphaEff,
        const scalarField& heSnGrad,
        const LiWallFunction::PatchData& data
    );

public:
    TypeName("compressible::alphatLiWallFunction");

    alphatLiWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    alphatLiWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    alphatLiWallFunctionFvPatchScalarField
    (
        const alphatLiWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    alphatLiWallFunctionFvPatchScalarField
    (
        const alphatLiWallFunctionFvPatchScalarField&
    );

    alphatLiWallFunctionFvPatchScalarField
    (
        const alphatLiWallFunctionFvPatchScalarField&,
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
    virtual void write(Ostream&) const;
};

} // End namespace compressible
} // End namespace Foam

#endif

// ************************************************************************* //
