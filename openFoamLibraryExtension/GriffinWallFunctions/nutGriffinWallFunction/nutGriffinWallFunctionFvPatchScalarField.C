#include "nutGriffinWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{

tmp<scalarField> nutGriffinWallFunctionFvPatchScalarField::calcNut() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const GriffinWallFunction::PatchData data =
        GriffinWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.nutShear);
}


nutGriffinWallFunctionFvPatchScalarField::
nutGriffinWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    coeffs_()
{}


nutGriffinWallFunctionFvPatchScalarField::
nutGriffinWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{}


nutGriffinWallFunctionFvPatchScalarField::
nutGriffinWallFunctionFvPatchScalarField
(
    const nutGriffinWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{}


nutGriffinWallFunctionFvPatchScalarField::
nutGriffinWallFunctionFvPatchScalarField
(
    const nutGriffinWallFunctionFvPatchScalarField& ptf
)
:
    nutWallFunctionFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{}


nutGriffinWallFunctionFvPatchScalarField::
nutGriffinWallFunctionFvPatchScalarField
(
    const nutGriffinWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{}


tmp<scalarField> nutGriffinWallFunctionFvPatchScalarField::yPlus() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const GriffinWallFunction::PatchData data =
        GriffinWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.yPlus);
}


void nutGriffinWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    nutGriffinWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
