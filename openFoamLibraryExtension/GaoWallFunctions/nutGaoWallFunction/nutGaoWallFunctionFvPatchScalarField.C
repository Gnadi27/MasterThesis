#include "nutGaoWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{

tmp<scalarField> nutGaoWallFunctionFvPatchScalarField::calcNut() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const GaoWallFunction::PatchData data =
        GaoWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.nutShear);
}


nutGaoWallFunctionFvPatchScalarField::nutGaoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    coeffs_()
{}


nutGaoWallFunctionFvPatchScalarField::nutGaoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{}


nutGaoWallFunctionFvPatchScalarField::nutGaoWallFunctionFvPatchScalarField
(
    const nutGaoWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{}


nutGaoWallFunctionFvPatchScalarField::nutGaoWallFunctionFvPatchScalarField
(
    const nutGaoWallFunctionFvPatchScalarField& ptf
)
:
    nutWallFunctionFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{}


nutGaoWallFunctionFvPatchScalarField::nutGaoWallFunctionFvPatchScalarField
(
    const nutGaoWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{}


tmp<scalarField> nutGaoWallFunctionFvPatchScalarField::yPlus() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const GaoWallFunction::PatchData data =
        GaoWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.yPlus);
}


void nutGaoWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    nutGaoWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
