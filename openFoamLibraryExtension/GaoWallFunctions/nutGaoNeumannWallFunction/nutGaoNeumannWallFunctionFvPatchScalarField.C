#include "nutGaoNeumannWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{

tmp<scalarField> nutGaoNeumannWallFunctionFvPatchScalarField::calcNut() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const GaoNeumannWallFunction::PatchData data =
        GaoNeumannWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.nutShear);
}


nutGaoNeumannWallFunctionFvPatchScalarField::
nutGaoNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    coeffs_()
{}


nutGaoNeumannWallFunctionFvPatchScalarField::
nutGaoNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{}


nutGaoNeumannWallFunctionFvPatchScalarField::
nutGaoNeumannWallFunctionFvPatchScalarField
(
    const nutGaoNeumannWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{}


nutGaoNeumannWallFunctionFvPatchScalarField::
nutGaoNeumannWallFunctionFvPatchScalarField
(
    const nutGaoNeumannWallFunctionFvPatchScalarField& ptf
)
:
    nutWallFunctionFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{}


nutGaoNeumannWallFunctionFvPatchScalarField::
nutGaoNeumannWallFunctionFvPatchScalarField
(
    const nutGaoNeumannWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{}


tmp<scalarField> nutGaoNeumannWallFunctionFvPatchScalarField::yPlus() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const GaoNeumannWallFunction::PatchData data =
        GaoNeumannWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.yPlus);
}


void nutGaoNeumannWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    nutGaoNeumannWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
