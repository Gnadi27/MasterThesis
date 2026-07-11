#include "nutLiWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{

tmp<scalarField> nutLiWallFunctionFvPatchScalarField::calcNut() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const LiWallFunction::PatchData data =
        LiWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.nutShear);
}


nutLiWallFunctionFvPatchScalarField::nutLiWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    coeffs_()
{}


nutLiWallFunctionFvPatchScalarField::nutLiWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{}


nutLiWallFunctionFvPatchScalarField::nutLiWallFunctionFvPatchScalarField
(
    const nutLiWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{}


nutLiWallFunctionFvPatchScalarField::nutLiWallFunctionFvPatchScalarField
(
    const nutLiWallFunctionFvPatchScalarField& ptf
)
:
    nutWallFunctionFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{}


nutLiWallFunctionFvPatchScalarField::nutLiWallFunctionFvPatchScalarField
(
    const nutLiWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{}


tmp<scalarField> nutLiWallFunctionFvPatchScalarField::yPlus() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const LiWallFunction::PatchData data =
        LiWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.yPlus);
}


void nutLiWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    nutLiWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
