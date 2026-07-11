#include "nutWangNeumannWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{

tmp<scalarField> nutWangNeumannWallFunctionFvPatchScalarField::calcNut() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const WangNeumannWallFunction::PatchData data =
        WangNeumannWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.nutShear);
}


nutWangNeumannWallFunctionFvPatchScalarField::
nutWangNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    coeffs_()
{}


nutWangNeumannWallFunctionFvPatchScalarField::
nutWangNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{}


nutWangNeumannWallFunctionFvPatchScalarField::
nutWangNeumannWallFunctionFvPatchScalarField
(
    const nutWangNeumannWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{}


nutWangNeumannWallFunctionFvPatchScalarField::
nutWangNeumannWallFunctionFvPatchScalarField
(
    const nutWangNeumannWallFunctionFvPatchScalarField& ptf
)
:
    nutWallFunctionFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{}


nutWangNeumannWallFunctionFvPatchScalarField::
nutWangNeumannWallFunctionFvPatchScalarField
(
    const nutWangNeumannWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{}


tmp<scalarField> nutWangNeumannWallFunctionFvPatchScalarField::yPlus() const
{
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const WangNeumannWallFunction::PatchData data =
        WangNeumannWallFunction::calcPatch(turbModel, patch(), coeffs_);

    return tmp<scalarField>::New(data.yStar);
}


void nutWangNeumannWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    nutWangNeumannWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
