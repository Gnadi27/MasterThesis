#include "kGaoWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "wallFvPatch.H"

namespace Foam
{

void kGaoWallFunctionFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


kGaoWallFunctionFvPatchScalarField::kGaoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    coeffs_()
{
    checkType();
}


kGaoWallFunctionFvPatchScalarField::kGaoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{
    checkType();
}


kGaoWallFunctionFvPatchScalarField::kGaoWallFunctionFvPatchScalarField
(
    const kGaoWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


kGaoWallFunctionFvPatchScalarField::kGaoWallFunctionFvPatchScalarField
(
    const kGaoWallFunctionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


kGaoWallFunctionFvPatchScalarField::kGaoWallFunctionFvPatchScalarField
(
    const kGaoWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


void kGaoWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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

    scalarField& kw = *this;

    typedef DimensionedField<scalar, volMesh> FieldType;
    FieldType& kInternal = const_cast<FieldType&>(internalField());

    const labelUList& faceCells = patch().faceCells();

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        kw[facei] = data.k[facei];
        kInternal[celli] = data.k[facei];
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void kGaoWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void kGaoWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    kGaoWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
