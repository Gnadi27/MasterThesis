#include "kGaoNeumannWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "wallFvPatch.H"

namespace Foam
{

void kGaoNeumannWallFunctionFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


kGaoNeumannWallFunctionFvPatchScalarField::kGaoNeumannWallFunctionFvPatchScalarField
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


kGaoNeumannWallFunctionFvPatchScalarField::kGaoNeumannWallFunctionFvPatchScalarField
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


kGaoNeumannWallFunctionFvPatchScalarField::kGaoNeumannWallFunctionFvPatchScalarField
(
    const kGaoNeumannWallFunctionFvPatchScalarField& ptf,
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


kGaoNeumannWallFunctionFvPatchScalarField::kGaoNeumannWallFunctionFvPatchScalarField
(
    const kGaoNeumannWallFunctionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


kGaoNeumannWallFunctionFvPatchScalarField::kGaoNeumannWallFunctionFvPatchScalarField
(
    const kGaoNeumannWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


void kGaoNeumannWallFunctionFvPatchScalarField::updateCoeffs()
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

    const GaoNeumannWallFunction::PatchData data =
        GaoNeumannWallFunction::calcPatch(turbModel, patch(), coeffs_);

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


void kGaoNeumannWallFunctionFvPatchScalarField::manipulateMatrix
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


void kGaoNeumannWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    kGaoNeumannWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
