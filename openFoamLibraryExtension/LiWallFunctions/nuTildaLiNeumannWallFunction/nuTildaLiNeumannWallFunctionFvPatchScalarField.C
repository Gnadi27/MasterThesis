#include "nuTildaLiNeumannWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "wallFvPatch.H"

namespace Foam
{

void nuTildaLiNeumannWallFunctionFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


nuTildaLiNeumannWallFunctionFvPatchScalarField::
nuTildaLiNeumannWallFunctionFvPatchScalarField
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


nuTildaLiNeumannWallFunctionFvPatchScalarField::
nuTildaLiNeumannWallFunctionFvPatchScalarField
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


nuTildaLiNeumannWallFunctionFvPatchScalarField::
nuTildaLiNeumannWallFunctionFvPatchScalarField
(
    const nuTildaLiNeumannWallFunctionFvPatchScalarField& ptf,
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


nuTildaLiNeumannWallFunctionFvPatchScalarField::
nuTildaLiNeumannWallFunctionFvPatchScalarField
(
    const nuTildaLiNeumannWallFunctionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


nuTildaLiNeumannWallFunctionFvPatchScalarField::
nuTildaLiNeumannWallFunctionFvPatchScalarField
(
    const nuTildaLiNeumannWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


void nuTildaLiNeumannWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // The SA variable remains zero at the physical wall. Equation (39) is
    // imposed in the first off-wall cell by manipulateMatrix().
    operator==(scalar(0));
    fixedValueFvPatchScalarField::updateCoeffs();
}


void nuTildaLiNeumannWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
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

    const LiNeumannWallFunction::PatchData data =
        LiNeumannWallFunction::calcPatch(turbModel, patch(), coeffs_);

    matrix.setValues(patch().faceCells(), data.nuTilda1);
    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void nuTildaLiNeumannWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    nuTildaLiNeumannWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
