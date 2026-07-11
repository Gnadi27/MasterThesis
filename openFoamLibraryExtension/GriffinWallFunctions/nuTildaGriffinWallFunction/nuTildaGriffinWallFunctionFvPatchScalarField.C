#include "nuTildaGriffinWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "wallFvPatch.H"

namespace Foam
{

void nuTildaGriffinWallFunctionFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


nuTildaGriffinWallFunctionFvPatchScalarField::
nuTildaGriffinWallFunctionFvPatchScalarField
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


nuTildaGriffinWallFunctionFvPatchScalarField::
nuTildaGriffinWallFunctionFvPatchScalarField
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


nuTildaGriffinWallFunctionFvPatchScalarField::
nuTildaGriffinWallFunctionFvPatchScalarField
(
    const nuTildaGriffinWallFunctionFvPatchScalarField& ptf,
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


nuTildaGriffinWallFunctionFvPatchScalarField::
nuTildaGriffinWallFunctionFvPatchScalarField
(
    const nuTildaGriffinWallFunctionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


nuTildaGriffinWallFunctionFvPatchScalarField::
nuTildaGriffinWallFunctionFvPatchScalarField
(
    const nuTildaGriffinWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{
    checkType();
}


void nuTildaGriffinWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // The SA variable remains zero at the physical wall. The model eddy
    // viscosity is imposed in the first off-wall cell by
    // manipulateMatrix().
    operator==(scalar(0));
    fixedValueFvPatchScalarField::updateCoeffs();
}


void nuTildaGriffinWallFunctionFvPatchScalarField::manipulateMatrix
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

    const GriffinWallFunction::PatchData data =
        GriffinWallFunction::calcPatch(turbModel, patch(), coeffs_);

    matrix.setValues(patch().faceCells(), data.nuTilda1);
    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void nuTildaGriffinWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    nuTildaGriffinWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
