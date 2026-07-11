#include "omegaGaoNeumannWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{

void omegaGaoNeumannWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel&,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();
    const labelUList& faceCells = patch.faceCells();

    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const GaoNeumannWallFunction::PatchData data =
        GaoNeumannWallFunction::calcPatch(turbModel, patch, coeffs_);

    const scalar Cmu25 = pow025(coeffs_.Cmu);
    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];

        omega0[celli] += cornerWeights[facei]*data.omega[facei];

        if (data.k[facei] > SMALL)
        {
            G0[celli] +=
                cornerWeights[facei]
               *(data.nutShear[facei] + nuw[facei])
               *magGradUw[facei]
               *Cmu25*sqrt(data.k[facei])
               /(coeffs_.kappa*max(y[facei], VSMALL));
        }
    }
}


omegaGaoNeumannWallFunctionFvPatchScalarField::omegaGaoNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    omegaWallFunctionFvPatchScalarField(p, iF),
    coeffs_()
{}


omegaGaoNeumannWallFunctionFvPatchScalarField::omegaGaoNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    omegaWallFunctionFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{}


omegaGaoNeumannWallFunctionFvPatchScalarField::omegaGaoNeumannWallFunctionFvPatchScalarField
(
    const omegaGaoNeumannWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    omegaWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{}


omegaGaoNeumannWallFunctionFvPatchScalarField::omegaGaoNeumannWallFunctionFvPatchScalarField
(
    const omegaGaoNeumannWallFunctionFvPatchScalarField& ptf
)
:
    omegaWallFunctionFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{}


omegaGaoNeumannWallFunctionFvPatchScalarField::omegaGaoNeumannWallFunctionFvPatchScalarField
(
    const omegaGaoNeumannWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    omegaWallFunctionFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{}


void omegaGaoNeumannWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    coeffs_.writeEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    omegaGaoNeumannWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
