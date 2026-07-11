#include "omegaWangNeumannWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{

void omegaWangNeumannWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulence,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const auto& turbModel =
        refCast<const compressible::turbulenceModel>(turbulence);

    const WangNeumannWallFunction::PatchData data =
        WangNeumannWallFunction::calcPatch(turbModel, patch, coeffs_);

    const volScalarField& rho = turbModel.rho();
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const labelUList& faceCells = patch.faceCells();

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        const scalar w = cornerWeights[facei];
        const scalar rhoc = max(rho[celli], VSMALL);
        const scalar kc = max(k[celli], VSMALL);

        // omega translation of the fixed nodal epsilon (deviation; not in
        // the papers): omega_P = epsilon_P/(Cmu kP)
        omega0[celli] += w*data.epsP[facei]/(coeffs_.Cmu*kc);

        // With the k-sink beta* omega k and beta* = Cmu, the implicit
        // near-wall sink equals epsilon_P; deliver the cell-averaged
        // production and dissipation exactly as in the epsilon variant.
        G0[celli] +=
            w
           *(
                data.Pk[facei]/rhoc
              + data.epsP[facei]
              - data.epsBar[facei]
            );
    }
}


omegaWangNeumannWallFunctionFvPatchScalarField::
omegaWangNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    omegaWallFunctionFvPatchScalarField(p, iF),
    coeffs_()
{}


omegaWangNeumannWallFunctionFvPatchScalarField::
omegaWangNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    omegaWallFunctionFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{}


omegaWangNeumannWallFunctionFvPatchScalarField::
omegaWangNeumannWallFunctionFvPatchScalarField
(
    const omegaWangNeumannWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    omegaWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{}


omegaWangNeumannWallFunctionFvPatchScalarField::
omegaWangNeumannWallFunctionFvPatchScalarField
(
    const omegaWangNeumannWallFunctionFvPatchScalarField& ptf
)
:
    omegaWallFunctionFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{}


omegaWangNeumannWallFunctionFvPatchScalarField::
omegaWangNeumannWallFunctionFvPatchScalarField
(
    const omegaWangNeumannWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    omegaWallFunctionFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{}


void omegaWangNeumannWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    omegaWallFunctionFvPatchScalarField::write(os);
    coeffs_.writeEntries(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    omegaWangNeumannWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
