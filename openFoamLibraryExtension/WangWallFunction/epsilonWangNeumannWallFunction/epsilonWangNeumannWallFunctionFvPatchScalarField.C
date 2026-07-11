#include "epsilonWangNeumannWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"

namespace Foam
{

void epsilonWangNeumannWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulence,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const auto& turbModel =
        refCast<const compressible::turbulenceModel>(turbulence);

    const WangNeumannWallFunction::PatchData data =
        WangNeumannWallFunction::calcPatch(turbModel, patch, coeffs_);

    const volScalarField& rho = turbModel.rho();
    const labelUList& faceCells = patch.faceCells();

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        const scalar w = cornerWeights[facei];
        const scalar rhoc = max(rho[celli], VSMALL);

        // Nodal epsilon fixed in the near-wall cell (Craft 2002 Eq. 4.9)
        epsilon0[celli] += w*data.epsP[facei];

        // Cell-averaged production and dissipation of k (WF13/WF16),
        // delivered so that the k equation's net source is Pk/rho - epsBar
        // while its implicit sink uses the fixed nodal epsilon_P.
        G0[celli] +=
            w
           *(
                data.Pk[facei]/rhoc
              + data.epsP[facei]
              - data.epsBar[facei]
            );
    }
}


epsilonWangNeumannWallFunctionFvPatchScalarField::
epsilonWangNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF),
    coeffs_()
{}


epsilonWangNeumannWallFunctionFvPatchScalarField::
epsilonWangNeumannWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF, dict),
    coeffs_(dict)
{}


epsilonWangNeumannWallFunctionFvPatchScalarField::
epsilonWangNeumannWallFunctionFvPatchScalarField
(
    const epsilonWangNeumannWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    epsilonWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_)
{}


epsilonWangNeumannWallFunctionFvPatchScalarField::
epsilonWangNeumannWallFunctionFvPatchScalarField
(
    const epsilonWangNeumannWallFunctionFvPatchScalarField& ptf
)
:
    epsilonWallFunctionFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_)
{}


epsilonWangNeumannWallFunctionFvPatchScalarField::
epsilonWangNeumannWallFunctionFvPatchScalarField
(
    const epsilonWangNeumannWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_)
{}


void epsilonWangNeumannWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    epsilonWallFunctionFvPatchScalarField::write(os);
    coeffs_.writeEntries(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    epsilonWangNeumannWallFunctionFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
