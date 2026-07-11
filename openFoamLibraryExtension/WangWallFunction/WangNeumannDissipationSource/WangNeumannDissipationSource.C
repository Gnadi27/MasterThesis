#include "WangNeumannDissipationSource.H"

#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"
#include "fvMatrices.H"
#include "turbulentFluidThermoModel.H"
#include "volFields.H"
#include "wallFvPatch.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(WangNeumannDissipationSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        WangNeumannDissipationSource,
        dictionary
    );
}
}


Foam::fv::WangNeumannDissipationSource::WangNeumannDissipationSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    wfCoeffs_(coeffs_),
    patchNames_()
{
    coeffs_.readEntry("patches", patchNames_);

    // Applies to the thermo energy variable (he)
    const auto* thermoPtr =
        mesh_.findObject<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        fieldNames_.resize(1, thermoPtr->he().name());
    }

    if (fieldNames_.empty())
    {
        coeffs().readEntry("fields", fieldNames_);
    }

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "This option applies to exactly one (energy) field, "
            << "settings are: " << fieldNames_ << exit(FatalError);
    }

    fv::option::resetApplied();
}


void Foam::fv::WangNeumannDissipationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    const auto& turbModel =
        mesh_.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    auto tS = volScalarField::New
    (
        name_ + ":S",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimPower/dimVolume, Zero)
    );
    scalarField& S = tS.ref().primitiveFieldRef();

    for (const fvPatch& patch : mesh_.boundary())
    {
        if (!isA<wallFvPatch>(patch) || !patchNames_.match(patch.name()))
        {
            continue;
        }

        const WangNeumannWallFunction::PatchData data =
            WangNeumannWallFunction::calcPatch(turbModel, patch, wfCoeffs_);

        const labelUList& faceCells = patch.faceCells();

        forAll(faceCells, facei)
        {
            // Replace the solver's wall-normal viscous dissipation by the
            // analytical cell average: add the difference.
            S[faceCells[facei]] =
                data.dissAvg[facei] - data.dissLocal[facei];
        }
    }

    eqn += tS;
}


bool Foam::fv::WangNeumannDissipationSource::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        wfCoeffs_ = WangNeumannWallFunction::Coeffs(coeffs_);
        coeffs_.readEntry("patches", patchNames_);
        return true;
    }

    return false;
}

// ************************************************************************* //
