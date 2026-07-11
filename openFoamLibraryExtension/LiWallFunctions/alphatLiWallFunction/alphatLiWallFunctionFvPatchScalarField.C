#include "alphatLiWallFunctionFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Time.H"
#include "UPstream.H"
#include "volFields.H"
#include "wallFvPatch.H"

namespace Foam
{
namespace compressible
{
namespace
{

void fieldStats
(
    const scalarField& values,
    scalar& minValue,
    scalar& maxValue,
    scalar& meanValue
)
{
    minValue = GREAT;
    maxValue = -GREAT;
    meanValue = 0;

    if (values.empty())
    {
        minValue = 0;
        maxValue = 0;
        return;
    }

    forAll(values, i)
    {
        minValue = min(minValue, values[i]);
        maxValue = max(maxValue, values[i]);
        meanValue += values[i];
    }

    meanValue /= values.size();
}

} // End anonymous namespace


void alphatLiWallFunctionFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


void alphatLiWallFunctionFvPatchScalarField::checkThermalFluxTarget() const
{
    if
    (
        thermalFluxTarget_ != "auto"
     && thermalFluxTarget_ != "qw"
     && thermalFluxTarget_ != "qHalf"
    )
    {
        FatalErrorInFunction
            << "Unsupported thermalFluxTarget " << thermalFluxTarget_ << nl
            << "Supported values are auto, qw, and qHalf" << nl
            << exit(FatalError);
    }
}


bool alphatLiWallFunctionFvPatchScalarField::adiabatic
(
    const fvPatchScalarField& Tpatch
) const
{
    if (coeffs_.thermalMode == "adiabatic")
    {
        return true;
    }
    if (coeffs_.thermalMode == "isothermal")
    {
        return false;
    }
    return !Tpatch.fixesValue();
}


void alphatLiWallFunctionFvPatchScalarField::writeDiagnostics
(
    const word& selectedTarget,
    const scalarField& alphaEff,
    const scalarField& heSnGrad,
    const LiWallFunction::PatchData& data
)
{
    if (!diagnostics_)
    {
        return;
    }

    const Time& runTime = db().time();
    if (diagnosticTimeIndex_ == runTime.timeIndex())
    {
        return;
    }
    if
    (
        diagnosticWriteInterval_ > 0
     && runTime.timeIndex() % diagnosticWriteInterval_ != 0
    )
    {
        return;
    }
    diagnosticTimeIndex_ = runTime.timeIndex();

    const scalarField& qPaper =
        selectedTarget == "qHalf" ? data.qHalf : data.qw;
    const scalarField qTargetOpenFOAM(-qPaper);
    const scalarField qOpenFOAM(alphaEff*heSnGrad);
    const scalarField qError(qOpenFOAM - qTargetOpenFOAM);

    if (diagnosticLog_)
    {
        scalar qTargetMin, qTargetMax, qTargetMean;
        scalar qFoamMin, qFoamMax, qFoamMean;
        scalar qErrMin, qErrMax, qErrMean;
        scalar uTauMin, uTauMax, uTauMean;

        fieldStats(qTargetOpenFOAM, qTargetMin, qTargetMax, qTargetMean);
        fieldStats(qOpenFOAM, qFoamMin, qFoamMax, qFoamMean);
        fieldStats(qError, qErrMin, qErrMax, qErrMean);
        fieldStats(data.uTau, uTauMin, uTauMax, uTauMean);

        Info<< "Li wall-function diagnostics, patch " << patch().name()
            << ", time " << runTime.timeName() << nl
            << "    thermalFluxTarget = " << selectedTarget << nl
            << "    qTargetOF min/max/mean = "
            << qTargetMin << " / " << qTargetMax << " / " << qTargetMean << nl
            << "    qOpenFOAM min/max/mean = "
            << qFoamMin << " / " << qFoamMax << " / " << qFoamMean << nl
            << "    qError min/max/mean = "
            << qErrMin << " / " << qErrMax << " / " << qErrMean << nl
            << "    uTau min/max/mean = "
            << uTauMin << " / " << uTauMax << " / " << uTauMean
            << endl;
    }

    if (!diagnosticWrite_)
    {
        return;
    }

    fileName dir
    (
        runTime.globalPath()
       /"postProcessing"
       /"LiWallFunctionDiagnostics"
    );
    if (UPstream::parRun())
    {
        dir = dir/("processor" + Foam::name(UPstream::myProcNo()));
    }
    dir = dir/patch().name();
    mkDir(dir);

    OFstream os(dir/(runTime.timeName() + ".dat"));
    os.precision(diagnosticWritePrecision_);
    os  << "# thermalFluxTarget = " << selectedTarget << nl
        << "# face x y z refCell status valid converged iterations residual"
        << " uTau tauw yPlus uPlus rg beta Gamma Tw hw h1 hRef hrg"
        << " qw qHalf uHalf mut1 nuTilda1 nut alphat alphaEff"
        << " heSnGrad qTargetOF qOpenFOAM qError" << nl;

    const vectorField& Cf = patch().Cf();
    const scalarField& alphat = *this;
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField nut(turbModel.nut(patch().index()));

    forAll(qOpenFOAM, facei)
    {
        os  << facei << ' '
            << Cf[facei].x() << ' ' << Cf[facei].y() << ' ' << Cf[facei].z() << ' '
            << data.referenceCell[facei] << ' ' << data.status[facei] << ' '
            << data.valid[facei] << ' ' << data.converged[facei] << ' '
            << data.iterations[facei] << ' ' << data.residual[facei] << ' '
            << data.uTau[facei] << ' ' << data.tauw[facei] << ' '
            << data.yPlus[facei] << ' ' << data.uPlus[facei] << ' '
            << data.rg[facei] << ' ' << data.beta[facei] << ' '
            << data.Gamma[facei] << ' ' << data.Tw[facei] << ' '
            << data.hw[facei] << ' ' << data.h1[facei] << ' '
            << data.hRef[facei] << ' ' << data.hrg[facei] << ' '
            << data.qw[facei] << ' ' << data.qHalf[facei] << ' '
            << data.uHalf[facei] << ' ' << data.mut1[facei] << ' '
            << data.nuTilda1[facei] << ' ' << nut[facei] << ' '
            << alphat[facei] << ' ' << alphaEff[facei] << ' '
            << heSnGrad[facei] << ' ' << qTargetOpenFOAM[facei] << ' '
            << qOpenFOAM[facei] << ' ' << qError[facei] << nl;
    }
}


alphatLiWallFunctionFvPatchScalarField::
alphatLiWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    coeffs_(),
    thermalFluxTarget_("auto"),
    diagnostics_(false),
    diagnosticLog_(true),
    diagnosticWrite_(true),
    diagnosticWriteInterval_(1),
    diagnosticWritePrecision_(12),
    diagnosticTimeIndex_(-1)
{
    checkType();
    checkThermalFluxTarget();
}


alphatLiWallFunctionFvPatchScalarField::
alphatLiWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    coeffs_(dict),
    thermalFluxTarget_(dict.getOrDefault<word>("thermalFluxTarget", "auto")),
    diagnostics_(dict.getOrDefault<bool>("diagnostics", false)),
    diagnosticLog_(dict.getOrDefault<bool>("diagnosticLog", true)),
    diagnosticWrite_(dict.getOrDefault<bool>("diagnosticWrite", true)),
    diagnosticWriteInterval_(dict.getOrDefault<label>("diagnosticWriteInterval", 1)),
    diagnosticWritePrecision_(dict.getOrDefault<label>("diagnosticWritePrecision", 12)),
    diagnosticTimeIndex_(-1)
{
    checkType();
    checkThermalFluxTarget();
}


alphatLiWallFunctionFvPatchScalarField::
alphatLiWallFunctionFvPatchScalarField
(
    const alphatLiWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_),
    thermalFluxTarget_(ptf.thermalFluxTarget_),
    diagnostics_(ptf.diagnostics_),
    diagnosticLog_(ptf.diagnosticLog_),
    diagnosticWrite_(ptf.diagnosticWrite_),
    diagnosticWriteInterval_(ptf.diagnosticWriteInterval_),
    diagnosticWritePrecision_(ptf.diagnosticWritePrecision_),
    diagnosticTimeIndex_(-1)
{
    checkType();
    checkThermalFluxTarget();
}


alphatLiWallFunctionFvPatchScalarField::
alphatLiWallFunctionFvPatchScalarField
(
    const alphatLiWallFunctionFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_),
    thermalFluxTarget_(ptf.thermalFluxTarget_),
    diagnostics_(ptf.diagnostics_),
    diagnosticLog_(ptf.diagnosticLog_),
    diagnosticWrite_(ptf.diagnosticWrite_),
    diagnosticWriteInterval_(ptf.diagnosticWriteInterval_),
    diagnosticWritePrecision_(ptf.diagnosticWritePrecision_),
    diagnosticTimeIndex_(-1)
{
    checkType();
    checkThermalFluxTarget();
}


alphatLiWallFunctionFvPatchScalarField::
alphatLiWallFunctionFvPatchScalarField
(
    const alphatLiWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_),
    thermalFluxTarget_(ptf.thermalFluxTarget_),
    diagnostics_(ptf.diagnostics_),
    diagnosticLog_(ptf.diagnosticLog_),
    diagnosticWrite_(ptf.diagnosticWrite_),
    diagnosticWriteInterval_(ptf.diagnosticWriteInterval_),
    diagnosticWritePrecision_(ptf.diagnosticWritePrecision_),
    diagnosticTimeIndex_(-1)
{
    checkType();
    checkThermalFluxTarget();
}


void alphatLiWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();
    const auto& turbModel = db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const auto& thermo = turbModel.transport();
    const fvPatchScalarField& Tpatch = thermo.T().boundaryField()[patchi];
    const bool isAdiabatic = adiabatic(Tpatch);

    word selectedTarget = thermalFluxTarget_;
    if (selectedTarget == "auto")
    {
        selectedTarget = isAdiabatic ? word("qw") : word("qHalf");
    }
    if (isAdiabatic && selectedTarget == "qHalf")
    {
        FatalErrorInFunction
            << "thermalFluxTarget qHalf is incompatible with adiabatic mode. "
            << "Use auto or qw so the physical boundary heat flux remains zero."
            << nl << exit(FatalError);
    }

    const LiWallFunction::PatchData data =
        LiWallFunction::calcPatch(turbModel, patch(), coeffs_);
    const scalarField& qPaper =
        selectedTarget == "qHalf" ? data.qHalf : data.qw;

    const tmp<scalarField> talpha = turbModel.alpha(patchi);
    const scalarField& alpha = talpha();
    const fvPatchScalarField& hePatch = thermo.he().boundaryField()[patchi];
    const scalarField heSnGrad(hePatch.snGrad());
    const tmp<scalarField> tCpByCpv = thermo.CpByCpv
    (
        thermo.p().boundaryField()[patchi],
        Tpatch,
        patchi
    );
    const scalarField& CpByCpv = tCpByCpv();

    scalarField& alphat = *this;
    forAll(alphat, facei)
    {
        scalar targetAlphaEff = CpByCpv[facei]*alpha[facei];
        if (mag(heSnGrad[facei]) > ROOTVSMALL)
        {
            targetAlphaEff = -qPaper[facei]/heSnGrad[facei];
        }

        alphat[facei] = max
        (
            scalar(0),
            targetAlphaEff/max(CpByCpv[facei], ROOTVSMALL) - alpha[facei]
        );
    }

    const tmp<scalarField> talphaEff = thermo.alphaEff(alphat, patchi);
    writeDiagnostics(selectedTarget, talphaEff(), heSnGrad, data);
    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatLiWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    coeffs_.writeEntries(os);
    os.writeEntryIfDifferent<word>
    (
        "thermalFluxTarget",
        word("auto"),
        thermalFluxTarget_
    );
    os.writeEntryIfDifferent<bool>("diagnostics", false, diagnostics_);
    os.writeEntryIfDifferent<bool>("diagnosticLog", true, diagnosticLog_);
    os.writeEntryIfDifferent<bool>("diagnosticWrite", true, diagnosticWrite_);
    os.writeEntryIfDifferent<label>
    (
        "diagnosticWriteInterval",
        1,
        diagnosticWriteInterval_
    );
    os.writeEntryIfDifferent<label>
    (
        "diagnosticWritePrecision",
        12,
        diagnosticWritePrecision_
    );
    fvPatchField<scalar>::writeValueEntry(os);
}


makePatchTypeField
(
    fvPatchScalarField,
    alphatLiWallFunctionFvPatchScalarField
);

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
