#include "alphatGaoWallFunctionFvPatchScalarField.H"

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

}


void alphatGaoWallFunctionFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


void alphatGaoWallFunctionFvPatchScalarField::checkThermalFluxTarget() const
{
    if (thermalFluxTarget_ != "qw" && thermalFluxTarget_ != "qHalf")
    {
        FatalErrorInFunction
            << "Unsupported thermalFluxTarget " << thermalFluxTarget_ << nl
            << "Supported values are: qw, qHalf" << nl
            << exit(FatalError);
    }
}


void alphatGaoWallFunctionFvPatchScalarField::writeDiagnostics
(
    const scalarField& alphaEff,
    const scalarField& heSnGrad,
    const GaoWallFunction::PatchData& data
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

    const scalarField& alphatw = *this;

    const scalarField& qTargetGao =
        thermalFluxTarget_ == "qHalf" ? data.qHalf : data.qw;
    const scalarField qTargetOpenFOAM(-qTargetGao);
    scalarField qOpenFOAM(alphaEff*heSnGrad);
    scalarField qError(qOpenFOAM - qTargetOpenFOAM);
    scalarField qRatio(qOpenFOAM.size(), Zero);

    forAll(qRatio, facei)
    {
        if (mag(qTargetOpenFOAM[facei]) > ROOTVSMALL)
        {
            qRatio[facei] = qOpenFOAM[facei]/qTargetOpenFOAM[facei];
        }
    }

    if (diagnosticLog_)
    {
        scalar qGaoMin, qGaoMax, qGaoMean;
        scalar qWallMin, qWallMax, qWallMean;
        scalar qTargetMin, qTargetMax, qTargetMean;
        scalar qFoamMin, qFoamMax, qFoamMean;
        scalar qErrMin, qErrMax, qErrMean;
        scalar ratioMin, ratioMax, ratioMean;
        scalar alphatMin, alphatMax, alphatMean;
        scalar alphaEffMin, alphaEffMax, alphaEffMean;

        fieldStats(data.qHalf, qGaoMin, qGaoMax, qGaoMean);
        fieldStats(data.qw, qWallMin, qWallMax, qWallMean);
        fieldStats(qTargetOpenFOAM, qTargetMin, qTargetMax, qTargetMean);
        fieldStats(qOpenFOAM, qFoamMin, qFoamMax, qFoamMean);
        fieldStats(qError, qErrMin, qErrMax, qErrMean);
        fieldStats(qRatio, ratioMin, ratioMax, ratioMean);
        fieldStats(alphatw, alphatMin, alphatMax, alphatMean);
        fieldStats(alphaEff, alphaEffMin, alphaEffMax, alphaEffMean);

        Info<< "Gao alphat diagnostics, patch " << patch().name()
            << ", time " << runTime.timeName() << nl
            << "    thermalFluxTarget = " << thermalFluxTarget_ << nl
            << "    qGaoHalf    min/max/mean = "
            << qGaoMin << " / " << qGaoMax << " / " << qGaoMean << nl
            << "    qGaoWall    min/max/mean = "
            << qWallMin << " / " << qWallMax << " / " << qWallMean << nl
            << "    qTargetOF   min/max/mean = "
            << qTargetMin << " / " << qTargetMax << " / " << qTargetMean << nl
            << "    qOpenFOAM   min/max/mean = "
            << qFoamMin << " / " << qFoamMax << " / " << qFoamMean << nl
            << "    qError      min/max/mean = "
            << qErrMin << " / " << qErrMax << " / " << qErrMean << nl
            << "    qRatio      min/max/mean = "
            << ratioMin << " / " << ratioMax << " / " << ratioMean << nl
            << "    alphat      min/max/mean = "
            << alphatMin << " / " << alphatMax << " / " << alphatMean << nl
            << "    alphaEff    min/max/mean = "
            << alphaEffMin << " / " << alphaEffMax << " / " << alphaEffMean
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
       /"GaoWallFunctionDiagnostics"
    );

    if (UPstream::parRun())
    {
        dir = dir/("processor" + Foam::name(UPstream::myProcNo()));
    }

    dir = dir/patch().name();
    mkDir(dir);

    OFstream os(dir/(runTime.timeName() + ".dat"));
    os.precision(diagnosticWritePrecision_);

    os  << "# thermalFluxTarget = " << thermalFluxTarget_ << nl;
    os  << "# qError and qRatio use OpenFOAM target flux = -"
        << thermalFluxTarget_ << nl;
    os  << "# face x y z"
        << " qGaoHalf qOpenFOAM qError qRatio"
        << " alphat alphaEff heSnGrad"
        << " qw uTau beta yPlus Gamma" << nl;

    const vectorField& Cf = patch().Cf();

    forAll(qOpenFOAM, facei)
    {
        os  << facei << ' '
            << Cf[facei].x() << ' '
            << Cf[facei].y() << ' '
            << Cf[facei].z() << ' '
            << data.qHalf[facei] << ' '
            << qOpenFOAM[facei] << ' '
            << qError[facei] << ' '
            << qRatio[facei] << ' '
            << alphatw[facei] << ' '
            << alphaEff[facei] << ' '
            << heSnGrad[facei] << ' '
            << data.qw[facei] << ' '
            << data.uTau[facei] << ' '
            << data.beta[facei] << ' '
            << data.yPlus[facei] << ' '
            << data.Gamma[facei] << nl;
    }
}


alphatGaoWallFunctionFvPatchScalarField::
alphatGaoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    coeffs_(),
    thermalFluxTarget_("qw"),
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


alphatGaoWallFunctionFvPatchScalarField::
alphatGaoWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    coeffs_(dict),
    thermalFluxTarget_(dict.getOrDefault<word>("thermalFluxTarget", "qw")),
    diagnostics_(dict.getOrDefault<bool>("diagnostics", false)),
    diagnosticLog_(dict.getOrDefault<bool>("diagnosticLog", true)),
    diagnosticWrite_(dict.getOrDefault<bool>("diagnosticWrite", true)),
    diagnosticWriteInterval_
    (
        dict.getOrDefault<label>("diagnosticWriteInterval", 1)
    ),
    diagnosticWritePrecision_
    (
        dict.getOrDefault<label>("diagnosticWritePrecision", 12)
    ),
    diagnosticTimeIndex_(-1)
{
    checkType();
    checkThermalFluxTarget();
}


alphatGaoWallFunctionFvPatchScalarField::
alphatGaoWallFunctionFvPatchScalarField
(
    const alphatGaoWallFunctionFvPatchScalarField& ptf,
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


alphatGaoWallFunctionFvPatchScalarField::
alphatGaoWallFunctionFvPatchScalarField
(
    const alphatGaoWallFunctionFvPatchScalarField& ptf
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


alphatGaoWallFunctionFvPatchScalarField::
alphatGaoWallFunctionFvPatchScalarField
(
    const alphatGaoWallFunctionFvPatchScalarField& ptf,
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


void alphatGaoWallFunctionFvPatchScalarField::updateCoeffs()
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

    const GaoWallFunction::PatchData data =
        GaoWallFunction::calcPatch(turbModel, patch(), coeffs_);

    const tmp<scalarField> talphaw = turbModel.alpha(patchi);
    const scalarField& alphaw = talphaw();

    const auto& thermo = turbModel.transport();

    const fvPatchScalarField& hew = thermo.he().boundaryField()[patchi];
    const scalarField heSnGrad(hew.snGrad());

    const tmp<scalarField> tCpByCpv =
        thermo.CpByCpv
        (
            thermo.p().boundaryField()[patchi],
            thermo.T().boundaryField()[patchi],
            patchi
        );
    const scalarField& CpByCpv = tCpByCpv();

    const scalarField& qTargetGao =
        thermalFluxTarget_ == "qHalf" ? data.qHalf : data.qw;

    scalarField& alphatw = *this;

    forAll(alphatw, facei)
    {
        scalar targetAlphaEff = CpByCpv[facei]*alphaw[facei];

        if (mag(heSnGrad[facei]) > ROOTVSMALL)
        {
            targetAlphaEff = -qTargetGao[facei]/heSnGrad[facei];
        }

        alphatw[facei] =
            max
            (
                scalar(0),
                targetAlphaEff/max(CpByCpv[facei], ROOTVSMALL)
              - alphaw[facei]
            );
    }

    const tmp<scalarField> talphaEff =
        turbModel.transport().alphaEff(alphatw, patchi);
    const scalarField& alphaEff = talphaEff();

    writeDiagnostics(alphaEff, heSnGrad, data);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatGaoWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    coeffs_.writeEntries(os);
    os.writeEntryIfDifferent<word>
    (
        "thermalFluxTarget",
        word("qw"),
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
    alphatGaoWallFunctionFvPatchScalarField
);

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
