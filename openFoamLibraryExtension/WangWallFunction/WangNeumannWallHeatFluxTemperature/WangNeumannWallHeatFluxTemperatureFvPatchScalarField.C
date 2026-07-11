#include "WangNeumannWallHeatFluxTemperatureFvPatchScalarField.H"

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


void WangNeumannWallHeatFluxTemperatureFvPatchScalarField::checkType() const
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


void WangNeumannWallHeatFluxTemperatureFvPatchScalarField::checkEntries() const
{
    if (!adiabaticIntent())
    {
        if (coeffs_.Tw <= 0)
        {
            FatalErrorInFunction
                << "WangNeumannWallHeatFluxTemperature on patch "
                << patch().name() << " requires a positive Tw entry: the "
                << "patch T value drifts along the prescribed gradient and "
                << "cannot supply the physical wall temperature." << nl
                << exit(FatalError);
        }

        if (coeffs_.wallTransportMode == "fromPatch")
        {
            FatalErrorInFunction
                << "WangNeumannWallHeatFluxTemperature on patch "
                << patch().name() << " requires a wallTransport sub-"
                << "dictionary with mode sutherland or constant: the patch "
                << "transport properties are evaluated at the drifted patch "
                << "T and are not the wall values." << nl << exit(FatalError);
        }
    }

    if (gradientRelaxation_ <= 0 || gradientRelaxation_ > 1)
    {
        FatalErrorInFunction
            << "gradientRelaxation must be in (0, 1]" << nl
            << exit(FatalError);
    }
}


void WangNeumannWallHeatFluxTemperatureFvPatchScalarField::writeDiagnostics
(
    const scalarField& alphaEff,
    const scalarField& heGradient,
    const scalarField& TFace,
    const WangNeumannWallFunction::PatchData& data
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

    const scalarField qTargetOpenFOAM(data.qw);
    const scalarField qOpenFOAM(alphaEff*heGradient);
    const scalarField qError(qOpenFOAM - qTargetOpenFOAM);

    if (diagnosticLog_)
    {
        scalar qTargetMin, qTargetMax, qTargetMean;
        scalar qFoamMin, qFoamMax, qFoamMean;
        scalar qErrMin, qErrMax, qErrMean;
        scalar tauwMin, tauwMax, tauwMean;
        scalar TFaceMin, TFaceMax, TFaceMean;

        fieldStats(qTargetOpenFOAM, qTargetMin, qTargetMax, qTargetMean);
        fieldStats(qOpenFOAM, qFoamMin, qFoamMax, qFoamMean);
        fieldStats(qError, qErrMin, qErrMax, qErrMean);
        fieldStats(data.tauw, tauwMin, tauwMax, tauwMean);
        fieldStats(TFace, TFaceMin, TFaceMax, TFaceMean);

        Info<< "Wang Neumann wall-function diagnostics, patch "
            << patch().name()
            << ", time " << runTime.timeName() << nl
            << "    qTargetOF min/max/mean = "
            << qTargetMin << " / " << qTargetMax << " / " << qTargetMean << nl
            << "    qOpenFOAM min/max/mean = "
            << qFoamMin << " / " << qFoamMax << " / " << qFoamMean << nl
            << "    qError min/max/mean = "
            << qErrMin << " / " << qErrMax << " / " << qErrMean << nl
            << "    TFace min/max/mean = "
            << TFaceMin << " / " << TFaceMax << " / " << TFaceMean << nl
            << "    tauw min/max/mean = "
            << tauwMin << " / " << tauwMax << " / " << tauwMean
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
       /"WangNeumannWallFunctionDiagnostics"
    );
    if (UPstream::parRun())
    {
        dir = dir/("processor" + Foam::name(UPstream::myProcNo()));
    }
    dir = dir/patch().name();
    mkDir(dir);

    OFstream os(dir/(runTime.timeName() + ".dat"));
    os.precision(diagnosticWritePrecision_);
    os  << "# face x y z thinCell tauw A1 yStar yStarN qw Tw Tv Ath1"
        << " Pk epsBar epsP dissAvg muw muv rhov cpw Pr bMu D C Dth"
        << " Un Tn yn nutShear"
        << " gradT TFace alphaEff heGradient qTargetOF qOpenFOAM qError"
        << nl;

    const vectorField& Cf = patch().Cf();
    const scalarField& gradT = gradient();

    forAll(qOpenFOAM, facei)
    {
        os  << facei << ' '
            << Cf[facei].x() << ' ' << Cf[facei].y() << ' '
            << Cf[facei].z() << ' '
            << data.thinCell[facei] << ' '
            << data.tauw[facei] << ' ' << data.A1[facei] << ' '
            << data.yStar[facei] << ' ' << data.yStarN[facei] << ' '
            << data.qw[facei] << ' ' << data.Tw[facei] << ' '
            << data.Tv[facei] << ' ' << data.Ath1[facei] << ' '
            << data.Pk[facei] << ' ' << data.epsBar[facei] << ' '
            << data.epsP[facei] << ' ' << data.dissAvg[facei] << ' '
            << data.muw[facei] << ' ' << data.muv[facei] << ' '
            << data.rhov[facei] << ' ' << data.cpw[facei] << ' '
            << data.Pr[facei] << ' ' << data.bMu[facei] << ' '
            << data.Dcoef[facei] << ' ' << data.Ccoef[facei] << ' '
            << data.DthCoef[facei] << ' '
            << data.Un[facei] << ' ' << data.Tn[facei] << ' '
            << data.yn[facei] << ' ' << data.nutShear[facei] << ' '
            << gradT[facei] << ' ' << TFace[facei] << ' '
            << alphaEff[facei] << ' ' << heGradient[facei] << ' '
            << qTargetOpenFOAM[facei] << ' '
            << qOpenFOAM[facei] << ' ' << qError[facei] << nl;
    }
}


WangNeumannWallHeatFluxTemperatureFvPatchScalarField::
WangNeumannWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    coeffs_(),
    gradientRelaxation_(1),
    TFaceMin_(1),
    TFaceMax_(GREAT),
    diagnostics_(false),
    diagnosticLog_(true),
    diagnosticWrite_(true),
    diagnosticWriteInterval_(1),
    diagnosticWritePrecision_(12),
    diagnosticTimeIndex_(-1)
{
    checkType();
}


WangNeumannWallHeatFluxTemperatureFvPatchScalarField::
WangNeumannWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    // Optional gradient entry: first start has no converged q_w, so an
    // absent gradient begins from the adiabatic (zero-gradient) state.
    fixedGradientFvPatchScalarField
    (
        p,
        iF,
        dict,
        IOobjectOption::READ_IF_PRESENT
    ),
    coeffs_(dict),
    gradientRelaxation_(dict.getOrDefault<scalar>("gradientRelaxation", 1)),
    TFaceMin_(dict.getOrDefault<scalar>("TFaceMin", 1)),
    TFaceMax_(dict.getOrDefault<scalar>("TFaceMax", GREAT)),
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
    checkEntries();
}


WangNeumannWallHeatFluxTemperatureFvPatchScalarField::
WangNeumannWallHeatFluxTemperatureFvPatchScalarField
(
    const WangNeumannWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    coeffs_(ptf.coeffs_),
    gradientRelaxation_(ptf.gradientRelaxation_),
    TFaceMin_(ptf.TFaceMin_),
    TFaceMax_(ptf.TFaceMax_),
    diagnostics_(ptf.diagnostics_),
    diagnosticLog_(ptf.diagnosticLog_),
    diagnosticWrite_(ptf.diagnosticWrite_),
    diagnosticWriteInterval_(ptf.diagnosticWriteInterval_),
    diagnosticWritePrecision_(ptf.diagnosticWritePrecision_),
    diagnosticTimeIndex_(-1)
{
    checkType();
}


WangNeumannWallHeatFluxTemperatureFvPatchScalarField::
WangNeumannWallHeatFluxTemperatureFvPatchScalarField
(
    const WangNeumannWallHeatFluxTemperatureFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    coeffs_(ptf.coeffs_),
    gradientRelaxation_(ptf.gradientRelaxation_),
    TFaceMin_(ptf.TFaceMin_),
    TFaceMax_(ptf.TFaceMax_),
    diagnostics_(ptf.diagnostics_),
    diagnosticLog_(ptf.diagnosticLog_),
    diagnosticWrite_(ptf.diagnosticWrite_),
    diagnosticWriteInterval_(ptf.diagnosticWriteInterval_),
    diagnosticWritePrecision_(ptf.diagnosticWritePrecision_),
    diagnosticTimeIndex_(-1)
{
    checkType();
}


WangNeumannWallHeatFluxTemperatureFvPatchScalarField::
WangNeumannWallHeatFluxTemperatureFvPatchScalarField
(
    const WangNeumannWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    coeffs_(ptf.coeffs_),
    gradientRelaxation_(ptf.gradientRelaxation_),
    TFaceMin_(ptf.TFaceMin_),
    TFaceMax_(ptf.TFaceMax_),
    diagnostics_(ptf.diagnostics_),
    diagnosticLog_(ptf.diagnosticLog_),
    diagnosticWrite_(ptf.diagnosticWrite_),
    diagnosticWriteInterval_(ptf.diagnosticWriteInterval_),
    diagnosticWritePrecision_(ptf.diagnosticWritePrecision_),
    diagnosticTimeIndex_(-1)
{
    checkType();
}


void WangNeumannWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
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

    const bool isAdiabatic =
        WangNeumannWallFunction::isAdiabatic(coeffs_, *this);

    if (isAdiabatic && !diagnostics_)
    {
        // Reduces to zeroGradient behaviour identically.
        gradient() = 0;
        fixedGradientFvPatchScalarField::updateCoeffs();
        return;
    }

    const WangNeumannWallFunction::PatchData data =
        WangNeumannWallFunction::calcPatch(turbModel, patch(), coeffs_);

    // Sign mapping: the paper's q_w (WF43, -lambda dT/dy at the wall) is
    // positive when the wall heats the fluid; OpenFOAM's discrete wall
    // diffusion flux alphaEff*snGrad(he) is positive in exactly the same
    // situation (face enthalpy above the cell value). The conventions
    // coincide, so the target is +q_w. (The Gao library's -q_w mapping
    // reflects Nichols' opposite q_w convention and must NOT be copied.)
    const scalarField qTargetOF(data.qw);

    // Effective diffusivity of the energy assembly (CpByCpv*(alpha+alphat));
    // it appears once in this inversion and once in the discrete flux, so
    // its patch-state evaluation cancels from the delivered flux.
    const tmp<scalarField> talphaEff = turbModel.alphaEff(patchi);
    const scalarField& alphaEff = talphaEff();

    const scalarField& deltaC = patch().deltaCoeffs();
    const volScalarField& heVolField = thermo.he();
    const fvPatchScalarField& hePatch = heVolField.boundaryField()[patchi];
    const scalarField heInternal(hePatch.patchInternalField());
    const scalarField pWall(thermo.p().boundaryField()[patchi]);
    const scalarField TInternal(patchInternalField());

    // The discrete wall diffusion flux is alphaEff*deltaC*(he_face - he_1)
    // with he_face evaluated at the extrapolated face temperature. Place
    // the face enthalpy so that this equals qTargetOF exactly, then back
    // out the face temperature and the gradient that extrapolates to it.
    scalarField heFaceNeeded(heInternal);
    forAll(heFaceNeeded, facei)
    {
        heFaceNeeded[facei] +=
            qTargetOF[facei]
           /(max(alphaEff[facei], VSMALL)*max(deltaC[facei], VSMALL));
    }

    const tmp<scalarField> tTFace =
        thermo.THE(heFaceNeeded, pWall, *this, patchi);
    scalarField TFace(tTFace());
    forAll(TFace, facei)
    {
        TFace[facei] = min(max(TFace[facei], TFaceMin_), TFaceMax_);
    }

    scalarField& gradT = gradient();
    forAll(gradT, facei)
    {
        const scalar gNew =
            isAdiabatic
          ? scalar(0)
          : (TFace[facei] - TInternal[facei])*deltaC[facei];

        gradT[facei] =
            gradientRelaxation_*gNew
          + (1 - gradientRelaxation_)*gradT[facei];
    }

    // Achieved face state after evaluate(); equals the target unless the
    // gradient was clamped or relaxed.
    scalarField TFaceAchieved(TInternal);
    forAll(TFaceAchieved, facei)
    {
        TFaceAchieved[facei] += gradT[facei]/max(deltaC[facei], VSMALL);
    }

    const tmp<scalarField> theFace =
        thermo.he(pWall, TFaceAchieved, patchi);
    const scalarField heGradient((theFace() - heInternal)*deltaC);

    // Mirror the delivered energy gradient onto the he boundary condition
    // (gradientEnergy, paired by heBoundaryTypes) so that consumers of
    // snGrad(he) -- the wallHeatFlux function object in particular --
    // report the imposed flux even when the solver never re-evaluates the
    // he boundary condition itself.
    {
        volScalarField& heRef =
            db().lookupObjectRef<volScalarField>(heVolField.name());
        fvPatchScalarField& hePatchRef =
            heRef.boundaryFieldRef()[patchi];
        if (isA<fixedGradientFvPatchScalarField>(hePatchRef))
        {
            refCast<fixedGradientFvPatchScalarField>(hePatchRef).gradient() =
                heGradient;
        }
    }

    writeDiagnostics(alphaEff, heGradient, TFaceAchieved, data);
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void WangNeumannWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    coeffs_.writeEntries(os);
    os.writeEntryIfDifferent<scalar>
    (
        "gradientRelaxation",
        1,
        gradientRelaxation_
    );
    os.writeEntryIfDifferent<scalar>("TFaceMin", 1, TFaceMin_);
    os.writeEntryIfDifferent<scalar>("TFaceMax", GREAT, TFaceMax_);
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
    WangNeumannWallHeatFluxTemperatureFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
