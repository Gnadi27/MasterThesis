#include "GriffinWallFunction.H"
#include "GriffinBoundaryLayerEdge.H"
#include "GriffinWallHeatFluxTemperatureFvPatchScalarField.H"

#include "dictionary.H"
#include "fvMesh.H"
#include "fvPatchField.H"
#include "Ostream.H"
#include "volFields.H"

#include <cmath>

namespace Foam
{
namespace
{

inline scalar positive(const scalar value)
{
    return max(value, VSMALL);
}


inline vector tangential(const vector& value, const vector& normal)
{
    return value - normal*(value & normal);
}

// Calibration of the incompressible law of the wall on the Lee & Moser
// (2015) channel DNS at Re_tau = {182, 543, 1000, 1990, 5190} (paper
// Sec. 2). Fixed by the model, not user-adjustable.
constexpr label nCalibration = 5;
const scalar calibrationReTau[nCalibration] = {182, 543, 1000, 1990, 5190};
const scalar calibrationKappa[nCalibration] =
    {0.400, 0.408, 0.400, 0.391, 0.391};
const scalar calibrationAPlus[nCalibration] = {18.2, 17.4, 17.0, 16.5, 16.5};

} // End anonymous namespace


namespace GriffinWallFunction
{

Coeffs::Coeffs()
:
    s(1.14),
    Pr(0.72),
    cv1(7.1),
    maxIter(60),
    tolerance(1e-6),
    integrationNodes(48),
    subIterations(5),
    thermalMode("auto"),
    strict(false),
    Tw(-1),
    viscosityMode("sutherland"),
    viscosityAs(1.458e-06),
    viscosityTs(110.4),
    viscosityMuRef(0),
    viscosityTRef(0),
    viscosityN(0),
    viscosityMu(0),
    edgeMode("reconstruct"),
    Ue(-1),
    Te(-1),
    delta(-1),
    maxEdgeCells(200),
    maxEdgeDistance(GREAT),
    edgeUpdateInterval(1),
    edgeThreshold(0.99)
{}


Coeffs::Coeffs(const dictionary& dict)
:
    s(dict.getOrDefault<scalar>("s", 1.14)),
    Pr(dict.getOrDefault<scalar>("Pr", 0.72)),
    cv1(dict.getOrDefault<scalar>("cv1", 7.1)),
    maxIter(dict.getOrDefault<label>("maxIter", 60)),
    tolerance(dict.getOrDefault<scalar>("tolerance", 1e-6)),
    integrationNodes(dict.getOrDefault<label>("integrationNodes", 48)),
    subIterations(dict.getOrDefault<label>("subIterations", 5)),
    thermalMode(dict.getOrDefault<word>("thermalMode", "auto")),
    strict(dict.getOrDefault<bool>("strict", false)),
    Tw(dict.getOrDefault<scalar>("Tw", -1)),
    viscosityMode("sutherland"),
    viscosityAs(1.458e-06),
    viscosityTs(110.4),
    viscosityMuRef(0),
    viscosityTRef(0),
    viscosityN(0),
    viscosityMu(0),
    edgeMode(dict.getOrDefault<word>("edgeMode", "reconstruct")),
    Ue(dict.getOrDefault<scalar>("Ue", -1)),
    Te(dict.getOrDefault<scalar>("Te", -1)),
    delta(dict.getOrDefault<scalar>("delta", -1)),
    maxEdgeCells(dict.getOrDefault<label>("maxEdgeCells", 200)),
    maxEdgeDistance(dict.getOrDefault<scalar>("maxEdgeDistance", GREAT)),
    edgeUpdateInterval(dict.getOrDefault<label>("edgeUpdateInterval", 1)),
    edgeThreshold(dict.getOrDefault<scalar>("edgeThreshold", 0.99))
{
    if
    (
        thermalMode != "auto"
     && thermalMode != "isothermal"
     && thermalMode != "adiabatic"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unsupported thermalMode " << thermalMode << nl
            << "Supported values are auto, isothermal, and adiabatic"
            << exit(FatalIOError);
    }

    const dictionary* visc = dict.findDict("viscosity");
    if (visc)
    {
        viscosityMode = visc->getOrDefault<word>("mode", "sutherland");

        if (viscosityMode == "sutherland")
        {
            viscosityAs = visc->getOrDefault<scalar>("As", 1.458e-06);
            viscosityTs = visc->getOrDefault<scalar>("Ts", 110.4);
        }
        else if (viscosityMode == "powerLaw")
        {
            viscosityMuRef = visc->get<scalar>("muRef");
            viscosityTRef = visc->get<scalar>("TRef");
            viscosityN = visc->get<scalar>("n");
        }
        else if (viscosityMode == "constant")
        {
            viscosityMu = visc->get<scalar>("mu");
        }
        else
        {
            FatalIOErrorInFunction(*visc)
                << "Unsupported viscosity mode " << viscosityMode << nl
                << "Supported values are sutherland, powerLaw, and constant"
                << exit(FatalIOError);
        }
    }

    if (edgeMode == "prescribed")
    {
        if (Ue <= 0 || Te <= 0 || delta <= 0)
        {
            FatalIOErrorInFunction(dict)
                << "edgeMode prescribed requires positive Ue, Te, and delta "
                << "entries" << exit(FatalIOError);
        }
    }
    else if (edgeMode != "reconstruct")
    {
        FatalIOErrorInFunction(dict)
            << "Unsupported edgeMode " << edgeMode << nl
            << "Supported values are reconstruct and prescribed"
            << exit(FatalIOError);
    }

    if (edgeThreshold <= 0 || edgeThreshold >= 1)
    {
        FatalIOErrorInFunction(dict)
            << "edgeThreshold must be in (0, 1)" << exit(FatalIOError);
    }
}


void Coeffs::writeEntries(Ostream& os) const
{
    os.writeEntryIfDifferent<scalar>("s", 1.14, s);
    os.writeEntryIfDifferent<scalar>("Pr", 0.72, Pr);
    os.writeEntryIfDifferent<scalar>("cv1", 7.1, cv1);
    os.writeEntryIfDifferent<label>("maxIter", 60, maxIter);
    os.writeEntryIfDifferent<scalar>("tolerance", 1e-6, tolerance);
    os.writeEntryIfDifferent<label>("integrationNodes", 48, integrationNodes);
    os.writeEntryIfDifferent<label>("subIterations", 5, subIterations);
    os.writeEntryIfDifferent<word>("thermalMode", word("auto"), thermalMode);
    os.writeEntryIfDifferent<bool>("strict", false, strict);
    os.writeEntryIfDifferent<scalar>("Tw", -1, Tw);

    os.beginBlock("viscosity");
    os.writeEntry("mode", viscosityMode);
    if (viscosityMode == "sutherland")
    {
        os.writeEntry("As", viscosityAs);
        os.writeEntry("Ts", viscosityTs);
    }
    else if (viscosityMode == "powerLaw")
    {
        os.writeEntry("muRef", viscosityMuRef);
        os.writeEntry("TRef", viscosityTRef);
        os.writeEntry("n", viscosityN);
    }
    else if (viscosityMode == "constant")
    {
        os.writeEntry("mu", viscosityMu);
    }
    os.endBlock();

    os.writeEntryIfDifferent<word>("edgeMode", word("reconstruct"), edgeMode);
    if (edgeMode == "prescribed")
    {
        os.writeEntry("Ue", Ue);
        os.writeEntry("Te", Te);
        os.writeEntry("delta", delta);
    }
    os.writeEntryIfDifferent<label>("maxEdgeCells", 200, maxEdgeCells);
    os.writeEntryIfDifferent<scalar>("maxEdgeDistance", GREAT, maxEdgeDistance);
    os.writeEntryIfDifferent<label>("edgeUpdateInterval", 1, edgeUpdateInterval);
    os.writeEntryIfDifferent<scalar>("edgeThreshold", 0.99, edgeThreshold);
}


scalar Coeffs::mu(const scalar T) const
{
    const scalar TSafe = positive(T);

    if (viscosityMode == "powerLaw")
    {
        return viscosityMuRef*pow(TSafe/positive(viscosityTRef), viscosityN);
    }
    else if (viscosityMode == "constant")
    {
        return viscosityMu;
    }

    return viscosityAs*TSafe*sqrt(TSafe)/(TSafe + viscosityTs);
}


scalar Coeffs::dmudT(const scalar T) const
{
    const scalar TSafe = positive(T);

    if (viscosityMode == "powerLaw")
    {
        return viscosityN*mu(TSafe)/TSafe;
    }
    else if (viscosityMode == "constant")
    {
        return 0;
    }

    return
        viscosityAs*sqrt(TSafe)*(0.5*TSafe + 1.5*viscosityTs)
       /sqr(TSafe + viscosityTs);
}


PatchData::PatchData(const label size)
:
    uTau(size, Zero),
    tauw(size, Zero),
    yPlus(size, Zero),
    yStar(size, Zero),
    uPlus(size, Zero),
    ReTauStar(size, Zero),
    kappaCal(size, Zero),
    APlusCal(size, Zero),
    Tw(size, Zero),
    Tr(size, Zero),
    T1(size, Zero),
    U1(size, Zero),
    Ue(size, Zero),
    Te(size, Zero),
    delta(size, Zero),
    qw(size, Zero),
    nutShear(size, Zero),
    mut1(size, Zero),
    nuTilda1(size, Zero),
    residual(size, Zero),
    iterations(size, Zero),
    converged(size, Zero),
    valid(size, Zero),
    edgeStatus(size, Zero),
    status(size, Zero)
{}


void calibration(const scalar ReTauStar, scalar& kappa, scalar& APlus)
{
    if (ReTauStar <= calibrationReTau[0])
    {
        kappa = calibrationKappa[0];
        APlus = calibrationAPlus[0];
        return;
    }

    if (ReTauStar >= calibrationReTau[nCalibration - 1])
    {
        kappa = calibrationKappa[nCalibration - 1];
        APlus = calibrationAPlus[nCalibration - 1];
        return;
    }

    for (label i = 1; i < nCalibration; ++i)
    {
        if (ReTauStar <= calibrationReTau[i])
        {
            const scalar w =
                (ReTauStar - calibrationReTau[i - 1])
               /(calibrationReTau[i] - calibrationReTau[i - 1]);
            kappa = (1 - w)*calibrationKappa[i - 1] + w*calibrationKappa[i];
            APlus = (1 - w)*calibrationAPlus[i - 1] + w*calibrationAPlus[i];
            return;
        }
    }
}


scalar mutPlusIncompressible
(
    const scalar yStar,
    const scalar kappa,
    const scalar APlus
)
{
    if (yStar <= 0)
    {
        return 0;
    }

    const scalar damping = 1 - exp(-yStar/positive(APlus));
    return kappa*yStar*sqr(damping);
}


scalar temperatureFromU
(
    const scalar U,
    const scalar UM,
    const scalar TM,
    const scalar Ue,
    const scalar Tw,
    const scalar Tr,
    const scalar s,
    const scalar Pr
)
{
    const scalar xe = U/positive(Ue);
    const scalar xm = U/positive(UM);

    return Tw + s*Pr*(Tr - Tw)*xe*(1 - xm) + sqr(xm)*(TM - Tw);
}


scalar dTemperaturedU
(
    const scalar U,
    const scalar UM,
    const scalar TM,
    const scalar Ue,
    const scalar Tw,
    const scalar Tr,
    const scalar s,
    const scalar Pr
)
{
    const scalar UeSafe = positive(Ue);
    const scalar UMSafe = positive(UM);

    return
        s*Pr*(Tr - Tw)*(1/UeSafe - 2*U/(UeSafe*UMSafe))
      + 2*U/sqr(UMSafe)*(TM - Tw);
}


scalar recoveryTemperature
(
    const scalar Te,
    const scalar Ue,
    const scalar Cp,
    const scalar Pr
)
{
    return Te + cbrt(Pr)*sqr(Ue)/(2*positive(Cp));
}


scalar solveNuTilda
(
    const scalar nutTarget,
    const scalar nu,
    const scalar cv1,
    const label maxIter,
    const scalar tolerance
)
{
    if (nutTarget <= SMALL)
    {
        return 0;
    }

    const scalar nuEff = positive(nu);
    const auto residual = [&](const scalar nuTilda) -> scalar
    {
        const scalar chi = nuTilda/nuEff;
        const scalar chi3 = pow3(chi);
        const scalar fv1 = chi3/(chi3 + pow3(cv1));
        return nuTilda*fv1 - nutTarget;
    };

    scalar lo = 0;
    scalar hi = max(nutTarget, nuEff);

    for (label iter = 0; iter < 80 && residual(hi) < 0; ++iter)
    {
        hi *= 2;
    }

    for (label iter = 0; iter < maxIter; ++iter)
    {
        const scalar mid = 0.5*(lo + hi);

        if (residual(mid) < 0)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }

        if (mag(hi - lo)/positive(mid) < tolerance)
        {
            break;
        }
    }

    return 0.5*(lo + hi);
}


InnerResult solveInner(const InnerInput& in, const Coeffs& coeffs)
{
    InnerResult out;
    out.uTau = 0;
    out.tauw = 0;
    out.qw = 0;
    out.yPlus = 0;
    out.yStar = 0;
    out.ReTauStar = 0;
    out.kappa = 0;
    out.APlus = 0;
    out.mut1 = 0;
    out.residual = GREAT;
    out.iterations = 0;
    out.converged = 0;

    const scalar s = coeffs.s;
    const scalar Pr = coeffs.Pr;
    const scalar Cp = positive(in.Cp);
    const scalar Ue = positive(in.Ue);

    out.Tr = recoveryTemperature(in.Te, in.Ue, Cp, Pr);
    const scalar Tw = in.adiabatic ? out.Tr : in.Tw;
    out.Tw = Tw;

    // Thin-boundary-layer assumption: p constant across the inner layer,
    // so the isobaric ideal gas gives rho(T) = rho_M T_M/T.
    const scalar rhoTM = in.rhoM*in.TM;
    const scalar rhoW = rhoTM/positive(Tw);
    const scalar muW = positive(coeffs.mu(Tw));
    const scalar rhoE = rhoTM/positive(in.Te);
    const scalar muE = positive(coeffs.mu(in.Te));

    // Clamp envelope for the algebraic T(U) along the profile.
    const scalar TLo =
        0.2*min(min(positive(Tw), positive(in.TM)), positive(in.Te));
    const scalar THi =
        5.0*max(max(Tw, in.TM), max(in.Te, out.Tr));

    scalar kappa = 0.41;
    scalar APlus = 17;

    // Integrate the inverse-transformation ODE (paper Eq. 2.5) from the
    // wall to the matching point for a given tau_w; returns U(yM).
    const auto integrate =
        [&](const scalar tauw, scalar& dUdyM, scalar& muM) -> scalar
    {
        const scalar ReTauStar =
            sqrt(max(rhoE*tauw, VSMALL))*in.delta/muE;
        calibration(ReTauStar, kappa, APlus);
        out.ReTauStar = ReTauStar;

        const scalar sqrtRhoWTauw = max(sqrt(rhoW*tauw), VSMALL);

        // dU/dy at (y, U) with per-node property sub-iteration: the
        // property-gradient terms of Eq. 2.5 need dT/dy = dT/dU * dU/dy,
        // which depends on the dU/dy being computed (paper Sec. 2.3).
        const auto rhs =
            [&](const scalar yEval, const scalar UEval, const scalar seed)
            -> scalar
        {
            scalar T =
                temperatureFromU(UEval, in.UM, in.TM, Ue, Tw, out.Tr, s, Pr);
            T = min(max(T, TLo), THi);
            const scalar rho = rhoTM/T;
            const scalar muT = positive(coeffs.mu(T));
            const scalar dTdU =
                dTemperaturedU(UEval, in.UM, in.TM, Ue, Tw, out.Tr, s, Pr);
            const scalar dmudTT = coeffs.dmudT(T);
            const scalar yStarEval = yEval*sqrt(rho*max(tauw, VSMALL))/muT;
            const scalar St =
                1/(1 + mutPlusIncompressible(yStarEval, kappa, APlus));
            const scalar viscousTerm = (muW/muT)*(1/St - 1);
            const scalar sqrtRhoPlus = sqrt(rho/rhoW);
            const scalar prefactor = (tauw/muT)*sqrtRhoPlus;

            scalar dUdy = seed;
            for (label k = 0; k < max(coeffs.subIterations, label(1)); ++k)
            {
                const scalar dTdy = dTdU*dUdy;
                // Isobaric ideal gas: (1/rho) drho/dy = -(1/T) dT/dy.
                const scalar bracket =
                    yEval*(-0.5*dTdy/T - dmudTT*dTdy/muT);
                const scalar Phi =
                    min(max(1 + bracket, scalar(0.05)), scalar(20));
                const scalar D = viscousTerm + sqrtRhoPlus*Phi;
                const scalar dUdyNew = prefactor*Phi/max(D, VSMALL);

                const scalar change =
                    mag(dUdyNew - dUdy)/max(mag(dUdyNew), VSMALL);
                dUdy = dUdyNew;
                if (change < coeffs.tolerance)
                {
                    break;
                }
            }

            return dUdy;
        };

        // Viscous seed point: inside y+ <= 0.05 the profile is linear to
        // round-off; start there and march on a log-spaced grid.
        scalar y = min(1e-3*in.yM, 0.05*muW/sqrtRhoWTauw);
        y = max(y, 1e-12*in.yM);
        scalar U = tauw*y/muW;
        scalar dUdy = tauw/muW;

        const label N = max(coeffs.integrationNodes, label(8));
        const scalar logRatio = log(in.yM/y)/N;

        const scalar y0 = y;
        for (label i = 1; i <= N; ++i)
        {
            const scalar yNext = (i == N) ? in.yM : y0*exp(i*logRatio);
            const scalar hStep = yNext - y;
            const scalar k1 = rhs(y, U, dUdy);
            const scalar k2 = rhs(y + 0.5*hStep, U + 0.5*hStep*k1, k1);
            U += hStep*k2;
            y = yNext;
            dUdy = k2;
        }

        dUdyM = rhs(in.yM, U, dUdy);

        scalar TAtM = temperatureFromU(U, in.UM, in.TM, Ue, Tw, out.Tr, s, Pr);
        TAtM = min(max(TAtM, TLo), THi);
        muM = positive(coeffs.mu(TAtM));

        return U;
    };

    scalar dUdyM = 0;
    scalar muM = muW;

    const auto shootResidual = [&](const scalar tauw) -> scalar
    {
        return integrate(tauw, dUdyM, muM) - in.UM;
    };

    // Bracket around the initial guess: U(yM; tau_w) increases
    // monotonically with tau_w.
    scalar lo = max(in.tauwGuess, SMALL);
    scalar hi = lo;
    scalar fLo = shootResidual(lo);
    scalar fHi = fLo;

    if (fLo > 0)
    {
        for (label iter = 0; iter < 60 && fLo > 0; ++iter)
        {
            hi = lo;
            fHi = fLo;
            lo *= 0.25;
            fLo = shootResidual(lo);
        }
    }
    else
    {
        for (label iter = 0; iter < 60 && fHi < 0; ++iter)
        {
            lo = hi;
            fLo = fHi;
            hi *= 4;
            fHi = shootResidual(hi);
        }
    }

    if (!std::isfinite(fLo) || !std::isfinite(fHi) || fLo > 0 || fHi < 0)
    {
        out.residual = min(mag(fLo), mag(fHi))/positive(in.UM);
        out.tauw = mag(fLo) < mag(fHi) ? lo : hi;
        out.uTau = sqrt(max(out.tauw/rhoW, VSMALL));
        return out;
    }

    scalar tauw = 0.5*(lo + hi);

    for (label iter = 0; iter < coeffs.maxIter; ++iter)
    {
        tauw = 0.5*(lo + hi);
        const scalar fMid = shootResidual(tauw);

        if (fMid > 0)
        {
            hi = tauw;
            fHi = fMid;
        }
        else
        {
            lo = tauw;
            fLo = fMid;
        }

        out.residual = mag(fMid)/positive(in.UM);
        out.iterations = iter + 1;

        if
        (
            mag(hi - lo)/positive(tauw) < coeffs.tolerance
         || out.residual < coeffs.tolerance
        )
        {
            out.converged = 1;
            break;
        }
    }

    if (!out.converged)
    {
        return out;
    }

    out.tauw = tauw;
    out.uTau = sqrt(max(tauw/rhoW, VSMALL));
    out.kappa = kappa;
    out.APlus = APlus;
    out.yPlus = rhoW*out.uTau*in.yM/muW;

    scalar TM =
        temperatureFromU(in.UM, in.UM, in.TM, Ue, Tw, out.Tr, s, Pr);
    TM = min(max(TM, TLo), THi);
    const scalar rhoAtM = rhoTM/TM;
    out.yStar = in.yM*sqrt(rhoAtM*tauw)/positive(coeffs.mu(TM));

    // Constant-stress layer: tau_R(yM) = tau_w - mu dU/dy|yM, and the
    // compressible eddy viscosity (PNAS Eq. 6) is tau_R/(dU/dy).
    if (dUdyM > VSMALL)
    {
        out.mut1 = max(scalar(0), (tauw - muM*dUdyM)/dUdyM);
    }

    // Paper Sec. 2.2: q_w = s tau_w Cp (T_w - T_r)/U_e; positive for heat
    // entering the fluid (Griffin B_q sign = OpenFOAM outward-normal
    // boundary-flux sign).
    out.qw = in.adiabatic ? scalar(0) : s*tauw*Cp*(Tw - out.Tr)/Ue;

    return out;
}


bool isAdiabatic
(
    const Coeffs& coeffs,
    const fvPatchScalarField& Tpatch
)
{
    if (coeffs.thermalMode == "adiabatic")
    {
        return true;
    }

    if (coeffs.thermalMode == "isothermal")
    {
        return false;
    }

    if (Tpatch.fixesValue())
    {
        return false;
    }

    // The Neumann heat-flux boundary condition does not fix a value but is
    // isothermal in intent unless it was itself set to adiabatic.
    if (isA<GriffinWallHeatFluxTemperatureFvPatchScalarField>(Tpatch))
    {
        return refCast<const GriffinWallHeatFluxTemperatureFvPatchScalarField>
        (
            Tpatch
        ).adiabaticIntent();
    }

    return true;
}


PatchData calcPatch
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const Coeffs& coeffs
)
{
    const label patchi = patch.index();
    const label nFaces = patch.size();
    const labelUList& faceCells = patch.faceCells();
    const fvMesh& mesh = turbModel.U().mesh();

    PatchData data(nFaces);

    const fluidThermo& thermo = turbModel.transport();
    const volVectorField& U = turbModel.U();
    const fvPatchVectorField& Uw = U.boundaryField()[patchi];
    const vectorField gradUw(Uw.snGrad());
    const scalarField& y = turbModel.y()[patchi];
    const vectorField nOut(patch.nf());

    const tmp<volScalarField> tmu = turbModel.mu();
    const volScalarField& mu = tmu();
    const tmp<scalarField> tmuPatch = turbModel.mu(patchi);
    const scalarField& muPatch = tmuPatch();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const volScalarField& rho = turbModel.rho();
    const fvPatchScalarField& rhoPatch = rho.boundaryField()[patchi];
    const volScalarField& T = thermo.T();
    const volScalarField& p = thermo.p();

    const fvPatchScalarField& Tpatch = T.boundaryField()[patchi];
    const bool adiabatic = isAdiabatic(coeffs, Tpatch);
    const bool neumannTBC =
        isA<GriffinWallHeatFluxTemperatureFvPatchScalarField>(Tpatch);

    // Physical wall temperature. Under the Neumann coupling the patch T
    // value drifts along the prescribed tangent and cannot supply Tw.
    scalarField TwField(Tpatch);
    if (!adiabatic)
    {
        if (coeffs.Tw > 0)
        {
            TwField = coeffs.Tw;
        }
        else if (neumannTBC)
        {
            FatalErrorInFunction
                << "Patch " << patch.name() << " uses "
                << "GriffinWallHeatFluxTemperature in isothermal mode "
                << "but no positive Tw was given in the wall-function "
                << "coefficients. The patch T value drifts under this "
                << "coupling and cannot supply the wall temperature."
                << nl << exit(FatalError);
        }
    }

    const scalarField pWall(p.boundaryField()[patchi]);
    const scalarField TInternal(Tpatch.patchInternalField());

    // Calorically-perfect-gas form of the model: Cp evaluated at the
    // matching (first-cell) state.
    const tmp<scalarField> tcp = thermo.Cp(pWall, TInternal, patchi);
    const scalarField& cp = tcp();

    // Boundary-layer edge quantities (shared, cached per patch).
    const GriffinBoundaryLayerEdge::PatchEdgeData& edge =
        GriffinBoundaryLayerEdge::New(mesh).edgeData(turbModel, patch, coeffs);

    data.Ue = edge.Ue;
    data.Te = edge.Te;
    data.delta = edge.delta;
    data.edgeStatus = edge.status;

    label invalidCount = 0;
    label outerMatchCount = 0;

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        const vector normal = nOut[facei];
        const vector wallVelocity = Uw[facei];
        const scalar U1 =
            mag(tangential(U[celli] - wallVelocity, normal));
        const scalar T1 = positive(T[celli]);
        const scalar rho1 = positive(rho[celli]);
        const scalar mu1 = positive(mu[celli]);
        const scalar muStress = positive(muPatch[facei]);
        const scalar rhoStress = positive(rhoPatch[facei]);
        const scalar yFace = positive(y[facei]);

        data.T1[facei] = T1;
        data.U1[facei] = U1;
        data.Tw[facei] = TwField[facei];

        if (U1 <= SMALL)
        {
            data.valid[facei] = 1;
            data.converged[facei] = 1;
            data.status[facei] = 1;
            continue;
        }

        const scalar UeFace = data.Ue[facei];
        const scalar TeFace = data.Te[facei];
        const scalar deltaFace = data.delta[facei];

        if
        (
            UeFace <= U1*(1 + 1e-6)
         || TeFace <= SMALL
         || deltaFace <= yFace
        )
        {
            data.status[facei] = 2;
            ++invalidCount;
            continue;
        }

        if (yFace > 0.3*deltaFace)
        {
            ++outerMatchCount;
        }

        // Initial tau_w from the current effective stress at the patch.
        const vector gradU = tangential(gradUw[facei], normal);
        const scalar magGradU = mag(gradU);
        scalar tauwGuess = rhoStress*(nutw[facei] + nuw[facei])*magGradU;
        if (tauwGuess <= SMALL)
        {
            tauwGuess = muStress*U1/yFace;
        }

        InnerInput in;
        in.yM = yFace;
        in.UM = U1;
        in.TM = T1;
        in.rhoM = rho1;
        in.Tw = TwField[facei];
        in.adiabatic = adiabatic;
        in.Ue = UeFace;
        in.Te = TeFace;
        in.delta = deltaFace;
        in.Cp = cp[facei];
        in.tauwGuess = tauwGuess;

        const InnerResult res = solveInner(in, coeffs);

        data.residual[facei] = res.residual;
        data.iterations[facei] = res.iterations;
        data.Tr[facei] = res.Tr;
        data.Tw[facei] = res.Tw;

        if
        (
            !res.converged
         || res.uTau <= SMALL
         || !std::isfinite(res.tauw)
         || !std::isfinite(res.qw)
        )
        {
            data.status[facei] = 3;
            ++invalidCount;
            continue;
        }

        // Invert OpenFOAM's discrete stress relation with the patch-field
        // viscosity and density that the stress operator actually uses;
        // under the Neumann coupling these differ from the physical wall
        // values because the patch T drifts.
        scalar nutShear = 0;
        if (magGradU > ROOTVSMALL)
        {
            nutShear =
                max
                (
                    scalar(0),
                    (res.tauw - muStress*magGradU)/(rhoStress*magGradU)
                );
        }
        else
        {
            nutShear = res.mut1/rhoStress;
        }

        data.uTau[facei] = res.uTau;
        data.tauw[facei] = res.tauw;
        data.yPlus[facei] = res.yPlus;
        data.yStar[facei] = res.yStar;
        data.uPlus[facei] = U1/positive(res.uTau);
        data.ReTauStar[facei] = res.ReTauStar;
        data.kappaCal[facei] = res.kappa;
        data.APlusCal[facei] = res.APlus;
        data.qw[facei] = res.qw;
        data.nutShear[facei] = nutShear;
        data.mut1[facei] = res.mut1;
        data.nuTilda1[facei] = solveNuTilda
        (
            res.mut1/rho1,
            mu1/rho1,
            coeffs.cv1,
            coeffs.maxIter,
            coeffs.tolerance
        );
        data.converged[facei] = res.converged;
        data.valid[facei] = 1;
        data.status[facei] = 0;
    }

    if (invalidCount)
    {
        if (coeffs.strict)
        {
            FatalErrorInFunction
                << invalidCount << " face(s) on patch " << patch.name()
                << " could not evaluate the Griffin wall function. "
                << "Inspect status diagnostics; strict=false enables the "
                << "laminar fallback." << nl
                << exit(FatalError);
        }

        static bool warned = false;
        if (!warned)
        {
            WarningInFunction
                << "Griffin wall function used its laminar fallback on "
                << invalidCount << " face(s) of patch " << patch.name()
                << ". Enable diagnostics or strict=true to locate invalid "
                << "edge/matching states." << endl;
            warned = true;
        }
    }

    if (outerMatchCount)
    {
        static bool warnedOuter = false;
        if (!warnedOuter)
        {
            WarningInFunction
                << "Griffin wall function: the matching point (first-cell "
                << "centre) lies beyond 0.3 delta on " << outerMatchCount
                << " face(s) of patch " << patch.name()
                << ". The model assumes the matching data lies in the "
                << "inner layer (paper Sec. 5); refine the wall-normal "
                << "mesh." << endl;
            warnedOuter = true;
        }
    }

    return data;
}

} // End namespace GriffinWallFunction
} // End namespace Foam

// ************************************************************************* //
