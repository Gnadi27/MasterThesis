#include "WangNeumannWallFunction.H"
#include "WangNeumannWallHeatFluxTemperatureFvPatchScalarField.H"

#include "dictionary.H"
#include "fvPatchField.H"
#include "fvMesh.H"
#include "Ostream.H"
#include "surfaceFields.H"
#include "volFields.H"

#include <cmath>
#include <functional>

namespace Foam
{
namespace
{

inline scalar positive(const scalar x)
{
    return max(x, VSMALL);
}


//- Guarded natural logarithm
inline scalar lnPos(const scalar x)
{
    return std::log(max(x, SMALL));
}


//- Composite Simpson integration of f over [a, b] with n (even) intervals
scalar simpson
(
    const scalar a,
    const scalar b,
    const label n,
    const std::function<scalar(scalar)>& f
)
{
    if (b <= a)
    {
        return 0;
    }

    const label nn = max(label(2), n + (n % 2));
    const scalar h = (b - a)/nn;

    scalar sum = f(a) + f(b);
    for (label i = 1; i < nn; ++i)
    {
        sum += (i % 2 ? 4.0 : 2.0)*f(a + i*h);
    }

    return sum*h/3.0;
}


//- All per-face inputs of the analytical solution
struct FaceInput
{
    scalar kP;          // nodal turbulence energy
    scalar yP;          // wall distance of the node
    scalar yn;          // wall-normal extent of the near-wall cell
    scalar Un;          // wall-parallel velocity at the north face
    scalar Tn;          // temperature at the north face
    scalar conv;        // nodal momentum convection [kg/(m^2 s^2)]
    scalar dpdx;        // nodal wall-parallel pressure gradient [Pa/m]
    scalar divE;        // nodal div((rho E + p) U) [W/m^3]
    scalar muv;         // molecular viscosity at Tv
    scalar rhov;        // density at (p_w, Tv)
    scalar muw;         // molecular viscosity at Tw
    scalar cp;          // cp at (p_w, Tw)
    scalar Pr;          // molecular Prandtl number
    scalar Tw;          // current wall temperature
    bool adiabatic;
};


//- All per-face outputs
struct FaceOutput
{
    scalar A1, tauw, qw, Ath1;
    scalar TvNew, TwNew;
    scalar Pk, epsBar, epsP, dissAvg;
    scalar yStar, yStarN, bMu, D, C, Dth;
    bool thin;
};


//- Solve the Wang AWF^Cmp on one face (paper Eqs. WF30-WF45; thesis
//  Appendix E when y*n < y*v). All formulas are assembled from the two
//  first-integral identities
//      mu       dU1/dy* = D y*^3/3 + C y* + A1        (sublayer)
//      mu_v Y   dU2/dy* = D y*^3/3 + C y* + A2        (turbulent layer)
//  and the profile/continuity conditions, which reproduces the printed
//  monolithic N (WF35), B2 (WF34), BT (WF41), Ath1 (WF44) and Tw (WF42)
//  expressions exactly with D1=D2=D, C1=C2=C, Dth1=Dth2=Dth.
FaceOutput solveFace
(
    const FaceInput& in,
    const WangNeumannWallFunction::Coeffs& coeffs
)
{
    FaceOutput out{};

    // Floor at SMALL (not VSMALL): kP appears in negative powers and an
    // uninitialised k = 0 would overflow the scaled coefficients.
    const scalar kP = max(in.kP, SMALL);
    const scalar sqrtk = std::sqrt(kP);
    const scalar muv = positive(in.muv);
    const scalar rhov = positive(in.rhov);
    const scalar muw = positive(in.muw);
    const scalar cp = positive(in.cp);
    const scalar Pr = positive(in.Pr);
    const scalar yn = positive(in.yn);

    // y* = y rho_v sqrt(kP)/mu_v (paper Eq. WF26)
    const scalar scale = rhov*sqrtk/muv;
    const scalar ysP = scale*positive(in.yP);
    const scalar ysN = max(scale*yn, SMALL);

    const scalar alpha = coeffs.alpha();
    const scalar alphat = max(alpha*Pr/positive(coeffs.Prt), SMALL);

    const bool thin = (ysN <= coeffs.yStarV);
    const scalar yv = thin ? ysN : coeffs.yStarV;

    // Scaled convection and pressure-gradient coefficients (WF18/WF23)
    const scalar denomK = positive(sqr(rhov)*kP);
    const scalar D = sqr(muv)/(denomK*sqr(ysN))*in.conv;
    const scalar C = sqr(muv)/denomK*in.dpdx;
    const scalar Dth = sqr(muv)/(cp*denomK*sqr(ysN))*in.divE;

    // Hyperbolic molecular-viscosity variation (WF28; thesis E.1 in the
    // thin-cell branch, where the "v" state coincides with the "n" state)
    scalar bMu = 0;
    if (coeffs.viscosityVariation == "hyperbolic")
    {
        bMu = (muw - muv)/(muw*yv);
    }

    // ---------------- Momentum (WF30-WF36 / E.1-E.2) ----------------

    // Gm(s) of the outer-layer double integral (bracket of WF31)
    const auto lnY = [&](const scalar s) -> scalar
    {
        return lnPos(1 + alpha*(s - yv));
    };
    const scalar oneMinusAYv = 1 - alpha*yv;
    const auto Gm = [&](const scalar s) -> scalar
    {
        return
            (sqr(alpha)/3)*pow3(s)
          - (alpha*oneMinusAYv/2)*sqr(s)
          + sqr(oneMinusAYv)*s
          - (pow3(oneMinusAYv)/alpha)*lnY(s);
    };

    scalar A1 = 0;
    scalar A2 = 0;
    scalar B2 = 0;
    scalar lnYn = 0;

    if (!thin)
    {
        lnYn = lnY(ysN);

        // Sublayer polynomial part of U1 at y*v (continuity terms)
        const scalar sub =
            D*pow4(yv)/12 + C*sqr(yv)/2
          - bMu*(D*pow5(yv)/60 + C*pow3(yv)/6);

        const scalar outerV = (D/(3*pow3(alpha)))*Gm(yv) + (C/alpha)*yv;

        // N of WF35 (with D1=D2, C1=C2)
        const scalar N =
            (D/(3*pow3(alpha)))*Gm(ysN)
          + (C/alpha)*(ysN - (oneMinusAYv/alpha)*lnYn)
          + sub
          - outerV;

        scalar denomA1 = lnYn/alpha + yv - 0.5*bMu*sqr(yv);
        denomA1 = (denomA1 >= 0 ? positive(denomA1) : min(denomA1, -VSMALL));

        A1 = (muv*in.Un - N)/denomA1;                       // WF32
        A2 = A1;                                            // WF33, D1=D2
        B2 = sub + A1*(yv - 0.5*bMu*sqr(yv)) - outerV;      // WF34
    }
    else
    {
        // Thesis E.2 (whole cell inside the sublayer; mu_n == mu_v here)
        const scalar poly =
            D*pow4(ysN)/12 + C*sqr(ysN)/2
          - bMu*(D*pow5(ysN)/60 + C*pow3(ysN)/6);

        scalar denomA1 = ysN - 0.5*bMu*sqr(ysN);
        denomA1 = (denomA1 >= 0 ? positive(denomA1) : min(denomA1, -VSMALL));

        A1 = (muv*in.Un - poly)/denomA1;
        A2 = A1;
        B2 = 0;
    }

    // The papers print tau_wall = -(rho_v sqrt(kP)/mu_v) A1 (T8.12), the
    // reaction force on the fluid (negative for attached flow, cf. Craft
    // 2002 Fig. 6). OpenFOAM's wall-stress path reproduces the
    // conventional mu dU/dy|_0 = +(rho_v sqrt(kP)/mu_v) A1, positive for
    // attached flow, so that sign is stored here.
    const scalar tauw = scale*A1;

    // Analytical profiles (WF30/WF31); thin branch: yv = y*n
    const auto U1 = [&](const scalar s) -> scalar
    {
        return
            (
                D*pow4(s)/12 + C*sqr(s)/2 + A1*s
              + bMu*
                (
                    D*pow5(s)/15 - D*yv*pow4(s)/12 + C*pow3(s)/3
                  + 0.5*(A1 - C*yv)*sqr(s) - A1*yv*s
                )
            )/muv;
    };
    const auto U2 = [&](const scalar s) -> scalar
    {
        return
            (
                (D/(3*pow3(alpha)))*Gm(s)
              + (C/alpha)*(s - (oneMinusAYv/alpha)*lnY(s))
              + (A2/alpha)*lnY(s)
              + B2
            )/muv;
    };

    // First integrals (exact): mu dU1/dy* and mu_v Y dU2/dy*
    const auto tau1 = [&](const scalar s) -> scalar
    {
        return D*pow3(s)/3 + C*s + A1;
    };
    const auto tau2 = [&](const scalar s) -> scalar
    {
        return D*pow3(s)/3 + C*s + A2;
    };

    // ------------- k-equation quantities (WF13, WF16, Craft 4.9) -------------

    scalar Pk = 0;
    if (!thin)
    {
        // integrand = mu_t (dU2/dy*)^2 = alpha(s-yv) tau2^2/(mu_v Y^2)
        Pk =
            (scale/yn)
           *simpson
            (
                yv,
                ysN,
                coeffs.nIntegrationPoints,
                [&](const scalar s) -> scalar
                {
                    const scalar Y = max(1 + alpha*(s - yv), SMALL);
                    return alpha*(s - yv)*sqr(tau2(s))/(muv*sqr(Y));
                }
            );
    }

    const scalar k32 = kP*sqrtk;
    const scalar yd = coeffs.yStarD();
    scalar epsBar = 0;
    if (ysN > yd)
    {
        epsBar = (2*k32/yd + (k32/coeffs.cl)*lnPos(ysN/yd))/yn;    // WF16
    }
    else
    {
        const scalar ydPhys = positive(yd/scale);
        epsBar = 2*(muv/rhov)*kP/sqr(ydPhys);                      // Craft 4.8
    }

    const scalar epsP = k32/(coeffs.cl*positive(in.yP));           // Craft 4.9

    // ---------------- Thermal (WF37-WF45 / E.3-E.4) ----------------

    const scalar oneMinusAtYv = 1 - alphat*yv;
    const auto lnYT = [&](const scalar s) -> scalar
    {
        return lnPos(1 + alphat*(s - yv));
    };
    const auto Fth = [&](const scalar s) -> scalar
    {
        return
            (Dth/(3*pow3(alphat)))
           *(
                (sqr(alphat)/3)*pow3(s)
              - (alphat*oneMinusAtYv/2)*sqr(s)
              + sqr(oneMinusAtYv)*s
              - (pow3(oneMinusAtYv)/alphat)*lnYT(s)
            );
    };

    // Viscous-dissipation integrals of the analytical temperature
    // (paper: "evaluated numerically"); the first-integral identities
    // reduce the integrands to smooth closed-form functions.
    const auto I1 = [&](const scalar b) -> scalar
    {
        return
            simpson
            (
                0,
                b,
                coeffs.nIntegrationPoints,
                [&](const scalar s) -> scalar
                {
                    return (1 + bMu*(s - yv))*U1(s)*tau1(s);
                }
            )/cp;
    };
    const auto I2 = [&](const scalar b) -> scalar
    {
        return
            simpson
            (
                yv,
                b,
                coeffs.nIntegrationPoints,
                [&](const scalar s) -> scalar
                {
                    const scalar YT = max(1 + alphat*(s - yv), SMALL);
                    return U2(s)*tau2(s)/YT;
                }
            )/cp;
    };

    scalar Ath1 = 0;
    scalar qw = 0;
    scalar TwNew = in.Tw;
    scalar TvNew = in.Tn;

    if (!thin)
    {
        const scalar lnYTn = lnYT(ysN);
        const scalar I1v = I1(yv);
        const scalar I2n = I2(ysN);

        // BT (WF41)
        const scalar BT =
            Dth*pow4(yv)/12 - (Dth/60)*bMu*pow5(yv) - I1v - Fth(yv);

        scalar denomTh = lnYTn/alphat + yv - 0.5*bMu*sqr(yv);
        denomTh = (denomTh >= 0 ? positive(denomTh) : min(denomTh, -VSMALL));

        if (in.adiabatic)
        {
            Ath1 = 0;
            qw = 0;
            // Tw from WF42 with A_th1 = 0
            TwNew = in.Tn - (Pr/muv)*(Fth(ysN) + BT - I2n);
        }
        else
        {
            // Ath1 from WF44
            Ath1 =
                ((muv/Pr)*(in.Tn - in.Tw) - Fth(ysN) - BT + I2n)/denomTh;
            qw = -(rhov*cp*sqrtk/muv)*Ath1;                        // WF43
        }

        // Tv = T1(y*v) (continuity value; polynomial evaluated at y*v)
        TvNew =
            (Pr/muv)
           *(
                Dth*pow4(yv)/12 + Ath1*yv
              - bMu*(Dth*pow5(yv)/60 + 0.5*Ath1*sqr(yv))
              - I1v
            )
          + TwNew;
    }
    else
    {
        // Thesis Appendix E (paper-consistent: no Cth terms)
        const scalar I1n = I1(ysN);

        scalar denomTh = ysN - 0.5*bMu*sqr(ysN);
        denomTh = (denomTh >= 0 ? positive(denomTh) : min(denomTh, -VSMALL));

        if (in.adiabatic)
        {
            Ath1 = 0;
            qw = 0;
            // E.3
            TwNew =
                in.Tn
              - (Pr/muv)
               *(
                    Dth*pow4(ysN)/12
                  - bMu*Dth*pow5(ysN)/60
                  - I1n
                );
        }
        else
        {
            // E.4
            Ath1 =
                (
                    (muv/Pr)*(in.Tn - in.Tw)
                  - Dth*pow4(ysN)/12
                  + bMu*Dth*pow5(ysN)/60
                  + I1n
                )/denomTh;
            qw = -(rhov*cp*sqrtk/muv)*Ath1;
        }

        TvNew = in.Tn;
    }

    TwNew = positive(TwNew);
    TvNew = positive(TvNew);

    // -------- CMAWF-D cell-averaged wall-normal dissipation (WF24/T8.64) ----
    // (1/yn) [ (mu+mu_t) U dU/dy ]_0^{yn} with U(0) = 0:
    const scalar tauN = scale*tau2(ysN);
    const scalar dissAvg = in.Un*tauN/yn;

    out.A1 = A1;
    out.tauw = tauw;
    out.qw = qw;
    out.Ath1 = Ath1;
    out.TvNew = TvNew;
    out.TwNew = TwNew;
    out.Pk = Pk;
    out.epsBar = epsBar;
    out.epsP = epsP;
    out.dissAvg = dissAvg;
    out.yStar = ysP;
    out.yStarN = ysN;
    out.bMu = bMu;
    out.D = D;
    out.C = C;
    out.Dth = Dth;
    out.thin = thin;

    return out;
}

} // End anonymous namespace


namespace WangNeumannWallFunction
{

Coeffs::Coeffs()
:
    yStarV(10.8),
    cl(2.55),
    Cmu(0.09),
    Prt(0.9),
    nIntegrationPoints(40),
    TvIterations(3),
    TvRelaxation(0.5),
    viscosityVariation("hyperbolic"),
    thermalMode("auto"),
    Tw(-1),
    wallTransportMode("fromPatch"),
    wallTransportAs(0),
    wallTransportTs(0),
    wallTransportPr(-1),
    wallTransportMu(0),
    wallTransportKappa(0)
{}


Coeffs::Coeffs(const dictionary& dict)
:
    yStarV(dict.getOrDefault<scalar>("yStarV", 10.8)),
    cl(dict.getOrDefault<scalar>("cl", 2.55)),
    Cmu(dict.getOrDefault<scalar>("Cmu", 0.09)),
    Prt(dict.getOrDefault<scalar>("Prt", 0.9)),
    nIntegrationPoints(dict.getOrDefault<label>("nIntegrationPoints", 40)),
    TvIterations(dict.getOrDefault<label>("TvIterations", 3)),
    TvRelaxation(dict.getOrDefault<scalar>("TvRelaxation", 0.5)),
    viscosityVariation
    (
        dict.getOrDefault<word>("viscosityVariation", "hyperbolic")
    ),
    thermalMode(dict.getOrDefault<word>("thermalMode", "auto")),
    Tw(dict.getOrDefault<scalar>("Tw", -1)),
    wallTransportMode("fromPatch"),
    wallTransportAs(0),
    wallTransportTs(0),
    wallTransportPr(-1),
    wallTransportMu(0),
    wallTransportKappa(0)
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

    if
    (
        viscosityVariation != "hyperbolic"
     && viscosityVariation != "constant"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unsupported viscosityVariation " << viscosityVariation << nl
            << "Supported values are hyperbolic (paper AWF^Cmp) and "
            << "constant (b_mu = 0, CMAWF limit)"
            << exit(FatalIOError);
    }

    const dictionary* wt = dict.findDict("wallTransport");
    if (wt)
    {
        wallTransportMode = wt->getOrDefault<word>("mode", "fromPatch");

        if (wallTransportMode == "sutherland")
        {
            wallTransportAs = wt->get<scalar>("As");
            wallTransportTs = wt->get<scalar>("Ts");
            wallTransportPr = wt->getOrDefault<scalar>("Pr", -1);
        }
        else if (wallTransportMode == "constant")
        {
            wallTransportMu = wt->get<scalar>("mu");
            wallTransportKappa = wt->get<scalar>("kappa");
        }
        else if (wallTransportMode != "fromPatch")
        {
            FatalIOErrorInFunction(*wt)
                << "Unsupported wallTransport mode " << wallTransportMode
                << nl
                << "Supported values are fromPatch, sutherland, and constant"
                << exit(FatalIOError);
        }
    }
}


void Coeffs::writeEntries(Ostream& os) const
{
    os.writeEntryIfDifferent<scalar>("yStarV", 10.8, yStarV);
    os.writeEntryIfDifferent<scalar>("cl", 2.55, cl);
    os.writeEntryIfDifferent<scalar>("Cmu", 0.09, Cmu);
    os.writeEntryIfDifferent<scalar>("Prt", 0.9, Prt);
    os.writeEntryIfDifferent<label>
    (
        "nIntegrationPoints", 40, nIntegrationPoints
    );
    os.writeEntryIfDifferent<label>("TvIterations", 3, TvIterations);
    os.writeEntryIfDifferent<scalar>("TvRelaxation", 0.5, TvRelaxation);
    os.writeEntryIfDifferent<word>
    (
        "viscosityVariation", word("hyperbolic"), viscosityVariation
    );
    os.writeEntryIfDifferent<word>("thermalMode", word("auto"), thermalMode);
    os.writeEntryIfDifferent<scalar>("Tw", -1, Tw);

    if (wallTransportMode != "fromPatch")
    {
        os.beginBlock("wallTransport");
        os.writeEntry("mode", wallTransportMode);
        if (wallTransportMode == "sutherland")
        {
            os.writeEntry("As", wallTransportAs);
            os.writeEntry("Ts", wallTransportTs);
            if (wallTransportPr > 0)
            {
                os.writeEntry("Pr", wallTransportPr);
            }
        }
        else if (wallTransportMode == "constant")
        {
            os.writeEntry("mu", wallTransportMu);
            os.writeEntry("kappa", wallTransportKappa);
        }
        os.endBlock();
    }
}


PatchData::PatchData(const label size)
:
    tauw(size, Zero),
    A1(size, Zero),
    nutShear(size, Zero),
    yStar(size, Zero),
    yStarN(size, Zero),
    uTau(size, Zero),
    qw(size, Zero),
    Tw(size, Zero),
    Tv(size, Zero),
    Ath1(size, Zero),
    Pk(size, Zero),
    epsBar(size, Zero),
    epsP(size, Zero),
    dissAvg(size, Zero),
    dissLocal(size, Zero),
    muw(size, Zero),
    muv(size, Zero),
    rhov(size, Zero),
    cpw(size, Zero),
    Pr(size, Zero),
    bMu(size, Zero),
    Dcoef(size, Zero),
    Ccoef(size, Zero),
    DthCoef(size, Zero),
    Un(size, Zero),
    Tn(size, Zero),
    yn(size, Zero),
    thinCell(size, Zero)
{}


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

    if (isA<WangNeumannWallHeatFluxTemperatureFvPatchScalarField>(Tpatch))
    {
        return refCast
        <
            const WangNeumannWallHeatFluxTemperatureFvPatchScalarField
        >(Tpatch).adiabaticIntent();
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

    PatchData data(nFaces);

    const fvMesh& mesh = patch.boundaryMesh().mesh();
    const fluidThermo& thermo = turbModel.transport();

    const volVectorField& U = turbModel.U();
    const fvPatchVectorField& Uw = U.boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));

    const volScalarField& TField = thermo.T();
    const volScalarField& pField = thermo.p();
    const volScalarField& heField = thermo.he();
    const volScalarField& rho = turbModel.rho();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tmuPatch = turbModel.mu(patchi);
    const scalarField& muPatch = tmuPatch();
    const fvPatchScalarField& rhoPatch = rho.boundaryField()[patchi];

    const fvPatchScalarField& Tpatch = TField.boundaryField()[patchi];
    const scalarField Tcell(Tpatch.patchInternalField());
    const scalarField pWall(pField.boundaryField()[patchi]);

    const bool adiabatic = isAdiabatic(coeffs, Tpatch);
    const bool neumannTBC =
        isA<WangNeumannWallHeatFluxTemperatureFvPatchScalarField>(Tpatch);

    if (!adiabatic && coeffs.Tw <= 0)
    {
        if (neumannTBC)
        {
            FatalErrorInFunction
                << "Patch " << patch.name() << " uses "
                << "WangNeumannWallHeatFluxTemperature in isothermal mode "
                << "but no positive Tw was given in the wall-function "
                << "coefficients. The patch T value drifts under this "
                << "coupling and cannot supply the wall temperature."
                << nl << exit(FatalError);
        }
    }

    if (!adiabatic && coeffs.wallTransportMode == "fromPatch" && neumannTBC)
    {
        FatalErrorInFunction
            << "Patch " << patch.name() << " uses "
            << "WangNeumannWallHeatFluxTemperature in isothermal mode with "
            << "wallTransport mode fromPatch. The patch transport "
            << "properties are evaluated at the drifted patch T and are "
            << "not the wall values; select wallTransport sutherland or "
            << "constant." << nl << exit(FatalError);
    }

    // ----------------- transport-property helpers -----------------

    const auto muAt = [&](const scalarField& T) -> tmp<scalarField>
    {
        if (coeffs.wallTransportMode == "sutherland")
        {
            auto tres = tmp<scalarField>::New(T.size(), Zero);
            scalarField& res = tres.ref();
            forAll(res, i)
            {
                const scalar Ti = max(T[i], SMALL);
                res[i] =
                    coeffs.wallTransportAs*Ti*std::sqrt(Ti)
                   /(Ti + coeffs.wallTransportTs);
            }
            return tres;
        }
        else if (coeffs.wallTransportMode == "constant")
        {
            return tmp<scalarField>::New(T.size(), coeffs.wallTransportMu);
        }

        // fromPatch (adiabatic mode only): patch molecular viscosity
        return tmp<scalarField>::New(muPatch);
    };

    // Molecular Prandtl number: wallTransport Pr if given, else from the
    // wallTransport conductivity (constant mode) or modified Eucken
    // (sutherland without Pr), else from the patch fields (fromPatch).
    const auto prandtl = [&]
    (
        const scalarField& T,
        const scalarField& mu,
        const scalarField& cp
    ) -> tmp<scalarField>
    {
        auto tres = tmp<scalarField>::New(nFaces, Zero);
        scalarField& res = tres.ref();

        if (coeffs.wallTransportMode == "sutherland")
        {
            if (coeffs.wallTransportPr > 0)
            {
                res = coeffs.wallTransportPr;
            }
            else
            {
                const tmp<scalarField> tcv = thermo.Cv(pWall, T, patchi);
                const scalarField& cv = tcv();
                forAll(res, i)
                {
                    const scalar cvi = positive(cv[i]);
                    const scalar R = max(cp[i] - cvi, scalar(0));
                    const scalar kappaTh = mu[i]*cvi*(1.32 + 1.77*R/cvi);
                    res[i] = positive(cp[i]*mu[i]/positive(kappaTh));
                }
            }
        }
        else if (coeffs.wallTransportMode == "constant")
        {
            forAll(res, i)
            {
                res[i] =
                    positive
                    (
                        cp[i]*coeffs.wallTransportMu
                       /positive(coeffs.wallTransportKappa)
                    );
            }
        }
        else
        {
            const tmp<scalarField> tkappa = turbModel.kappa(patchi);
            const scalarField& kappaTh = tkappa();
            forAll(res, i)
            {
                res[i] = positive(cp[i]*mu[i]/positive(kappaTh[i]));
            }
        }

        return tres;
    };

    // ----------------- geometry / face-interpolation helpers -----------------

    const labelUList& own = mesh.faceOwner();
    const labelUList& nei = mesh.faceNeighbour();
    const surfaceScalarField& weights = mesh.surfaceInterpolation::weights();
    // Primitive-mesh geometry: valid for every face, including faces of
    // empty/wedge patches whose fv boundary fields have zero size.
    const vectorField& faceAreas = mesh.faceAreas();
    const vectorField& faceCentres = mesh.faceCentres();
    const scalarField& V = mesh.V();

    // Boundary patch and local index of a boundary face; local index is
    // set to -1 when the fv boundary field does not cover the face
    // (empty patches), in which case callers fall back to the owner-cell
    // value -- exact for empty patches, whose paired opposite faces
    // cancel in the Gauss sums.
    const auto facePatchAndIndex = [&](const label f, label& li) -> label
    {
        const label pI = mesh.boundaryMesh().whichPatch(f);
        li = f - mesh.boundaryMesh()[pI].start();
        return pI;
    };

    const auto faceSf = [&](const label f) -> vector
    {
        return faceAreas[f];
    };

    const auto faceCf = [&](const label f) -> vector
    {
        return faceCentres[f];
    };

    const auto faceScalar =
        [&](const label f, const volScalarField& fld) -> scalar
    {
        if (mesh.isInternalFace(f))
        {
            const scalar w = weights[f];
            return w*fld[own[f]] + (1 - w)*fld[nei[f]];
        }
        label li;
        const label pI = facePatchAndIndex(f, li);
        if (li < fld.boundaryField()[pI].size())
        {
            return fld.boundaryField()[pI][li];
        }
        return fld[own[f]];
    };

    const auto faceVector =
        [&](const label f, const volVectorField& fld) -> vector
    {
        if (mesh.isInternalFace(f))
        {
            const scalar w = weights[f];
            return w*fld[own[f]] + (1 - w)*fld[nei[f]];
        }
        label li;
        const label pI = facePatchAndIndex(f, li);
        if (li < fld.boundaryField()[pI].size())
        {
            return fld.boundaryField()[pI][li];
        }
        return fld[own[f]];
    };

    // Gauss gradient of a scalar field at one cell
    const auto cellGradScalar =
        [&](const label celli, const volScalarField& fld) -> vector
    {
        vector g(Zero);
        for (const label f : mesh.cells()[celli])
        {
            const vector sf = faceSf(f);
            const scalar sgn =
                (!mesh.isInternalFace(f) || own[f] == celli) ? 1 : -1;
            g += sgn*sf*faceScalar(f, fld);
        }
        return g/V[celli];
    };

    // Gauss gradient of a vector field at one cell; (grad)_ij = d u_j/d x_i
    const auto cellGradVector =
        [&](const label celli, const volVectorField& fld) -> tensor
    {
        tensor g(Zero);
        for (const label f : mesh.cells()[celli])
        {
            const vector sf = faceSf(f);
            const scalar sgn =
                (!mesh.isInternalFace(f) || own[f] == celli) ? 1 : -1;
            g += sgn*(sf*faceVector(f, fld));
        }
        return g/V[celli];
    };

    // Internal energy per cell/boundary face; total energy E = e + |u|^2/2.
    // If the thermo energy variable is an enthalpy, e = he - p/rho.
    const bool heIsEnthalpy =
        (heField.member().size() && heField.member()[0] == 'h');

    const auto cellE = [&](const label celli) -> scalar
    {
        const scalar e =
            heIsEnthalpy
          ? heField[celli] - pField[celli]/positive(rho[celli])
          : heField[celli];
        return e + 0.5*magSqr(U[celli]);
    };

    const auto boundaryE = [&](const label pI, const label li) -> scalar
    {
        const scalar heB = heField.boundaryField()[pI][li];
        const scalar pB = pField.boundaryField()[pI][li];
        const scalar rhoB = positive(rho.boundaryField()[pI][li]);
        const scalar e = heIsEnthalpy ? heB - pB/rhoB : heB;
        return e + 0.5*magSqr(U.boundaryField()[pI][li]);
    };

    // Gauss divergence of (rho E + p) u at one cell [W/m^3]
    const auto cellDivPhiE = [&](const label celli) -> scalar
    {
        scalar div = 0;
        for (const label f : mesh.cells()[celli])
        {
            const vector sf = faceSf(f);
            const scalar sgn =
                (!mesh.isInternalFace(f) || own[f] == celli) ? 1 : -1;

            scalar rhoEpP;
            vector uf;
            if (mesh.isInternalFace(f))
            {
                const scalar w = weights[f];
                const scalar Ef = w*cellE(own[f]) + (1 - w)*cellE(nei[f]);
                const scalar rhof = w*rho[own[f]] + (1 - w)*rho[nei[f]];
                const scalar pf = w*pField[own[f]] + (1 - w)*pField[nei[f]];
                uf = w*U[own[f]] + (1 - w)*U[nei[f]];
                rhoEpP = rhof*Ef + pf;
            }
            else
            {
                label li;
                const label pI = facePatchAndIndex(f, li);
                if (li < rho.boundaryField()[pI].size())
                {
                    const scalar Ef = boundaryE(pI, li);
                    const scalar rhof = rho.boundaryField()[pI][li];
                    const scalar pf = pField.boundaryField()[pI][li];
                    uf = U.boundaryField()[pI][li];
                    rhoEpP = rhof*Ef + pf;
                }
                else
                {
                    // Empty-patch face: owner-cell values (paired faces
                    // cancel in the Gauss sum)
                    const scalar Ef = cellE(own[f]);
                    uf = U[own[f]];
                    rhoEpP = rho[own[f]]*Ef + pField[own[f]];
                }
            }

            div += sgn*(sf & (rhoEpP*uf));
        }
        return div/V[celli];
    };

    // ----------------- per-face geometric/nodal precomputation -----------------

    const vectorField nf(patch.nf());          // outward patch normal
    const vectorField& CfWall = patch.Cf();

    scalarField Un(nFaces, Zero);
    scalarField Tn(nFaces, Zero);
    scalarField ynField(nFaces, Zero);
    scalarField conv(nFaces, Zero);
    scalarField dpdx(nFaces, Zero);
    scalarField divE(nFaces, Zero);
    scalarField dUtdn(nFaces, Zero);

    const tmp<volScalarField> tmuEff = turbModel.muEff();
    const volScalarField& muEff = tmuEff();

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        const label wallGlobal = patch.start() + facei;
        const vector nHat = -nf[facei];        // into the fluid

        // Wall-parallel unit vector along the near-wall cell velocity
        const vector urel = U[celli] - Uw[facei];
        vector tHat = urel - (urel & nHat)*nHat;
        const scalar magT = mag(tHat);
        if (magT > ROOTVSMALL)
        {
            tHat /= magT;
        }
        else
        {
            tHat = Zero;
        }

        // North face: the cell face most anti-parallel to the wall face
        label bestFace = -1;
        scalar bestAlign = 0.5;
        for (const label f : mesh.cells()[celli])
        {
            if (f == wallGlobal)
            {
                continue;
            }
            vector sf = faceSf(f);
            if (mesh.isInternalFace(f) && own[f] != celli)
            {
                sf = -sf;
            }
            const scalar align = (sf/positive(mag(sf))) & nHat;
            if (align > bestAlign)
            {
                bestAlign = align;
                bestFace = f;
            }
        }

        scalar yn = 2*y[facei];
        scalar UnFace = urel & tHat;
        scalar TnFace = TField[celli];
        if (bestFace >= 0)
        {
            const scalar ynCand = (faceCf(bestFace) - CfWall[facei]) & nHat;
            if (ynCand > SMALL)
            {
                yn = ynCand;
                UnFace = (faceVector(bestFace, U) - Uw[facei]) & tHat;
                TnFace = faceScalar(bestFace, TField);
            }
        }

        Un[facei] = UnFace;
        Tn[facei] = positive(TnFace);
        ynField[facei] = yn;

        // Nodal convection/pressure/energy-convection terms
        if (magT > ROOTVSMALL)
        {
            const tensor gU = cellGradVector(celli, U);
            const vector gradUt = gU & tHat;      // grad of (u . tHat)
            const scalar Ut = urel & tHat;
            const scalar UnComp = urel & nHat;

            conv[facei] =
                rho[celli]*(Ut*(tHat & gradUt) + UnComp*(nHat & gradUt));
            dpdx[facei] = tHat & cellGradScalar(celli, pField);
            dUtdn[facei] = nHat & gradUt;
        }

        divE[facei] = cellDivPhiE(celli);
    }

    // ----------------- wall/sublayer state fixed-point loop -----------------

    scalarField Tw(nFaces, Zero);
    forAll(Tw, facei)
    {
        Tw[facei] =
            (!adiabatic && coeffs.Tw > 0)
          ? coeffs.Tw
          : positive(Tcell[facei]);
    }

    scalarField Tv(Tcell);
    forAll(Tv, facei)
    {
        Tv[facei] = positive(Tv[facei]);
    }

    const labelList cells(faceCells);

    List<FaceOutput> sol(nFaces);

    const label nIter = max(label(1), coeffs.TvIterations);
    for (label iter = 0; iter < nIter; ++iter)
    {
        const tmp<scalarField> tmuv = muAt(Tv);
        const scalarField& muv = tmuv();
        const tmp<scalarField> tmuw = muAt(Tw);
        const scalarField& muw = tmuw();
        const tmp<scalarField> trhov = thermo.rhoEoS(pWall, Tv, cells);
        const scalarField& rhov = trhov();
        const tmp<scalarField> tcpw = thermo.Cp(pWall, Tw, patchi);
        const scalarField& cpw = tcpw();
        const tmp<scalarField> tPr = prandtl(Tv, muv, cpw);
        const scalarField& Pr = tPr();

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];

            FaceInput in;
            in.kP = k[celli];
            in.yP = y[facei];
            in.yn = ynField[facei];
            in.Un = Un[facei];
            in.Tn = Tn[facei];
            in.conv = conv[facei];
            in.dpdx = dpdx[facei];
            in.divE = divE[facei];
            in.muv = muv[facei];
            in.rhov = rhov[facei];
            in.muw = muw[facei];
            in.cp = cpw[facei];
            in.Pr = Pr[facei];
            in.Tw = Tw[facei];
            in.adiabatic = adiabatic;

            sol[facei] = solveFace(in, coeffs);

            const scalar relax = coeffs.TvRelaxation;
            Tv[facei] =
                positive
                (
                    relax*sol[facei].TvNew + (1 - relax)*Tv[facei]
                );
            if (adiabatic)
            {
                Tw[facei] =
                    positive
                    (
                        relax*sol[facei].TwNew + (1 - relax)*Tw[facei]
                    );
            }
        }

        if (iter == nIter - 1)
        {
            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                const FaceOutput& s = sol[facei];

                data.tauw[facei] = s.tauw;
                data.A1[facei] = s.A1;
                data.yStar[facei] = s.yStar;
                data.yStarN[facei] = s.yStarN;
                data.qw[facei] = s.qw;
                data.Tw[facei] = Tw[facei];
                data.Tv[facei] = Tv[facei];
                data.Ath1[facei] = s.Ath1;
                data.Pk[facei] = s.Pk;
                data.epsBar[facei] = s.epsBar;
                data.epsP[facei] = s.epsP;
                data.dissAvg[facei] = s.dissAvg;
                // Discrete estimate of the same wall-normal dissipation
                // from the resolved cell gradient (wall-face part is zero
                // since U = 0 at the wall)
                data.dissLocal[facei] =
                    muEff[celli]*Un[facei]*dUtdn[facei]
                   /positive(ynField[facei]);
                data.muw[facei] = muw[facei];
                data.muv[facei] = muv[facei];
                data.rhov[facei] = rhov[facei];
                data.cpw[facei] = cpw[facei];
                data.Pr[facei] = Pr[facei];
                data.bMu[facei] = s.bMu;
                data.Dcoef[facei] = s.D;
                data.Ccoef[facei] = s.C;
                data.DthCoef[facei] = s.Dth;
                data.Un[facei] = Un[facei];
                data.Tn[facei] = Tn[facei];
                data.yn[facei] = ynField[facei];
                data.thinCell[facei] = s.thin ? 1 : 0;

                data.uTau[facei] =
                    std::sqrt(mag(s.tauw)/positive(rhov[facei]));

                // Invert OpenFOAM's discrete stress relation with the
                // patch-field viscosity and density that the stress
                // operator actually uses (they differ from the physical
                // wall values under the Neumann coupling).
                const scalar tauwEff = max(s.tauw, scalar(0));
                if (magGradU[facei] > ROOTVSMALL)
                {
                    data.nutShear[facei] =
                        max
                        (
                            scalar(0),
                            (tauwEff - muPatch[facei]*magGradU[facei])
                           /(positive(rhoPatch[facei])*magGradU[facei])
                        );
                }
                else
                {
                    data.nutShear[facei] = 0;
                }
            }
        }
    }

    return data;
}

} // End namespace WangNeumannWallFunction
} // End namespace Foam

// ************************************************************************* //
