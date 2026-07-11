#include "GaoWallFunction.H"

#include "dictionary.H"
#include "fvPatchField.H"
#include "Ostream.H"
#include "surfaceFields.H"
#include "volFields.H"

#include <cmath>

namespace Foam
{
namespace
{

inline scalar clampUnit(const scalar x)
{
    return max(scalar(-1), min(scalar(1), x));
}


inline scalar positive(const scalar x)
{
    return max(x, VSMALL);
}


inline bool isAdiabaticMode
(
    const word& mode,
    const fvPatchScalarField& Tpatch
)
{
    if (mode == "adiabatic")
    {
        return true;
    }

    if (mode == "isothermal")
    {
        return false;
    }

    return !Tpatch.fixesValue();
}


scalar yWhitePlus
(
    const scalar uPlus,
    const scalar beta,
    const scalar Gamma,
    const scalar kappa,
    const scalar y0
)
{
    const scalar GammaEff = positive(Gamma);
    const scalar Q = positive(sqrt(sqr(beta) + 4*GammaEff));
    const scalar phi = asin(clampUnit(-beta/Q));
    const scalar arg = clampUnit((2*GammaEff*uPlus - beta)/Q);

    const scalar exponent =
        min
        (
            scalar(50),
            max(scalar(-50), kappa/sqrt(GammaEff)*(asin(arg) - phi))
        );

    return exp(exponent)*y0;
}


scalar lawYPlus
(
    const scalar uPlus,
    const scalar beta,
    const scalar Gamma,
    const scalar kappa,
    const scalar B
)
{
    const scalar y0 = exp(-kappa*B);
    const scalar ku = kappa*uPlus;

    return
        uPlus
      + yWhitePlus(uPlus, beta, Gamma, kappa, y0)
      - y0*(1 + ku + 0.5*sqr(ku) + (1.0/6.0)*pow3(ku));
}


scalar betaFromTemperature
(
    const bool adiabatic,
    const scalar Tcell,
    const scalar Tw,
    const scalar uPlus,
    const scalar Gamma
)
{
    if (adiabatic || uPlus <= SMALL)
    {
        return 0;
    }

    return (Tcell/positive(Tw) - 1 + Gamma*sqr(uPlus))/uPlus;
}


scalar solveUPlusForYPlus
(
    const scalar yPlusTarget,
    const scalar uPlusCell,
    const scalar beta,
    const scalar Gamma,
    const scalar kappa,
    const scalar B
)
{
    if (yPlusTarget <= SMALL || uPlusCell <= SMALL)
    {
        return 0;
    }

    scalar lo = 0;
    scalar hi = max(uPlusCell, SMALL);

    for (label i = 0; i < 12 && lawYPlus(hi, beta, Gamma, kappa, B) < yPlusTarget; ++i)
    {
        hi *= 2;
    }

    for (label i = 0; i < 60; ++i)
    {
        const scalar mid = 0.5*(lo + hi);

        if (lawYPlus(mid, beta, Gamma, kappa, B) < yPlusTarget)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    return 0.5*(lo + hi);
}


scalar solveUPlus
(
    const bool adiabatic,
    const scalar U,
    const scalar yPlusUTau,
    const scalar Tcell,
    const scalar Tw,
    const scalar r,
    const scalar cp,
    const scalar kappa,
    const scalar B,
    const scalar initialUPlus,
    const label maxIter,
    const scalar tolerance,
    scalar& residual,
    label& iterations,
    label& converged
)
{
    if (U <= SMALL || yPlusUTau <= SMALL)
    {
        residual = 0;
        iterations = 0;
        converged = 1;
        return 0;
    }

    const auto lawResidual = [&](const scalar uPlus) -> scalar
    {
        const scalar uTau = U/positive(uPlus);
        const scalar Gamma = r*sqr(uTau)/(2*cp*Tw);
        const scalar beta =
            betaFromTemperature(adiabatic, Tcell, Tw, uPlus, Gamma);

        return
            lawYPlus(uPlus, beta, Gamma, kappa, B)
          - yPlusUTau/uPlus;
    };

    scalar lo = max(SMALL, min(scalar(0.1), max(initialUPlus, scalar(1))*1e-6));
    scalar hi = max(scalar(1), initialUPlus);

    scalar fLo = lawResidual(lo);
    scalar fHi = lawResidual(hi);

    for (label i = 0; i < 80 && fLo > 0; ++i)
    {
        hi = lo;
        fHi = fLo;
        lo *= 0.1;
        fLo = lawResidual(lo);
    }

    for (label i = 0; i < 80 && fHi < 0; ++i)
    {
        lo = hi;
        fLo = fHi;
        hi *= 2;
        fHi = lawResidual(hi);
    }

    if (fLo > 0 || fHi < 0)
    {
        residual = min(mag(fLo), mag(fHi));
        iterations = 0;
        converged = 0;
        return mag(fLo) < mag(fHi) ? lo : hi;
    }

    scalar mid = 0.5*(lo + hi);

    for (label i = 0; i < maxIter; ++i)
    {
        mid = 0.5*(lo + hi);
        const scalar fMid = lawResidual(mid);

        if (fMid > 0)
        {
            hi = mid;
            fHi = fMid;
        }
        else
        {
            lo = mid;
            fLo = fMid;
        }

        residual = mag(fMid);
        iterations = i + 1;

        if
        (
            mag(hi - lo)/positive(mid) < tolerance
         || residual < tolerance
        )
        {
            converged = 1;
            return mid;
        }
    }

    converged = 0;
    return 0.5*(lo + hi);
}

} // End anonymous namespace


namespace GaoWallFunction
{

Coeffs::Coeffs()
:
    kappa(0.41),
    B(5.0),
    Prt(0.9),
    Cmu(0.09),
    maxIter(50),
    tolerance(1e-7),
    thermalMode("auto")
{}


Coeffs::Coeffs(const dictionary& dict)
:
    kappa(dict.getOrDefault<scalar>("kappa", 0.41)),
    B(dict.getOrDefault<scalar>("B", 5.0)),
    Prt(dict.getOrDefault<scalar>("Prt", 0.9)),
    Cmu(dict.getOrDefault<scalar>("Cmu", 0.09)),
    maxIter(dict.getOrDefault<label>("maxIter", 50)),
    tolerance(dict.getOrDefault<scalar>("tolerance", 1e-7)),
    thermalMode(dict.getOrDefault<word>("thermalMode", "auto"))
{}


void Coeffs::writeEntries(Ostream& os) const
{
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa);
    os.writeEntryIfDifferent<scalar>("B", 5.0, B);
    os.writeEntryIfDifferent<scalar>("Prt", 0.9, Prt);
    os.writeEntryIfDifferent<scalar>("Cmu", 0.09, Cmu);
    os.writeEntryIfDifferent<label>("maxIter", 50, maxIter);
    os.writeEntryIfDifferent<scalar>("tolerance", 1e-7, tolerance);
    os.writeEntryIfDifferent<word>("thermalMode", word("auto"), thermalMode);
}


PatchData::PatchData(const label size)
:
    uTau(size, Zero),
    tauw(size, Zero),
    yPlus(size, Zero),
    uPlus(size, Zero),
    Gamma(size, Zero),
    beta(size, Zero),
    Pr(size, Zero),
    recovery(size, Zero),
    Tw(size, Zero),
    cp(size, Zero),
    rhoW(size, Zero),
    qw(size, Zero),
    qHalf(size, Zero),
    uHalf(size, Zero),
    nutShear(size, Zero),
    mutEq8(size, Zero),
    omega(size, Zero),
    k(size, Zero),
    residual(size, Zero),
    iterations(size, Zero),
    converged(size, Zero)
{}


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

    const fluidThermo& thermo = turbModel.transport();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradU(mag(Uw.snGrad()));

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const tmp<volScalarField> tmu = turbModel.mu();
    const volScalarField& mu = tmu();

    const volScalarField& rho = turbModel.rho();
    const fvPatchScalarField& rhoPatch = rho.boundaryField()[patchi];

    const tmp<scalarField> tkw = turbModel.kappa(patchi);
    const scalarField& kw = tkw();

    const fvPatchScalarField& Tpatch = thermo.T().boundaryField()[patchi];
    const scalarField Tcell(Tpatch.patchInternalField());

    const fvPatchScalarField& ppatch = thermo.p().boundaryField()[patchi];
    const scalarField pWall(ppatch);

    const bool adiabatic = isAdiabaticMode(coeffs.thermalMode, Tpatch);

    scalarField Tw(nFaces, Zero);
    forAll(Tw, facei)
    {
        Tw[facei] = adiabatic ? positive(Tcell[facei]) : positive(Tpatch[facei]);
    }

    scalarField cpw(nFaces, Zero);

    for (label iter = 0; iter < 3; ++iter)
    {
        const tmp<scalarField> tcp = thermo.Cp(pWall, Tw, patchi);
        const scalarField& cpTmp = tcp();

        forAll(Tw, facei)
        {
            cpw[facei] = positive(cpTmp[facei]);

            const scalar Pr =
                positive(cpw[facei]*positive(muw[facei])/positive(kw[facei]));
            const scalar r = pow(Pr, 1.0/3.0);

            if (adiabatic)
            {
                Tw[facei] =
                    positive(Tcell[facei] + r*sqr(magUp[facei])/(2*cpw[facei]));
            }
        }
    }

    const tmp<scalarField> tcp = thermo.Cp(pWall, Tw, patchi);
    const scalarField& cpFinal = tcp();

    const labelList cells(faceCells);
    const tmp<scalarField> trhow = thermo.rhoEoS(pWall, Tw, cells);
    const scalarField& rhoW = trhow();

    const scalar y0 = exp(-coeffs.kappa*coeffs.B);
    const scalar sqrtCmu = sqrt(coeffs.Cmu);

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];

        const scalar muWall = positive(muw[facei]);
        const scalar rhoWall = positive(rhoW[facei]);
        const scalar rhoStress = positive(rhoPatch[facei]);
        const scalar rhoCell = positive(rho[celli]);
        const scalar kThermal = positive(kw[facei]);
        const scalar cpFace = positive(cpFinal[facei]);
        const scalar TwFace = positive(Tw[facei]);
        const scalar yFace = positive(y[facei]);
        const scalar U = magUp[facei];

        const scalar Pr = positive(cpFace*muWall/kThermal);
        const scalar r = pow(Pr, 1.0/3.0);

        scalar ut =
            sqrt(max((nutw[facei] + nuw[facei])*magGradU[facei], VSMALL));

        if (ut <= SMALL && U > SMALL)
        {
            ut = sqrt(max(nuw[facei]*U/yFace, VSMALL));
        }

        scalar finalResidual = GREAT;
        label nIter = 0;
        label converged = 0;

        const scalar uPlus =
            solveUPlus
            (
                adiabatic,
                U,
                rhoWall*U*yFace/muWall,
                Tcell[facei],
                TwFace,
                r,
                cpFace,
                coeffs.kappa,
                coeffs.B,
                U/positive(ut),
                coeffs.maxIter,
                coeffs.tolerance,
                finalResidual,
                nIter,
                converged
            );

        ut = U/positive(uPlus);
        data.iterations[facei] = nIter;
        data.converged[facei] = converged;

        const scalar yPlus = rhoWall*ut*yFace/muWall;
        const scalar Gamma = r*sqr(ut)/(2*cpFace*TwFace);
        const scalar beta =
            betaFromTemperature
            (
                adiabatic,
                Tcell[facei],
                TwFace,
                uPlus,
                Gamma
            );

        const scalar tauw = rhoWall*sqr(ut);
        const scalar qw =
            adiabatic
          ? scalar(0)
          : beta*rhoWall*TwFace*kThermal*ut/(1.2*muWall);

        const scalar uPlusHalf =
            solveUPlusForYPlus
            (
                0.5*yPlus,
                uPlus,
                beta,
                Gamma,
                coeffs.kappa,
                coeffs.B
            );

        const scalar uHalf = ut*uPlusHalf;
        const scalar qHalf = qw - r*uHalf*tauw/positive(coeffs.Prt);

        const scalar yWhite =
            yWhitePlus(uPlus, beta, Gamma, coeffs.kappa, y0);
        const scalar GammaEff = positive(Gamma);
        const scalar Q = positive(sqrt(sqr(beta) + 4*GammaEff));
        const scalar arg = clampUnit((2*GammaEff*uPlus - beta)/Q);
        const scalar bracket = max(scalar(1e-12), 1 - sqr(arg));
        const scalar dyWhiteDyPlus =
            2*yWhite*(coeffs.kappa*sqrt(GammaEff)/Q)
           /sqrt(bracket);

        const scalar ku = coeffs.kappa*uPlus;
        const scalar mutByMu =
            1
          + dyWhiteDyPlus
          - coeffs.kappa*y0*(1 + ku + 0.5*sqr(ku))
          - positive(mu[celli])/muWall;

        const scalar mutEq8 = max(scalar(0), muWall*mutByMu);

        scalar nutShear = 0;
        if (magGradU[facei] > ROOTVSMALL)
        {
            nutShear =
                max
                (
                    scalar(0),
                    (tauw - muWall*magGradU[facei])
                   /(rhoStress*magGradU[facei])
                );
        }
        else
        {
            nutShear = mutEq8/rhoStress;
        }

        const scalar omegaI =
            6*muWall/(0.075*rhoWall*sqr(yFace));
        const scalar omegaO =
            ut/(sqrtCmu*coeffs.kappa*yFace);
        const scalar omega = sqrt(sqr(omegaI) + sqr(omegaO));
        const scalar k = max(scalar(0), omega*mutEq8/rhoCell);

        data.uTau[facei] = ut;
        data.tauw[facei] = tauw;
        data.yPlus[facei] = yPlus;
        data.uPlus[facei] = uPlus;
        data.Gamma[facei] = Gamma;
        data.beta[facei] = beta;
        data.Pr[facei] = Pr;
        data.recovery[facei] = r;
        data.Tw[facei] = TwFace;
        data.cp[facei] = cpFace;
        data.rhoW[facei] = rhoWall;
        data.qw[facei] = qw;
        data.qHalf[facei] = qHalf;
        data.uHalf[facei] = uHalf;
        data.nutShear[facei] = nutShear;
        data.mutEq8[facei] = mutEq8;
        data.omega[facei] = omega;
        data.k[facei] = k;
        data.residual[facei] = finalResidual;
    }

    return data;
}

} // End namespace GaoWallFunction
} // End namespace Foam

// ************************************************************************* //
