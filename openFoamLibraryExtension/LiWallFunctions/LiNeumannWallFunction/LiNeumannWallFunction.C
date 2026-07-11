#include "LiNeumannWallFunction.H"
#include "LiNeumannWallHeatFluxTemperatureFvPatchScalarField.H"

#include "dictionary.H"
#include "fvMesh.H"
#include "fvPatchField.H"
#include "IOobject.H"
#include "Ostream.H"
#include "processorPolyPatch.H"
#include "surfaceFields.H"
#include "volFields.H"

#include <cmath>

namespace Foam
{
namespace
{

inline scalar clampUnit(const scalar value)
{
    return max(scalar(-1), min(scalar(1), value));
}


inline scalar positive(const scalar value)
{
    return max(value, VSMALL);
}


inline vector tangential(const vector& value, const vector& normal)
{
    return value - normal*(value & normal);
}


scalar spaldingYPlus
(
    const scalar uPlus,
    const scalar kappa,
    const scalar B
)
{
    const scalar y0 = exp(-kappa*B);
    const scalar ku = kappa*uPlus;

    return
        uPlus
      + y0
       *(
            exp(min(scalar(50), ku))
          - (1 + ku + 0.5*sqr(ku) + (1.0/6.0)*pow3(ku))
        );
}


scalar yWhitePlus
(
    const scalar uPlus,
    const scalar beta,
    const scalar Gamma,
    const scalar kappa,
    const scalar B
)
{
    const scalar discriminant = sqr(beta) + 4*Gamma;

    if (Gamma <= ROOTVSMALL || discriminant <= ROOTVSMALL)
    {
        return exp(-kappa*B + min(scalar(50), kappa*uPlus));
    }

    const scalar Q = sqrt(discriminant);
    const scalar phi = asin(clampUnit(-beta/Q));
    const scalar arg = clampUnit((2*Gamma*uPlus - beta)/Q);
    const scalar exponent =
        min
        (
            scalar(50),
            max
            (
                scalar(-50),
                kappa/sqrt(Gamma)*(asin(arg) - phi)
            )
        );

    return exp(exponent - kappa*B);
}


scalar solveUPlusForYPlus
(
    const scalar yPlusTarget,
    const scalar uPlusUpper,
    const scalar beta,
    const scalar Gamma,
    const LiNeumannWallFunction::Coeffs& coeffs
)
{
    if (yPlusTarget <= SMALL || uPlusUpper <= SMALL)
    {
        return 0;
    }

    scalar lo = 0;
    scalar hi = max(uPlusUpper, scalar(1));

    for
    (
        label iter = 0;
        iter < 20
     && LiNeumannWallFunction::lawYPlus(hi, beta, Gamma, coeffs) < yPlusTarget;
        ++iter
    )
    {
        hi *= 2;
    }

    for (label iter = 0; iter < 80; ++iter)
    {
        const scalar mid = 0.5*(lo + hi);

        if (LiNeumannWallFunction::lawYPlus(mid, beta, Gamma, coeffs) < yPlusTarget)
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


scalar solveCoupledUPlus
(
    const scalar U1,
    const scalar rhoWall,
    const scalar muWall,
    const scalar y,
    const scalar hrg,
    const scalar hw,
    const scalar hRef,
    const scalar URef,
    const scalar rg,
    const scalar initialUPlus,
    const LiNeumannWallFunction::Coeffs& coeffs,
    scalar& finalResidual,
    label& iterations,
    label& converged
)
{
    const auto residual = [&](const scalar uPlus) -> scalar
    {
        const scalar uTau = U1/positive(uPlus);
        const scalar beta = uTau*(hrg - hw)/(positive(URef)*positive(hw));
        const scalar Gamma = rg*sqr(uTau)/(2*positive(hw));
        const scalar target = rhoWall*uTau*y/muWall;

        return LiNeumannWallFunction::lawYPlus
        (
            uPlus,
            beta,
            Gamma,
            coeffs
        ) - target;
    };

    scalar lo = max(SMALL, min(scalar(0.1), max(initialUPlus, scalar(1))*1e-6));
    scalar hi = max(scalar(1), initialUPlus);
    scalar fLo = residual(lo);
    scalar fHi = residual(hi);

    for (label iter = 0; iter < 80 && fLo > 0; ++iter)
    {
        hi = lo;
        fHi = fLo;
        lo *= 0.1;
        fLo = residual(lo);
    }

    for (label iter = 0; iter < 80 && fHi < 0; ++iter)
    {
        lo = hi;
        fLo = fHi;
        hi *= 2;
        fHi = residual(hi);
    }

    if (!std::isfinite(fLo) || !std::isfinite(fHi) || fLo > 0 || fHi < 0)
    {
        finalResidual = min(mag(fLo), mag(fHi));
        iterations = 0;
        converged = 0;
        return mag(fLo) < mag(fHi) ? lo : hi;
    }

    scalar mid = 0.5*(lo + hi);

    for (label iter = 0; iter < coeffs.maxIter; ++iter)
    {
        mid = 0.5*(lo + hi);
        const scalar fMid = residual(mid);

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

        finalResidual = mag(fMid);
        iterations = iter + 1;

        if
        (
            mag(hi - lo)/positive(mid) < coeffs.tolerance
         || finalResidual < coeffs.tolerance
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


namespace LiNeumannWallFunction
{

Coeffs::Coeffs()
:
    kappa(0.41),
    B(6.0),
    Prt(0.9),
    omegaH(0.72),
    cv1(7.1),
    enthalpyOffset(0),
    maxIter(80),
    tolerance(1e-8),
    thermalMode("auto"),
    strict(false),
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
    kappa(dict.getOrDefault<scalar>("kappa", 0.41)),
    B(dict.getOrDefault<scalar>("B", 6.0)),
    Prt(dict.getOrDefault<scalar>("Prt", 0.9)),
    omegaH(dict.getOrDefault<scalar>("omegaH", 0.72)),
    cv1(dict.getOrDefault<scalar>("cv1", 7.1)),
    enthalpyOffset(dict.getOrDefault<scalar>("enthalpyOffset", 0)),
    maxIter(dict.getOrDefault<label>("maxIter", 80)),
    tolerance(dict.getOrDefault<scalar>("tolerance", 1e-8)),
    thermalMode(dict.getOrDefault<word>("thermalMode", "auto")),
    strict(dict.getOrDefault<bool>("strict", false)),
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
                << "Unsupported wallTransport mode " << wallTransportMode << nl
                << "Supported values are fromPatch, sutherland, and constant"
                << exit(FatalIOError);
        }
    }
}


void Coeffs::writeEntries(Ostream& os) const
{
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa);
    os.writeEntryIfDifferent<scalar>("B", 6.0, B);
    os.writeEntryIfDifferent<scalar>("Prt", 0.9, Prt);
    os.writeEntryIfDifferent<scalar>("omegaH", 0.72, omegaH);
    os.writeEntryIfDifferent<scalar>("cv1", 7.1, cv1);
    os.writeEntryIfDifferent<scalar>("enthalpyOffset", 0, enthalpyOffset);
    os.writeEntryIfDifferent<label>("maxIter", 80, maxIter);
    os.writeEntryIfDifferent<scalar>("tolerance", 1e-8, tolerance);
    os.writeEntryIfDifferent<word>("thermalMode", word("auto"), thermalMode);
    os.writeEntryIfDifferent<bool>("strict", false, strict);
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
    uTau(size, Zero),
    tauw(size, Zero),
    yPlus(size, Zero),
    uPlus(size, Zero),
    recovery(size, Zero),
    rg(size, Zero),
    beta(size, Zero),
    Gamma(size, Zero),
    Pr(size, Zero),
    Tw(size, Zero),
    hw(size, Zero),
    h1(size, Zero),
    hRef(size, Zero),
    hrg(size, Zero),
    qw(size, Zero),
    qHalf(size, Zero),
    uHalf(size, Zero),
    nutShear(size, Zero),
    magGradU(size, Zero),
    tauDeficit(size, Zero),
    mut1(size, Zero),
    nuTilda1(size, Zero),
    residual(size, Zero),
    iterations(size, Zero),
    converged(size, Zero),
    valid(size, Zero),
    referenceCell(size, -1),
    status(size, Zero)
{}


scalar lawYPlus
(
    const scalar uPlus,
    const scalar beta,
    const scalar Gamma,
    const Coeffs& coeffs
)
{
    if (Gamma <= ROOTVSMALL && mag(beta) <= ROOTVSMALL)
    {
        return spaldingYPlus(uPlus, coeffs.kappa, coeffs.B);
    }

    const scalar y0 = exp(-coeffs.kappa*coeffs.B);
    const scalar ku = coeffs.kappa*uPlus;
    const scalar uGM =
        uPlus
       *(
            1
          + 0.5*coeffs.omegaH*beta*uPlus
          - (1.0/3.0)*coeffs.omegaH*Gamma*sqr(uPlus)
        );

    return
        uGM
      + yWhitePlus
        (
            uPlus,
            beta,
            Gamma,
            coeffs.kappa,
            coeffs.B
        )
      - y0*(1 + ku + 0.5*sqr(ku) + (1.0/6.0)*pow3(ku));
}


scalar dYWhiteDuPlus
(
    const scalar uPlus,
    const scalar beta,
    const scalar Gamma,
    const Coeffs& coeffs
)
{
    const scalar discriminant = sqr(beta) + 4*Gamma;

    if (Gamma <= ROOTVSMALL || discriminant <= ROOTVSMALL)
    {
        return
            coeffs.kappa
           *exp(-coeffs.kappa*coeffs.B + min(scalar(50), coeffs.kappa*uPlus));
    }

    const scalar Q = sqrt(discriminant);
    const scalar arg = clampUnit((2*Gamma*uPlus - beta)/Q);
    const scalar bracket = max(scalar(1e-14), 1 - sqr(arg));

    return
        2
       *yWhitePlus(uPlus, beta, Gamma, coeffs.kappa, coeffs.B)
       *coeffs.kappa*sqrt(Gamma)/Q
       /sqrt(bracket);
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
    if (isA<LiNeumannWallHeatFluxTemperatureFvPatchScalarField>(Tpatch))
    {
        return refCast<const LiNeumannWallHeatFluxTemperatureFvPatchScalarField>
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
    const volScalarField& p = thermo.p();
    const volScalarField& he = thermo.he();
    const bool internalEnergy = he.name() == "e";

    const fvPatchScalarField& Tpatch = thermo.T().boundaryField()[patchi];
    const bool adiabatic = isAdiabatic(coeffs, Tpatch);
    const bool neumannTBC =
        isA<LiNeumannWallHeatFluxTemperatureFvPatchScalarField>(Tpatch);

    // Physical wall temperature. Under the Neumann coupling the patch T
    // value drifts along the prescribed tangent and cannot supply Tw.
    scalarField Tw(Tpatch);
    if (!adiabatic)
    {
        if (coeffs.Tw > 0)
        {
            Tw = coeffs.Tw;
        }
        else if (neumannTBC)
        {
            FatalErrorInFunction
                << "Patch " << patch.name() << " uses "
                << "LiNeumannWallHeatFluxTemperature in isothermal mode "
                << "but no positive Tw was given in the wall-function "
                << "coefficients. The patch T value drifts under this "
                << "coupling and cannot supply the wall temperature."
                << nl << exit(FatalError);
        }
    }

    const scalarField pWall(p.boundaryField()[patchi]);

    const tmp<scalarField> tcpw = thermo.Cp(pWall, Tw, patchi);
    const scalarField& cpw = tcpw();
    const labelList cells(faceCells);
    const tmp<scalarField> trhow = thermo.rhoEoS(pWall, Tw, cells);
    const scalarField& rhoW = trhow();
    const tmp<scalarField> theWall = thermo.he(pWall, Tw, patchi);
    const scalarField& heWall = theWall();

    // Physical wall transport properties at (p_wall, Tw).
    scalarField muw(nFaces, Zero);
    scalarField lambdaW(nFaces, Zero);

    if (coeffs.wallTransportMode == "sutherland")
    {
        forAll(muw, facei)
        {
            const scalar T = max(Tw[facei], SMALL);
            muw[facei] =
                coeffs.wallTransportAs*T*sqrt(T)
               /(T + coeffs.wallTransportTs);
        }

        if (coeffs.wallTransportPr > 0)
        {
            lambdaW = muw*cpw/coeffs.wallTransportPr;
        }
        else
        {
            // Modified Eucken correction.
            const tmp<scalarField> tcvw = thermo.Cv(pWall, Tw, patchi);
            const scalarField& cvw = tcvw();
            forAll(lambdaW, facei)
            {
                const scalar cv = positive(cvw[facei]);
                const scalar R = max(cpw[facei] - cv, scalar(0));
                lambdaW[facei] = muw[facei]*cv*(1.32 + 1.77*R/cv);
            }
        }
    }
    else if (coeffs.wallTransportMode == "constant")
    {
        muw = coeffs.wallTransportMu;
        lambdaW = coeffs.wallTransportKappa;
    }
    else
    {
        if (neumannTBC && !adiabatic)
        {
            FatalErrorInFunction
                << "Patch " << patch.name() << " uses "
                << "LiNeumannWallHeatFluxTemperature in isothermal mode "
                << "with wallTransport mode fromPatch. The patch transport "
                << "properties are evaluated at the drifted patch T and are "
                << "not the wall values; select wallTransport sutherland "
                << "or constant." << nl << exit(FatalError);
        }

        muw = muPatch;
        const tmp<scalarField> tlambdaW = turbModel.kappa(patchi);
        lambdaW = tlambdaW();
    }

    scalarField hsWall(heWall);
    if (internalEnergy)
    {
        forAll(hsWall, facei)
        {
            hsWall[facei] += pWall[facei]/positive(rhoW[facei]);
        }
    }

    // Cache values across processor boundaries for second-cell selection.
    const label nPatches = mesh.boundary().size();
    List<vectorField> procC(nPatches);
    List<vectorField> procU(nPatches);
    List<scalarField> procP(nPatches);
    List<scalarField> procRho(nPatches);
    List<scalarField> procHe(nPatches);

    forAll(mesh.boundary(), proci)
    {
        if (isA<processorPolyPatch>(mesh.boundaryMesh()[proci]))
        {
            const processorPolyPatch& pp =
                refCast<const processorPolyPatch>(mesh.boundaryMesh()[proci]);
            procC[proci] = pp.neighbFaceCellCentres();
            procU[proci] = U.boundaryField()[proci].patchNeighbourField();
            procP[proci] = p.boundaryField()[proci].patchNeighbourField();
            procRho[proci] = rho.boundaryField()[proci].patchNeighbourField();
            procHe[proci] = he.boundaryField()[proci].patchNeighbourField();
        }
    }

    label invalidCount = 0;

    forAll(faceCells, facei)
    {
        const label celli = faceCells[facei];
        const vector normal = nOut[facei];
        const vector inward = -normal;
        const vector wallVelocity = Uw[facei];
        const scalar U1 = mag(tangential(U[celli] - wallVelocity, normal));
        const scalar rho1 = positive(rho[celli]);
        const scalar mu1 = positive(mu[celli]);
        const scalar muWall = positive(muw[facei]);
        const scalar muStress = positive(muPatch[facei]);
        const scalar rhoWall = positive(rhoW[facei]);
        const scalar rhoStress = positive(rhoPatch[facei]);
        const scalar cpWall = positive(cpw[facei]);
        const scalar lambdaWall = positive(lambdaW[facei]);
        const scalar yFace = positive(y[facei]);
        const scalar hw = positive(cpWall*positive(Tw[facei]) + coeffs.enthalpyOffset);

        scalar hs1 = he[celli];
        if (internalEnergy)
        {
            hs1 += p[celli]/rho1;
        }
        const scalar h1 = hw + hs1 - hsWall[facei];

        data.Tw[facei] = Tw[facei];
        data.hw[facei] = hw;
        data.h1[facei] = h1;

        if (U1 <= SMALL)
        {
            data.valid[facei] = 1;
            data.converged[facei] = 1;
            data.status[facei] = 1;
            continue;
        }

        scalar bestScore = -GREAT;
        scalar URef = 0;
        scalar hsRef = 0;
        label refCell = -1;

        const labelList& neighbours = mesh.cellCells()[celli];
        forAll(neighbours, nbrI)
        {
            const label nbr = neighbours[nbrI];
            const scalar score = (mesh.C()[nbr] - mesh.C()[celli]) & inward;

            if (score > bestScore)
            {
                bestScore = score;
                URef = mag(tangential(U[nbr] - wallVelocity, normal));
                hsRef = he[nbr];
                if (internalEnergy)
                {
                    hsRef += p[nbr]/positive(rho[nbr]);
                }
                refCell = nbr;
            }
        }

        const cell& ownerFaces = mesh.cells()[celli];
        forAll(ownerFaces, ownerFaceI)
        {
            const label meshFacei = ownerFaces[ownerFaceI];
            if (meshFacei < mesh.nInternalFaces())
            {
                continue;
            }

            const label proci = mesh.boundaryMesh().whichPatch(meshFacei);
            if
            (
                proci < 0
             || !isA<processorPolyPatch>(mesh.boundaryMesh()[proci])
            )
            {
                continue;
            }

            const label localFacei = meshFacei - mesh.boundaryMesh()[proci].start();
            const scalar score =
                (procC[proci][localFacei] - mesh.C()[celli]) & inward;

            if (score > bestScore)
            {
                bestScore = score;
                URef = mag
                (
                    tangential
                    (
                        procU[proci][localFacei] - wallVelocity,
                        normal
                    )
                );
                hsRef = procHe[proci][localFacei];
                if (internalEnergy)
                {
                    hsRef +=
                        procP[proci][localFacei]
                       /positive(procRho[proci][localFacei]);
                }
                refCell = -meshFacei - 2;
            }
        }

        if (bestScore <= SMALL || URef <= U1*(1 + 1e-6))
        {
            data.status[facei] = 2;
            ++invalidCount;
            continue;
        }

        const scalar hRef = hw + hsRef - hsWall[facei];
        const scalar x = U1/URef;
        scalar hrg = hw;
        scalar rg = 2*(hw - hRef)/sqr(URef);

        if (!adiabatic)
        {
            const scalar denominator = x - sqr(x);
            if (mag(denominator) <= ROOTVSMALL)
            {
                data.status[facei] = 3;
                ++invalidCount;
                continue;
            }

            hrg =
                (h1 - hw*(1 - x) - hRef*sqr(x))/denominator;
            rg = 2*(hrg - hRef)/sqr(URef);
        }

        data.hRef[facei] = hRef;
        data.hrg[facei] = hrg;
        data.rg[facei] = rg;
        data.referenceCell[facei] = refCell;

        if
        (
            !std::isfinite(hw)
         || !std::isfinite(hRef)
         || !std::isfinite(hrg)
         || !std::isfinite(rg)
         || rg <= ROOTVSMALL
        )
        {
            data.status[facei] = 4;
            ++invalidCount;
            continue;
        }

        const vector gradU = tangential(gradUw[facei], normal);
        const scalar magGradU = mag(gradU);
        scalar initialUTau =
            sqrt(max((nutw[facei] + nuw[facei])*magGradU, VSMALL));
        if (initialUTau <= SMALL)
        {
            initialUTau = sqrt(max(nuw[facei]*U1/yFace, VSMALL));
        }

        scalar finalResidual = GREAT;
        label iterations = 0;
        label converged = 0;
        const scalar uPlus = solveCoupledUPlus
        (
            U1,
            rhoWall,
            muWall,
            yFace,
            hrg,
            hw,
            hRef,
            URef,
            rg,
            U1/positive(initialUTau),
            coeffs,
            finalResidual,
            iterations,
            converged
        );

        if (!converged || uPlus <= SMALL || !std::isfinite(uPlus))
        {
            data.status[facei] = 5;
            data.residual[facei] = finalResidual;
            data.iterations[facei] = iterations;
            ++invalidCount;
            continue;
        }

        const scalar uTau = U1/uPlus;
        const scalar beta = uTau*(hrg - hw)/(URef*hw);
        const scalar Gamma = rg*sqr(uTau)/(2*hw);
        const scalar yPlus = rhoWall*uTau*yFace/muWall;
        const scalar tauw = rhoWall*sqr(uTau);
        const scalar Pr = positive(cpWall*muWall/lambdaWall);
        const scalar recovery = pow(Pr, 1.0/3.0);
        const scalar qw =
            adiabatic
          ? scalar(0)
          :
            (
                rg - 2*(hw - hRef)/sqr(URef)
            )*URef*tauw/(2*Pr);

        const scalar uPlusHalf = solveUPlusForYPlus
        (
            0.5*yPlus,
            uPlus,
            beta,
            Gamma,
            coeffs
        );
        const scalar uHalf = uTau*uPlusHalf;
        const scalar qHalf =
            qw - recovery*uHalf*tauw/positive(coeffs.Prt);

        const scalar y0 = exp(-coeffs.kappa*coeffs.B);
        const scalar ku = coeffs.kappa*uPlus;
        const scalar mutByMu =
            1
          + coeffs.omegaH*beta*uPlus
          - coeffs.omegaH*Gamma*sqr(uPlus)
          + dYWhiteDuPlus(uPlus, beta, Gamma, coeffs)
          - coeffs.kappa*y0*(1 + ku + 0.5*sqr(ku))
          - mu1/muWall;
        const scalar mut1 = max(scalar(0), muWall*mutByMu);
        const scalar nuTilda1 = solveNuTilda
        (
            mut1/rho1,
            mu1/rho1,
            coeffs.cv1,
            coeffs.maxIter,
            coeffs.tolerance
        );

        // Invert OpenFOAM's discrete stress relation with the patch-field
        // viscosity and density that the stress operator actually uses;
        // under the Neumann coupling these differ from the physical wall
        // values because the patch T drifts.
        scalar nutShear = 0;
        scalar tauDeficit = 0;
        if (magGradU > ROOTVSMALL)
        {
            nutShear =
                max
                (
                    scalar(0),
                    (tauw - muStress*magGradU)/(rhoStress*magGradU)
                );
            // Negative where nut >= 0 cannot reach the target: the patch
            // viscosity alone already delivers more stress than tauw, so
            // the clipped nut over-delivers by exactly this margin.
            tauDeficit = min(scalar(0), tauw - muStress*magGradU);
        }
        else
        {
            nutShear = mut1/rhoStress;
        }

        data.uTau[facei] = uTau;
        data.tauw[facei] = tauw;
        data.yPlus[facei] = yPlus;
        data.uPlus[facei] = uPlus;
        data.recovery[facei] = recovery;
        data.beta[facei] = beta;
        data.Gamma[facei] = Gamma;
        data.Pr[facei] = Pr;
        data.qw[facei] = qw;
        data.qHalf[facei] = qHalf;
        data.uHalf[facei] = uHalf;
        data.nutShear[facei] = nutShear;
        data.magGradU[facei] = magGradU;
        data.tauDeficit[facei] = tauDeficit;
        data.mut1[facei] = mut1;
        data.nuTilda1[facei] = nuTilda1;
        data.residual[facei] = finalResidual;
        data.iterations[facei] = iterations;
        data.converged[facei] = converged;
        data.valid[facei] = 1;
        data.status[facei] = 0;
    }

    if (invalidCount)
    {
        if (coeffs.strict)
        {
            FatalErrorInFunction
                << invalidCount << " face(s) on patch " << patch.name()
                << " could not evaluate the Li wall function. "
                << "Inspect status diagnostics; strict=false enables the "
                << "laminar fallback." << nl
                << exit(FatalError);
        }

        static bool warned = false;
        if (!warned)
        {
            WarningInFunction
                << "Li (Neumann) wall function used its laminar fallback on "
                << invalidCount << " face(s) of patch " << patch.name()
                << ". Enable diagnostics or strict=true to locate invalid "
                << "second-cell/reference states." << endl;
            warned = true;
        }
    }

    return data;
}

} // End namespace LiNeumannWallFunction
} // End namespace Foam

// ************************************************************************* //
