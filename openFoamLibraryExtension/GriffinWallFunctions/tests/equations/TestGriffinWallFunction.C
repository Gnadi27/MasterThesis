#include "GriffinWallFunction.H"
#include "GriffinBoundaryLayerEdge.H"
#include "dictionary.H"
#include "IOstreams.H"

#include <cmath>

using namespace Foam;

namespace
{

bool close(const scalar actual, const scalar expected, const scalar tolerance)
{
    return mag(actual - expected) <= tolerance*max(scalar(1), mag(expected));
}


void require(const bool condition, const word& message, label& failures)
{
    if (!condition)
    {
        Info<< "FAIL: " << message << nl;
        ++failures;
    }
    else
    {
        Info<< "PASS: " << message << nl;
    }
}


// Independent fine-grid RK4 integration of the incompressible
// constant-property law of the wall, dU/dy = (tau_w/mu)/(1 + mut+[y*]),
// which the inverse transformation must reduce to when rho and mu are
// uniform.
scalar incompressibleU
(
    const scalar tauw,
    const scalar rho,
    const scalar mu,
    const scalar yM,
    const scalar kappa,
    const scalar APlus
)
{
    const auto dUdy = [&](const scalar y) -> scalar
    {
        const scalar yStar = y*Foam::sqrt(rho*tauw)/mu;
        return
            (tauw/mu)
           /(1 + GriffinWallFunction::mutPlusIncompressible(yStar, kappa, APlus));
    };

    const label n = 20000;
    const scalar h = yM/n;
    scalar U = 0;
    for (label i = 0; i < n; ++i)
    {
        const scalar y = i*h;
        const scalar k1 = dUdy(y);
        const scalar k2 = dUdy(y + 0.5*h);
        const scalar k4 = dUdy(y + h);
        U += h/6*(k1 + 4*k2 + k4);
    }
    return U;
}

} // End anonymous namespace


int main()
{
    label failures = 0;

    GriffinWallFunction::Coeffs coeffs;

    // ---- Calibration table (paper Sec. 2) -------------------------------
    {
        scalar kappa = 0, APlus = 0;

        GriffinWallFunction::calibration(182, kappa, APlus);
        require
        (
            close(kappa, 0.400, 1e-12) && close(APlus, 18.2, 1e-12),
            "calibration at Re_tau* = 182",
            failures
        );

        GriffinWallFunction::calibration(543, kappa, APlus);
        require
        (
            close(kappa, 0.408, 1e-12) && close(APlus, 17.4, 1e-12),
            "calibration at Re_tau* = 543",
            failures
        );

        // Linear interpolation midway between 543 and 1000.
        GriffinWallFunction::calibration(0.5*(543 + 1000), kappa, APlus);
        require
        (
            close(kappa, 0.5*(0.408 + 0.400), 1e-12)
         && close(APlus, 0.5*(17.4 + 17.0), 1e-12),
            "calibration linear interpolation",
            failures
        );

        // Constant extrapolation outside the table.
        GriffinWallFunction::calibration(50, kappa, APlus);
        require
        (
            close(kappa, 0.400, 1e-12) && close(APlus, 18.2, 1e-12),
            "calibration constant extrapolation below",
            failures
        );

        GriffinWallFunction::calibration(2e4, kappa, APlus);
        require
        (
            close(kappa, 0.391, 1e-12) && close(APlus, 16.5, 1e-12),
            "calibration constant extrapolation above",
            failures
        );
    }

    // ---- Van Driest-damped mixing length (paper Eq. 2.2) ----------------
    {
        require
        (
            GriffinWallFunction::mutPlusIncompressible(0, 0.41, 17) == 0,
            "mut+ vanishes at the wall",
            failures
        );

        const scalar yStar = 3000;
        require
        (
            close
            (
                GriffinWallFunction::mutPlusIncompressible(yStar, 0.41, 17),
                0.41*yStar,
                1e-3
            ),
            "mut+ approaches kappa y* in the log layer",
            failures
        );
    }

    // ---- Temperature-velocity relation (paper Eq. 2.8) ------------------
    {
        const scalar UM = 600, TM = 250, Ue = 900, Tw = 300, Tr = 421;
        const scalar s = 1.14, Pr = 0.72;

        require
        (
            close
            (
                GriffinWallFunction::temperatureFromU
                (
                    0, UM, TM, Ue, Tw, Tr, s, Pr
                ),
                Tw,
                1e-12
            ),
            "Eq. 2.8 endpoint U = 0 -> T_w",
            failures
        );

        require
        (
            close
            (
                GriffinWallFunction::temperatureFromU
                (
                    UM, UM, TM, Ue, Tw, Tr, s, Pr
                ),
                TM,
                1e-12
            ),
            "Eq. 2.8 endpoint U = U_m -> T_m",
            failures
        );

        // dT/dU against finite difference.
        const scalar U = 350;
        const scalar dU = 1e-4;
        const scalar numerical =
            (
                GriffinWallFunction::temperatureFromU
                (
                    U + dU, UM, TM, Ue, Tw, Tr, s, Pr
                )
              - GriffinWallFunction::temperatureFromU
                (
                    U - dU, UM, TM, Ue, Tw, Tr, s, Pr
                )
            )/(2*dU);
        const scalar analytic =
            GriffinWallFunction::dTemperaturedU
            (
                U, UM, TM, Ue, Tw, Tr, s, Pr
            );
        require
        (
            close(numerical, analytic, 1e-7),
            "Eq. 2.8 dT/dU agrees with finite difference",
            failures
        );
    }

    // ---- Recovery temperature -------------------------------------------
    {
        const scalar Te = 60, Ue = 900, Cp = 1005, Pr = 0.72;
        require
        (
            close
            (
                GriffinWallFunction::recoveryTemperature(Te, Ue, Cp, Pr),
                Te + Foam::cbrt(Pr)*sqr(Ue)/(2*Cp),
                1e-12
            ),
            "recovery temperature T_r = T_e + Pr^(1/3) U_e^2/(2 Cp)",
            failures
        );
    }

    // ---- Viscosity law and its derivative --------------------------------
    {
        GriffinWallFunction::Coeffs c;
        const scalar T = 231.5;
        const scalar dT = 1e-3;
        const scalar numerical = (c.mu(T + dT) - c.mu(T - dT))/(2*dT);
        require
        (
            close(numerical, c.dmudT(T), 1e-6),
            "Sutherland dmu/dT agrees with finite difference",
            failures
        );

        c.viscosityMode = "powerLaw";
        c.viscosityMuRef = 1.8e-5;
        c.viscosityTRef = 273;
        c.viscosityN = 0.75;
        const scalar numericalPL = (c.mu(T + dT) - c.mu(T - dT))/(2*dT);
        require
        (
            close(numericalPL, c.dmudT(T), 1e-6),
            "power-law dmu/dT agrees with finite difference",
            failures
        );
    }

    // ---- SA nuTilda inversion --------------------------------------------
    {
        const scalar nu = 1.7e-5;
        const scalar nutTarget = 3.2e-4;
        const scalar nuTilda = GriffinWallFunction::solveNuTilda
        (
            nutTarget,
            nu,
            coeffs.cv1,
            100,
            1e-12
        );
        const scalar chi = nuTilda/nu;
        const scalar fv1 = pow3(chi)/(pow3(chi) + pow3(coeffs.cv1));
        require
        (
            close(nuTilda*fv1, nutTarget, 1e-10),
            "SA nuTilda inversion round trip",
            failures
        );
    }

    // ---- Incompressible limit of the inverse transformation --------------
    // With uniform temperature (huge Cp kills the recovery term, constant
    // viscosity) the transformation must reduce to the incompressible law
    // of the wall; the shooting solver must then recover a manufactured
    // tau_w from the velocity it implies.
    {
        GriffinWallFunction::Coeffs c;
        c.viscosityMode = "constant";
        c.viscosityMu = 1.8e-5;
        c.integrationNodes = 96;
        c.tolerance = 1e-9;
        c.maxIter = 100;

        const scalar rho = 1.2;
        const scalar mu = c.viscosityMu;
        const scalar T0 = 300;
        const scalar delta = 0.05;
        const scalar yM = 0.004;
        const scalar tauwTrue = 0.36;   // uTau ~ 0.55 m/s

        const scalar ReTauStar = Foam::sqrt(rho*tauwTrue)*delta/mu;
        scalar kappa = 0, APlus = 0;
        GriffinWallFunction::calibration(ReTauStar, kappa, APlus);

        const scalar UM =
            incompressibleU(tauwTrue, rho, mu, yM, kappa, APlus);

        GriffinWallFunction::InnerInput in;
        in.yM = yM;
        in.UM = UM;
        in.TM = T0;
        in.rhoM = rho;
        in.Tw = T0;
        in.adiabatic = false;
        in.Ue = 2*UM;
        in.Te = T0;
        in.delta = delta;
        in.Cp = 1e12;   // recovery term negligible -> uniform T
        in.tauwGuess = 0.5*tauwTrue;

        const GriffinWallFunction::InnerResult res =
            GriffinWallFunction::solveInner(in, c);

        require(res.converged == 1, "incompressible shooting converged", failures);
        require
        (
            close(res.tauw, tauwTrue, 2e-3),
            "incompressible limit recovers manufactured tau_w",
            failures
        );
        require
        (
            close(res.ReTauStar, ReTauStar, 2e-3),
            "semi-local friction Reynolds number consistent",
            failures
        );

        // Grid convergence: doubling the nodes must not move tau_w.
        GriffinWallFunction::Coeffs cFine(c);
        cFine.integrationNodes = 192;
        const GriffinWallFunction::InnerResult resFine =
            GriffinWallFunction::solveInner(in, cFine);
        require
        (
            close(resFine.tauw, res.tauw, 5e-3),
            "tau_w grid-converged",
            failures
        );
    }

    // ---- Compressible cooled-wall case ------------------------------------
    // Hypersonic-like state: the wall is colder than recovery, so the model
    // must extract heat from the fluid (q_w < 0 in the OpenFOAM/Griffin
    // boundary-flux sign) and converge cleanly.
    {
        GriffinWallFunction::Coeffs c;
        c.integrationNodes = 64;

        GriffinWallFunction::InnerInput in;
        in.yM = 2e-4;
        in.UM = 550;
        in.TM = 210;
        in.rhoM = 0.12;
        in.Tw = 300;
        in.adiabatic = false;
        in.Ue = 900;
        in.Te = 65;
        in.delta = 0.01;
        in.Cp = 1005;
        in.tauwGuess = 30;

        const GriffinWallFunction::InnerResult res =
            GriffinWallFunction::solveInner(in, c);

        require(res.converged == 1, "compressible shooting converged", failures);
        require(res.tauw > 0, "compressible tau_w positive", failures);
        require
        (
            res.Tr > in.Te && res.Tr > in.Tw,
            "cooled wall: T_r above T_e and T_w",
            failures
        );
        require
        (
            res.qw < 0,
            "cooled wall extracts heat from the fluid (q_w < 0, OF sign)",
            failures
        );

        const scalar qwExpected =
            c.s*res.tauw*in.Cp*(in.Tw - res.Tr)/in.Ue;
        require
        (
            close(res.qw, qwExpected, 1e-10),
            "q_w = s tau_w Cp (T_w - T_r)/U_e",
            failures
        );

        // Adiabatic variant: q_w = 0 identically and T_w = T_r.
        GriffinWallFunction::InnerInput inAd(in);
        inAd.adiabatic = true;
        const GriffinWallFunction::InnerResult resAd =
            GriffinWallFunction::solveInner(inAd, c);
        require(resAd.converged == 1, "adiabatic shooting converged", failures);
        require(resAd.qw == 0, "adiabatic q_w = 0", failures);
        require
        (
            close(resAd.Tw, resAd.Tr, 1e-12),
            "adiabatic T_w = T_r",
            failures
        );
    }

    // ---- Boundary-layer edge finder (PRF 2021, compressible form) --------
    {
        // Uniform-temperature exponential profile: U = Ue (1 - exp(-5 y/d)),
        // V = 0, h constant. Then h0ref = h + U_last^2/2, U_I = U_last, and
        // the edge is where U = threshold*U_last.
        const label n = 400;
        const scalar Ue = 800;
        const scalar d = 0.01;
        const scalar h0 = 2e5;
        const scalar T0 = 200;
        const scalar yMax = 3*d;

        scalarField dist(n), Utan(n), Vn(n), h(n), T(n);
        forAll(dist, i)
        {
            const scalar y = (i + 1)*yMax/n;
            dist[i] = y;
            Utan[i] = Ue*(1 - Foam::exp(-5*y/d));
            Vn[i] = 0;
            h[i] = h0;
            T[i] = T0 + 1e4*y;   // linear ramp: checks Te interpolation
        }

        const scalar ULast = Utan[n - 1];
        const scalar deltaExpected = -d/5*Foam::log(1 - 0.99*ULast/Ue);

        scalar delta = 0, UeOut = 0, TeOut = 0;
        const label status = GriffinWallFunction::findEdge
        (
            dist, Utan, Vn, h, T, 0.99, delta, UeOut, TeOut
        );

        require(status == 0, "edge finder found a crossing", failures);
        require
        (
            close(delta, deltaExpected, 5e-3),
            "edge finder delta matches analytic 0.99 crossing",
            failures
        );
        require
        (
            close(UeOut, 0.99*ULast, 5e-3),
            "edge finder Ue = 0.99 U_I at the crossing",
            failures
        );
        require
        (
            close(TeOut, T0 + 1e4*delta, 5e-3),
            "edge finder Te interpolated at delta",
            failures
        );

        // No-crossing profile that also exercises the arming logic:
        // constant U with h decreasing outward puts the stagnation-
        // enthalpy maximum at the first sample, so U/U_I starts at 1
        // (above the threshold: search not yet armed), drops below the
        // threshold after a few samples, and never crosses again ->
        // status 1 with the outermost sample as the edge estimate.
        const label nc = 50;
        scalarField distC(nc), UtanC(nc), VnC(nc), hC(nc), TC(nc);
        forAll(distC, i)
        {
            distC[i] = (i + 1)*1e-4;
            UtanC[i] = 400;
            VnC[i] = 0;
            hC[i] = 2e5 - 800*(i + 1);
            TC[i] = 200;
        }

        const label statusCut = GriffinWallFunction::findEdge
        (
            distC, UtanC, VnC, hC, TC, 0.99, delta, UeOut, TeOut
        );
        require
        (
            statusCut == 1 && close(delta, distC[nc - 1], 1e-12),
            "monotone-below-threshold profile: status 1, outermost sample",
            failures
        );
    }

    if (failures)
    {
        Info<< nl << failures << " test(s) FAILED" << nl << endl;
        return 1;
    }

    Info<< nl << "All tests passed" << nl << endl;
    return 0;
}

// ************************************************************************* //
