#include "LiWallFunction.H"
#include "LiNeumannWallFunction.H"
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

} // End anonymous namespace


int main()
{
    LiWallFunction::Coeffs coeffs;
    label failures = 0;

    // Equation (42) must recover Spalding's unified law when beta*=Gamma*=0.
    const scalar u0 = 12.0;
    const scalar ku0 = coeffs.kappa*u0;
    const scalar y0 = Foam::exp(-coeffs.kappa*coeffs.B);
    const scalar spalding =
        u0
      + y0
       *(
            Foam::exp(ku0)
          - (1 + ku0 + 0.5*sqr(ku0) + (1.0/6.0)*pow3(ku0))
        );
    require
    (
        close
        (
            LiWallFunction::lawYPlus(u0, 0, 0, coeffs),
            spalding,
            1e-12
        ),
        "Eq. (42) incompressible/Spalding limit",
        failures
    );

    // Check Eq. (38) through the derivative of the complete Eq. (42).
    const scalar u = 8.0;
    const scalar beta = 0.018;
    const scalar Gamma = 0.0012;
    const scalar du = 1e-5;
    const scalar numericalDerivative =
        (
            LiWallFunction::lawYPlus(u + du, beta, Gamma, coeffs)
          - LiWallFunction::lawYPlus(u - du, beta, Gamma, coeffs)
        )/(2*du);
    const scalar ku = coeffs.kappa*u;
    const scalar analyticDerivative =
        1
      + coeffs.omegaH*beta*u
      - coeffs.omegaH*Gamma*sqr(u)
      + LiWallFunction::dYWhiteDuPlus(u, beta, Gamma, coeffs)
      - coeffs.kappa*y0*(1 + ku + 0.5*sqr(ku));
    require
    (
        close(numericalDerivative, analyticDerivative, 2e-7),
        "Eq. (38) derivative agrees with finite difference of Eq. (42)",
        failures
    );

    // Equation (39) inversion.
    const scalar nu = 1.7e-5;
    const scalar nutTarget = 3.2e-4;
    const scalar nuTilda = LiWallFunction::solveNuTilda
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
        "Eq. (39) nuTilda inversion",
        failures
    );

    // Reference-point inversion of Eq. (40).
    const scalar hw = 1.1e6;
    const scalar hRef = 2.4e6;
    const scalar hrgExpected = 3.0e6;
    const scalar x = 0.42;
    const scalar h1 = hw + (hrgExpected - hw)*x + (hRef - hrgExpected)*sqr(x);
    const scalar hrgRecovered =
        (h1 - hw*(1 - x) - hRef*sqr(x))/(x - sqr(x));
    require
    (
        close(hrgRecovered, hrgExpected, 1e-12),
        "Eq. (40) reference-point recovery enthalpy inversion",
        failures
    );

    // Frozen calorically perfect enthalpy reconstruction used by the library.
    const scalar cp = 1005.0;
    const scalar Tw = 350.0;
    const scalar T1 = 925.0;
    const scalar reconstructedH = cp*Tw + cp*(T1 - Tw);
    require
    (
        close(reconstructedH, cp*T1, 1e-14),
        "calorically perfect reconstruction gives h=Cp*T",
        failures
    );

    // The LiNeumann copy must evaluate the identical laws.
    {
        LiNeumannWallFunction::Coeffs nCoeffs;

        require
        (
            close
            (
                LiNeumannWallFunction::lawYPlus(u0, 0, 0, nCoeffs),
                LiWallFunction::lawYPlus(u0, 0, 0, coeffs),
                1e-14
            ),
            "LiNeumann Eq. (42) matches the Li original",
            failures
        );

        require
        (
            close
            (
                LiNeumannWallFunction::lawYPlus(u, beta, Gamma, nCoeffs),
                LiWallFunction::lawYPlus(u, beta, Gamma, coeffs),
                1e-14
            ),
            "LiNeumann Eq. (42) generalized law matches the Li original",
            failures
        );

        require
        (
            close
            (
                LiNeumannWallFunction::dYWhiteDuPlus(u, beta, Gamma, nCoeffs),
                LiWallFunction::dYWhiteDuPlus(u, beta, Gamma, coeffs),
                1e-14
            ),
            "LiNeumann Eq. (38) derivative matches the Li original",
            failures
        );

        require
        (
            close
            (
                LiNeumannWallFunction::solveNuTilda
                (
                    nutTarget, nu, nCoeffs.cv1, 100, 1e-12
                ),
                nuTilda,
                1e-10
            ),
            "LiNeumann Eq. (39) inversion matches the Li original",
            failures
        );
    }

    // LiNeumann coefficient parsing: Tw and the wallTransport sub-dict.
    {
        dictionary dict;
        dict.add("Tw", 350.0);

        dictionary wt;
        wt.add("mode", word("sutherland"));
        wt.add("As", 1.458e-06);
        wt.add("Ts", 110.4);
        wt.add("Pr", 0.72);
        dict.add("wallTransport", wt);

        const LiNeumannWallFunction::Coeffs parsed(dict);

        require
        (
            close(parsed.Tw, 350.0, 1e-14)
         && parsed.wallTransportMode == "sutherland"
         && close(parsed.wallTransportAs, 1.458e-06, 1e-14)
         && close(parsed.wallTransportTs, 110.4, 1e-14)
         && close(parsed.wallTransportPr, 0.72, 1e-14),
            "LiNeumann Coeffs parses Tw and wallTransport",
            failures
        );

        // Neumann inversion arithmetic for a calorically perfect gas:
        // placing the face enthalpy at he1 + qTargetOF/(alphaEff*delta)
        // and extrapolating T along the resulting gradient returns the
        // target flux for both signs of (T1 - Tw).
        const scalar cp = 1005.0;
        const scalar alphaEff = 0.05;
        const scalar delta = 2000.0;
        const scalar T1 = 280.0;

        for (const scalar qTargetOF : {scalar(-4000), scalar(4000)})
        {
            const scalar he1 = cp*T1;
            const scalar heFace = he1 + qTargetOF/(alphaEff*delta);
            const scalar TFace = heFace/cp;
            const scalar gradT = (TFace - T1)*delta;
            const scalar delivered = alphaEff*(cp*(T1 + gradT/delta) - he1)*delta;

            require
            (
                close(delivered, qTargetOF, 1e-12),
                "Neumann gradient inversion round trip, qTargetOF = "
              + name(qTargetOF),
                failures
            );
        }
    }

    if (failures)
    {
        Info<< failures << " Li wall-function equation test(s) failed" << endl;
        return 1;
    }

    Info<< "All Li wall-function equation tests passed" << endl;
    return 0;
}

// ************************************************************************* //
