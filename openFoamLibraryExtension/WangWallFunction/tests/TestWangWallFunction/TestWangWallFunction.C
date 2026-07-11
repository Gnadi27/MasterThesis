/*---------------------------------------------------------------------------*\
  Standalone verification driver for the Wang AWF^Cmp per-face solver.

  Includes the library translation unit so the anonymous-namespace
  solveFace/FaceInput become visible, runs a fixed table of synthetic
  inputs (mirrored in tests/reference_wang_awf.py), and prints one
  machine-readable line per case for comparison against the independent
  Python implementation of paper Eqs. 30-45.
\*---------------------------------------------------------------------------*/

#include "../../WangNeumannWallFunction/WangNeumannWallFunction.C"

#include <cstdio>

using namespace Foam;

int main()
{
    WangNeumannWallFunction::Coeffs coeffs;   // paper defaults

    // Synthetic cases; keep in sync with reference_wang_awf.py
    for (label i = 0; i < 8; ++i)
    {
        FaceInput in;

        const bool thin = (i >= 6);

        in.kP = thin ? 0.01 : (0.5 + 0.37*i);
        in.yP = thin ? 5e-05 : (0.0004 + 0.0001*i);
        in.yn = 2*in.yP;
        in.Un = 50 + 30*i;
        in.Tn = 250 + 40*i;
        in.conv = ((i % 2 == 0) ? 1.0 : -1.0)*2.0e5*i;
        in.dpdx = ((i % 2 == 0) ? -1.0 : 1.0)*3.0e4*i;
        in.divE = ((i % 2 == 0) ? 1.0 : -1.0)*5.0e7*i;
        in.muv = 1.8e-05 + 2.0e-06*i;
        in.muw = in.muv*(1.0 + 0.15*((i % 3) - 1));
        in.rhov = 1.2 - 0.05*i;
        in.cp = 1004.5;
        in.Pr = 0.72;
        in.Tw = 300 + 10*i;
        in.adiabatic = (i % 4 == 3);

        const FaceOutput out = solveFace(in, coeffs);

        std::printf
        (
            "case %ld %.16e %.16e %.16e %.16e %.16e %.16e %.16e "
            "%.16e %.16e %.16e %.16e %.16e %.16e %.16e %d\n",
            static_cast<long>(i),
            out.A1, out.tauw, out.qw, out.Ath1,
            out.TvNew, out.TwNew,
            out.Pk, out.epsBar, out.epsP, out.dissAvg,
            out.yStar, out.yStarN, out.bMu,
            out.D, // D coefficient; C and Dth recomputable but print D
            out.thin ? 1 : 0
        );
    }

    return 0;
}

// ************************************************************************* //
