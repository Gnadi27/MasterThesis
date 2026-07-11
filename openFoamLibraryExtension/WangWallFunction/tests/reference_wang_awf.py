#!/usr/bin/env python3
"""Independent reference implementation of the Wang et al. (2023) AWF^Cmp
(paper Eqs. 30-45, thesis Appendix E thin-cell branch), written directly
from the printed equations, for cross-checking the OpenFOAM library's
per-face solver (tests/TestWangWallFunction).

Also runs internal consistency checks on its own solution:
  - profile continuity at y*_v (U, mu dU/dy*, T),
  - boundary conditions U2(y*_n) = Un, T2(y*_n) = Tn,
  - reduction to the Craft (2002) forced-convection formulas for
    b_mu = 0, D = C = Dth = 0,
  - q_w/Tn-Tw inversion round trip.

Usage:
    TestWangWallFunction > cpp_output.txt
    python3 reference_wang_awf.py cpp_output.txt
"""

import math
import sys

# Paper constants (must match Coeffs defaults)
Y_STAR_V = 10.8
CL = 2.55
CMU = 0.09
PRT = 0.9
ALPHA = CMU*CL           # 0.2295
Y_STAR_D = 2*CL          # 5.1
NPTS = 40                # Simpson intervals (matches nIntegrationPoints)


def simpson(a, b, n, f):
    """Composite Simpson matching the library's rule exactly."""
    if b <= a:
        return 0.0
    nn = max(2, n + (n % 2))
    h = (b - a)/nn
    s = f(a) + f(b)
    for i in range(1, nn):
        s += (4.0 if i % 2 else 2.0)*f(a + i*h)
    return s*h/3.0


def solve_case(kP, yP, yn, Un, Tn, conv, dpdx, divE,
               muv, muw, rhov, cp, Pr, Tw, adiabatic):
    """Solve one face; returns a dict of outputs."""
    kP = max(kP, 1e-15)
    sqrtk = math.sqrt(kP)
    scale = rhov*sqrtk/muv               # 1/m; y* = scale*y (WF26)
    ysP = scale*yP
    ysN = scale*yn

    alphat = ALPHA*Pr/PRT
    thin = ysN <= Y_STAR_V
    yv = ysN if thin else Y_STAR_V

    # Scaled coefficients: paper WF18/WF23 scaling (1/y*_n^2), thesis
    # hyper-CMAWF property state (mu_v, rho_v)
    denomK = rhov*rhov*kP
    D = muv*muv/(denomK*ysN*ysN)*conv
    C = muv*muv/denomK*dpdx
    Dth = muv*muv/(cp*denomK*ysN*ysN)*divE

    # Hyperbolic viscosity (WF28 / E.1)
    bmu = (muw - muv)/(muw*yv)

    beta = 1 - ALPHA*yv                  # (1 - alpha y*_v)

    def lnY(s):
        return math.log(max(1 + ALPHA*(s - yv), 1e-15))

    def lnYT(s):
        return math.log(max(1 + alphat*(s - yv), 1e-15))

    # ---------------- momentum ----------------
    if not thin:
        # N as printed in WF35 (D1=D2=D, C1=C2=C -> difference terms zero)
        def brJ(s):   # bracket of WF31/WF35
            return (ALPHA**2/3*s**3 - ALPHA*beta/2*s**2 + beta**2*s
                    - beta**3/ALPHA*lnY(s))

        N = (D/3*(1/ALPHA**3)*brJ(ysN)
             + C/ALPHA*(ysN - (1/ALPHA - yv)*lnY(ysN))
             + (D*yv**4/12 + C*yv**2/2
                - bmu*(D*yv**5/60 + C*yv**3/6))
             - (D/3*(1/ALPHA**3)*(ALPHA**2/3*yv**3 - ALPHA*beta/2*yv**2
                                  + beta**2*yv)
                + C/ALPHA*yv))

        A1 = (muv*Un - N)/(lnY(ysN)/ALPHA + yv - 0.5*bmu*yv**2)   # WF32
        A2 = A1                                                    # WF33
        # B2 as printed in WF34
        B2 = (D*yv**4/12 + C*yv**2/2 + A1*yv
              - bmu*(D*yv**5/60 + C*yv**3/6 + 0.5*A1*yv**2)
              - (D/3*(1/ALPHA**3)*(ALPHA**2/3*yv**3 - ALPHA*beta/2*yv**2
                                   + beta**2*yv)
                 + C/ALPHA*yv))
    else:
        # Thesis E.2
        poly = (D*ysN**4/12 + C*ysN**2/2
                - bmu*ysN**3*(D*ysN**2/60 + C/6))
        A1 = (muv*Un - poly)/(ysN - 0.5*bmu*ysN**2)
        A2 = A1
        B2 = 0.0

    # Conventional wall shear stress (papers print the on-fluid sign,
    # tau_wall = -scale*A1; the library stores +scale*A1, see design doc)
    tauw = scale*A1

    # Profiles (WF30/WF31)
    def U1(s):
        return (D*s**4/12 + C*s**2/2 + A1*s
                + bmu*(D*s**5/15 - D*yv*s**4/12 + C*s**3/3
                       + 0.5*(A1 - C*yv)*s**2 - A1*yv*s))/muv

    def U2(s):
        return ((D/3*(1/ALPHA**3)*(ALPHA**2/3*s**3 - ALPHA*beta/2*s**2
                                   + beta**2*s - beta**3/ALPHA*lnY(s)))
                + C/ALPHA*(s - (1/ALPHA - yv)*lnY(s))
                + A2/ALPHA*lnY(s) + B2)/muv

    def tau1(s):   # mu dU1/dy*
        return D*s**3/3 + C*s + A1

    def tau2(s):   # mu_v Y dU2/dy*
        return D*s**3/3 + C*s + A2

    # ---------------- k-equation quantities ----------------
    Pk = 0.0
    if not thin:
        def pk_integrand(s):
            Y = 1 + ALPHA*(s - yv)
            return ALPHA*(s - yv)*tau2(s)**2/(muv*Y*Y)
        Pk = scale/yn*simpson(yv, ysN, NPTS, pk_integrand)

    k32 = kP*sqrtk
    if ysN > Y_STAR_D:
        epsBar = (2*k32/Y_STAR_D + k32/CL*math.log(ysN/Y_STAR_D))/yn  # WF16
    else:
        yd = Y_STAR_D/scale
        epsBar = 2*(muv/rhov)*kP/yd**2                               # C4.8
    epsP = k32/(CL*yP)                                               # C4.9

    # ---------------- thermal ----------------
    betat = 1 - alphat*yv

    def Fth(s):    # Dth2/3 (1/alphat^3) [ ... ] group of WF38/WF41/WF44
        return (Dth/3*(1/alphat**3)
                * (alphat**2/3*s**3 - alphat*betat/2*s**2 + betat**2*s
                   - betat**3/alphat*lnYT(s)))

    def I1(b):     # (1/Cp) int_0^b (1+bmu(s-yv)) U1 (mu dU1/dy*) ds
        return simpson(0.0, b, NPTS,
                       lambda s: (1 + bmu*(s - yv))*U1(s)*tau1(s))/cp

    def I2(b):     # (1/Cp) int_{yv}^b (1/YT) U2 (mu_v Y dU2/dy*) ds
        return simpson(yv, b, NPTS,
                       lambda s: U2(s)*tau2(s)
                       / (1 + alphat*(s - yv)))/cp

    if not thin:
        BT = (Dth*yv**4/12 - Dth/60*bmu*yv**5 - I1(yv)
              - Fth(yv))                                             # WF41
        denomTh = lnYT(ysN)/alphat + yv - 0.5*bmu*yv**2

        if adiabatic:
            Ath1 = 0.0
            qw = 0.0
            TwNew = Tn - Pr/muv*(Fth(ysN) + BT - I2(ysN))            # WF42
        else:
            Ath1 = ((muv/Pr)*(Tn - Tw) - Fth(ysN) - BT
                    + I2(ysN))/denomTh                               # WF44
            qw = -(rhov*cp*sqrtk/muv)*Ath1                           # WF43
            TwNew = Tw

        # Tv = T1(y*_v) (continuity value)
        TvNew = (Pr/muv*(Dth*yv**4/12 + Ath1*yv
                         - bmu*(Dth*yv**5/60 + 0.5*Ath1*yv**2)
                         - I1(yv))
                 + TwNew)
    else:
        denomTh = ysN - 0.5*bmu*ysN**2
        if adiabatic:
            Ath1 = 0.0
            qw = 0.0
            TwNew = Tn - Pr/muv*(Dth*ysN**4/12 - bmu*Dth*ysN**5/60
                                 - I1(ysN))                          # E.3
        else:
            Ath1 = ((muv/Pr)*(Tn - Tw) - Dth*ysN**4/12
                    + bmu*Dth*ysN**5/60 + I1(ysN))/denomTh           # E.4
            qw = -(rhov*cp*sqrtk/muv)*Ath1
            TwNew = Tw
        TvNew = Tn

    dissAvg = Un*scale*tau2(ysN)/yn                                  # WF24

    return dict(A1=A1, tauw=tauw, qw=qw, Ath1=Ath1,
                TvNew=TvNew, TwNew=TwNew,
                Pk=Pk, epsBar=epsBar, epsP=epsP, dissAvg=dissAvg,
                yStar=ysP, yStarN=ysN, bMu=bmu, D=D, thin=int(thin),
                # extras for internal checks
                _U1=U1, _U2=U2, _tau1=tau1, _tau2=tau2,
                _yv=yv, _ysN=ysN, _Un=Un, _Tn=Tn, _muv=muv,
                _Fth=Fth, _I1=I1, _I2=I2, _BT=None if thin else BT,
                _bmu=bmu, _Dth=Dth, _Pr=Pr, _alphat=alphat,
                _adiabatic=adiabatic, _Tw=TwNew)


def cases():
    """Synthetic inputs; keep in sync with TestWangWallFunction.C."""
    for i in range(8):
        thin = i >= 6
        kP = 0.01 if thin else 0.5 + 0.37*i
        yP = 5e-05 if thin else 0.0004 + 0.0001*i
        yield dict(
            kP=kP, yP=yP, yn=2*yP,
            Un=50 + 30*i, Tn=250 + 40*i,
            conv=(1.0 if i % 2 == 0 else -1.0)*2.0e5*i,
            dpdx=(-1.0 if i % 2 == 0 else 1.0)*3.0e4*i,
            divE=(1.0 if i % 2 == 0 else -1.0)*5.0e7*i,
            muv=1.8e-05 + 2.0e-06*i,
            muw=(1.8e-05 + 2.0e-06*i)*(1.0 + 0.15*((i % 3) - 1)),
            rhov=1.2 - 0.05*i, cp=1004.5, Pr=0.72,
            Tw=300 + 10*i, adiabatic=(i % 4 == 3))


COLUMNS = ["A1", "tauw", "qw", "Ath1", "TvNew", "TwNew",
           "Pk", "epsBar", "epsP", "dissAvg", "yStar", "yStarN",
           "bMu", "D", "thin"]


def internal_checks(sol, tol=1e-8):
    """Consistency checks on the Python solution itself."""
    fails = []
    yv, ysN = sol["_yv"], sol["_ysN"]
    U1, U2 = sol["_U1"], sol["_U2"]
    muv = sol["_muv"]

    if not sol["thin"]:
        # Continuity of U and of the stress first-integral at y*_v
        du = abs(U1(yv) - U2(yv))
        if du > tol*max(1.0, abs(U1(yv))):
            fails.append(f"U continuity at yv: {du:.3e}")
        # mu dU1/dy*(yv) = tau1(yv) vs mu_v Y dU2/dy*(yv) = tau2(yv), Y=1
        dt = abs(sol["_tau1"](yv) - sol["_tau2"](yv))
        if dt > tol*max(1.0, abs(sol["_tau1"](yv))):
            fails.append(f"stress continuity at yv: {dt:.3e}")
        # U2(y*_n) = Un
        dn = abs(U2(ysN) - sol["_Un"])
        if dn > tol*max(1.0, abs(sol["_Un"])):
            fails.append(f"U2(ysN) != Un: {dn:.3e}")

        # T2(y*_n) = Tn: rebuild T2 from WF38 and check
        Pr, alphat = sol["_Pr"], sol["_alphat"]
        Fth, I2, BT = sol["_Fth"], sol["_I2"], sol["_BT"]
        Ath1, Tw = sol["Ath1"], sol["_Tw"]
        bmu = sol["_bmu"]

        def lnYT(s):
            return math.log(1 + alphat*(s - yv))

        def T2(s):
            Bth2 = Ath1*yv - 0.5*bmu*Ath1*yv**2 + BT + muv/Pr*Tw
            return Pr/muv*(Fth(s) + Ath1/alphat*lnYT(s) + Bth2
                           - I2(s))

        dT = abs(T2(ysN) - sol["_Tn"])
        if dT > 1e-6*max(1.0, abs(sol["_Tn"])):
            fails.append(f"T2(ysN) != Tn: {dT:.3e}")

        # T continuity: T1(yv) = T2(yv)
        Dth = sol["_Dth"]
        I1 = sol["_I1"]
        T1v = (Pr/muv*(Dth*yv**4/12 + Ath1*yv
                       - bmu*(Dth*yv**5/60 + 0.5*Ath1*yv**2)
                       - I1(yv)) + Tw)
        dTv = abs(T1v - T2(yv))
        if dTv > 1e-6*max(1.0, abs(T1v)):
            fails.append(f"T continuity at yv: {dTv:.3e}")
    else:
        dn = abs(U1(ysN) - sol["_Un"])
        if dn > tol*max(1.0, abs(sol["_Un"])):
            fails.append(f"thin U1(ysN) != Un: {dn:.3e}")

    return fails


def craft_limit_check():
    """b_mu = 0, D = C = Dth = 0 must reduce to Craft (2002) A.11/B.5:
    A1 = mu_v Un / (ln(Yn)/alpha + y*_v),
    Ath1 = (Tn - Tw) (mu_v/Pr) / (ln(YTn)/alphat + y*_v)."""
    kP, yP = 1.0, 5e-4
    yn, Un, Tn, Tw = 1e-3, 100.0, 280.0, 320.0
    muv = 1.8e-5
    rhov, cp, Pr = 1.0, 1004.5, 0.72
    sol = solve_case(kP, yP, yn, Un, Tn, 0, 0, 0,
                     muv, muv, rhov, cp, Pr, Tw, False)
    scale = rhov*math.sqrt(kP)/muv
    ysN = scale*yn
    yv = Y_STAR_V
    alphat = ALPHA*Pr/PRT
    A1_craft = muv*Un/(math.log(1 + ALPHA*(ysN - yv))/ALPHA + yv)
    Ath1_craft = ((muv/Pr)*(Tn - Tw)
                  / (math.log(1 + alphat*(ysN - yv))/alphat + yv))
    e1 = abs(sol["A1"] - A1_craft)/abs(A1_craft)
    # Dissipation integrals are NOT in Craft's thermal AWF; with the
    # velocity nonzero they contribute here, so compare only the no-flow
    # limit for the thermal constant:
    sol0 = solve_case(kP, yP, yn, 0.0, Tn, 0, 0, 0,
                      muv, muv, rhov, cp, Pr, Tw, False)
    e2 = abs(sol0["Ath1"] - Ath1_craft)/abs(Ath1_craft)
    ok = e1 < 1e-12 and e2 < 1e-12
    print(f"Craft-limit check: A1 rel err {e1:.3e}, "
          f"Ath1 (U=0) rel err {e2:.3e} -> {'OK' if ok else 'FAIL'}")
    return ok


def qw_roundtrip_check():
    """Prescribing the q_w obtained from the isothermal solve as a heat
    flux must reproduce Tw through WF42 (adiabatic-form with Ath1)."""
    ins = dict(kP=1.3, yP=6e-4, yn=1.2e-3, Un=180.0, Tn=400.0,
               conv=1e5, dpdx=-2e4, divE=3e7,
               muv=2.1e-5, muw=2.4e-5, rhov=0.9, cp=1004.5, Pr=0.72,
               Tw=350.0, adiabatic=False)
    sol = solve_case(**ins)
    # WF42: Tw = Tn - Pr/muv (Fth(ysN) + lnYTn/alphat*Ath1 + Ath1 yv
    #             - 1/2 bmu Ath1 yv^2 + BT - I2n)   [Ath2 = Ath1]
    yv, ysN = sol["_yv"], sol["_ysN"]
    alphat = sol["_alphat"]
    lnYTn = math.log(1 + alphat*(ysN - yv))
    Tw_back = (ins["Tn"] - sol["_Pr"]/ins["muv"]
               * (sol["_Fth"](ysN)
                  + (lnYTn/alphat + yv - 0.5*sol["_bmu"]*yv**2)
                  * sol["Ath1"]
                  + sol["_BT"] - sol["_I2"](ysN)))
    err = abs(Tw_back - ins["Tw"])
    ok = err < 1e-8
    print(f"q_w/Tw round-trip check: |Tw_back - Tw| = {err:.3e} K "
          f"-> {'OK' if ok else 'FAIL'}")
    return ok


def main():
    all_ok = True

    # Internal consistency of the reference itself
    for i, ins in enumerate(cases()):
        sol = solve_case(**ins)
        fails = internal_checks(sol)
        if fails:
            all_ok = False
            print(f"case {i}: INTERNAL CHECK FAILURES: {fails}")
        else:
            print(f"case {i}: internal checks OK "
                  f"(thin={sol['thin']}, y*n={sol['yStarN']:.3f})")

    all_ok &= craft_limit_check()
    all_ok &= qw_roundtrip_check()

    # Cross-check against the C++ output, if provided
    if len(sys.argv) > 1:
        ref = {i: solve_case(**ins) for i, ins in enumerate(cases())}
        ncomp = nfail = 0
        with open(sys.argv[1]) as f:
            for line in f:
                parts = line.split()
                if not parts or parts[0] != "case":
                    continue
                i = int(parts[1])
                vals = [float(x) for x in parts[2:]]
                for col, cpp in zip(COLUMNS, vals):
                    py = ref[i][col]
                    denom = max(1.0, abs(py))
                    err = abs(cpp - py)/denom
                    ncomp += 1
                    if err > 1e-9:
                        nfail += 1
                        all_ok = False
                        print(f"case {i} {col}: cpp={cpp:.12e} "
                              f"py={py:.12e} rel={err:.3e}  MISMATCH")
        print(f"C++ cross-check: {ncomp} values compared, "
              f"{nfail} mismatches")

    print("ALL OK" if all_ok else "FAILURES PRESENT")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
