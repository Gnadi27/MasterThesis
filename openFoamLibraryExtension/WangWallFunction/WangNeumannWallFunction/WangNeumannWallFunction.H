/*---------------------------------------------------------------------------*\
  Wang et al. (2023) compressible Analytical Wall Function (AWF^Cmp) for
  OpenFOAM v2412, with Neumann (prescribed-gradient) energy coupling.

  Reference: X. Wang, J. Chen, T. Craft, H. Iacovides, "Extension of
  Analytical Wall Functions to Supersonic and Hypersonic Flows", Flow,
  Turbulence and Combustion 111:1-34 (2023) -- referred to as "paper";
  equation numbers WFxx below refer to the paper, Txx.yy to Wang's PhD
  thesis (Manchester, 2018).

  This namespace hosts the shared per-patch evaluator used by all
  WangNeumann* boundary conditions and the WangNeumannDissipationSource
  fvOption. For a given wall patch it solves the analytically integrated
  near-wall momentum and total-energy equations (paper Appendix, Eqs.
  WF30-WF45; hyperbolic molecular-viscosity variation across the
  zero-turbulent-viscosity sublayer) and returns:

    - wall shear stress tau_w                  (T8.12)
    - wall heat flux q_w (isothermal wall)     (WF43/WF44)
      or wall temperature T_w (adiabatic wall) (WF42)
    - cell-averaged turbulence-energy production Pk (WF13)
      and dissipation epsBar (WF16), nodal epsP (Craft 2002 Eq. 4.9)
    - cell-averaged wall-normal viscous dissipation of the discretized
      total-energy equation, dissAvg (paper Eq. WF24 / thesis Eq. T8.64;
      evaluated exactly as Un*tau(yn)/yn since U(0)=0)

  As in the Gao/Li Neumann families, the patch T value drifts along the
  prescribed gradient under the Neumann coupling and cannot carry the
  physical wall temperature: Tw and the wall transport properties are
  supplied through the coefficients (Tw, wallTransport).
\*---------------------------------------------------------------------------*/

#ifndef WangNeumannWallFunction_H
#define WangNeumannWallFunction_H

#include "fvPatch.H"
#include "labelField.H"
#include "scalarField.H"
#include "turbulentFluidThermoModel.H"
#include "word.H"

namespace Foam
{
class dictionary;
class Ostream;

namespace WangNeumannWallFunction
{

struct Coeffs
{
    // Model constants (defaults faithful to the paper)
    scalar yStarV;          // sublayer edge, 10.8 (paper Sect. 3.2)
    scalar cl;              // length-scale constant, 2.55
    scalar Cmu;             // 0.09; alpha = Cmu*cl = 0.2295
    scalar Prt;             // turbulent Prandtl number, 0.9

    // Numerics (not physics)
    label nIntegrationPoints;   // Simpson intervals per region (even)
    label TvIterations;         // fixed-point iterations for Tv (and
                                // adiabatic Tw)
    scalar TvRelaxation;        // under-relaxation of the Tv update

    // hyperbolic (paper Eq. WF28) | constant (b_mu = 0, CMAWF limit)
    word viscosityVariation;

    // auto | isothermal | adiabatic
    word thermalMode;

    // Physical wall temperature (isothermal walls); negative means unset
    scalar Tw;

    // Wall/sublayer transport evaluation at (p_wall, T):
    //   fromPatch  - only valid for adiabatic walls (patch T is physical)
    //   sutherland - mu = As*sqrt(T)/(1 + Ts/T); Pr constant if given,
    //                otherwise kappa from modified Eucken
    //   constant   - fixed mu and kappa
    word wallTransportMode;
    scalar wallTransportAs;
    scalar wallTransportTs;
    scalar wallTransportPr;
    scalar wallTransportMu;
    scalar wallTransportKappa;

    Coeffs();
    explicit Coeffs(const dictionary& dict);

    scalar alpha() const { return Cmu*cl; }
    scalar yStarD() const { return 2*cl; }

    void writeEntries(Ostream& os) const;
};


struct PatchData
{
    // Momentum
    scalarField tauw;       // wall shear stress [Pa], positive downstream
    scalarField A1;         // sublayer first-integral constant (WF32/E.2)
    scalarField nutShear;   // kinematic nut reproducing tauw through
                            // OpenFOAM's discrete wall-stress path
    scalarField yStar;      // nodal y* (paper Eq. WF26)
    scalarField yStarN;     // y* at the north face of the near-wall cell
    scalarField uTau;       // sqrt(|tauw|/rho_v), diagnostic only

    // Thermal
    scalarField qw;         // wall heat flux [W/m2], paper sign: positive
                            // from the wall into the fluid
    scalarField Tw;         // physical wall temperature (dictionary value
                            // or WF42 recovery value for adiabatic walls)
    scalarField Tv;         // temperature at the sublayer edge (T1(y*v))
    scalarField Ath1;       // thermal first-integral constant

    // Turbulence quantities for the near-wall cell
    scalarField Pk;         // cell-averaged k production [kg/(m s^3)] (WF13)
    scalarField epsBar;     // cell-averaged dissipation [m^2/s^3] (WF16)
    scalarField epsP;       // nodal dissipation [m^2/s^3] (Craft Eq. 4.9)

    // CMAWF-D cell-averaged wall-normal viscous dissipation [W/m^3]
    scalarField dissAvg;

    // Discrete estimate of the same term from the resolved fields:
    // [mu_eff U_t dU_t/dn]_northface / y_n (the wall-face contribution
    // vanishes since U = 0 there). The dissipation fvOption adds the
    // difference dissAvg - dissLocal because the solver's own version of
    // the term is embedded in its face-based viscous fluxes and cannot be
    // removed directly.
    scalarField dissLocal;

    // States and coefficients (diagnostics)
    scalarField muw;        // molecular viscosity at Tw
    scalarField muv;        // molecular viscosity at Tv
    scalarField rhov;       // density at (p_wall, Tv)
    scalarField cpw;        // cp at (p_wall, Tw)
    scalarField Pr;         // molecular Prandtl number used
    scalarField bMu;        // hyperbolic viscosity coefficient (WF28)
    scalarField Dcoef;      // scaled momentum convection (WF18)
    scalarField Ccoef;      // scaled pressure gradient (WF18)
    scalarField DthCoef;    // scaled energy convection (WF23)
    scalarField Un;         // wall-parallel velocity at the north face
    scalarField Tn;         // temperature at the north face
    scalarField yn;         // wall-normal extent of the near-wall cell
    labelField thinCell;    // 1 where y*n < y*v (thesis Appendix E branch)

    explicit PatchData(const label size);
};


//- Resolve the thermal mode against the actual T boundary condition,
//  recognising WangNeumannWallHeatFluxTemperature.
bool isAdiabatic
(
    const Coeffs& coeffs,
    const fvPatchScalarField& Tpatch
);


//- Evaluate the Wang AWF^Cmp on one wall patch
PatchData calcPatch
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const Coeffs& coeffs
);

} // End namespace WangNeumannWallFunction
} // End namespace Foam

#endif

// ************************************************************************* //
