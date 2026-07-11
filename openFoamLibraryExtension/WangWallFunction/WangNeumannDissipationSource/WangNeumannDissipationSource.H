/*---------------------------------------------------------------------------*\
  CMAWF-D near-wall-cell viscous-dissipation correction (paper Eq. WF24 /
  thesis Eq. T8.64) for the Wang et al. (2023) AWF^Cmp Neumann family.

  In the CMAWF-D form of the model, the wall-normal viscous dissipation
  term of the *discretized* total-energy equation over near-wall cells is
  replaced by a cell-averaged value evaluated from the analytical velocity
  profile:

      dissAvg = (1/y_n) [ (mu+mu_t) U dU/dy ]_0^{y_n} = Un tau(y_n)/y_n

  (exact, since U(0) = 0). Wang implemented this by modifying the solver;
  here the solver's own version of the term is embedded in its face-based
  viscous fluxes and cannot be removed, so this fvOption adds the
  difference

      S = dissAvg - dissLocal        [W/m^3]

  in the near-wall cells of the selected wall patches, where dissLocal is
  the equivalent discrete estimate built from the resolved fields
  ([mu_eff U_t dU_t/dn]_northface / y_n). This is the roughest paper-to-
  OpenFOAM mapping in the library -- see DESIGN_AND_DEVIATIONS.md. Omit
  the fvOption to run the plain CMAWF (dissipation only in the analytical
  temperature profile).

  Usage (constant/fvOptions):

      WangDissipation
      {
          type            WangNeumannDissipationSource;
          patches         (wall);
          // Coefficients identical to the WangNeumann* BCs on the patch:
          thermalMode     isothermal;
          Tw              300;
          wallTransport { mode sutherland; As 1.458e-06; Ts 110.4; }
      }
\*---------------------------------------------------------------------------*/

#ifndef WangNeumannDissipationSource_H
#define WangNeumannDissipationSource_H

#include "WangNeumannWallFunction.H"
#include "fvOption.H"
#include "wordRes.H"

namespace Foam
{
namespace fv
{

class WangNeumannDissipationSource
:
    public fv::option
{
    //- Wall-function coefficients (base-class coeffs_ is the dictionary)
    WangNeumannWallFunction::Coeffs wfCoeffs_;

    //- Selected wall patches
    wordRes patchNames_;

public:

    TypeName("WangNeumannDissipationSource");

    WangNeumannDissipationSource
    (
        const word& sourceName,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- No copy construct
    WangNeumannDissipationSource
    (
        const WangNeumannDissipationSource&
    ) = delete;

    //- No copy assignment
    void operator=(const WangNeumannDissipationSource&) = delete;

    virtual ~WangNeumannDissipationSource() = default;

    //- Add explicit contribution to the compressible energy equation
    virtual void addSup
    (
        const volScalarField& rho,
        fvMatrix<scalar>& eqn,
        const label fieldi
    );

    virtual bool read(const dictionary& dict);
};

} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
