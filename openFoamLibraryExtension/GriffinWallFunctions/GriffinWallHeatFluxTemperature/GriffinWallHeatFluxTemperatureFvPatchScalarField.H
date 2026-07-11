/*---------------------------------------------------------------------------*\
  Neumann (prescribed-gradient) energy coupling for the Griffin wall
  function.

  Boundary condition on T that imposes the Griffin wall heat flux
      q_w = s tau_w Cp (T_w - T_r)/U_e
  directly by prescribing the wall-normal temperature gradient, so that
  the discrete energy flux through the wall face equals the model q_w
  with plain non-negative diffusivities, independent of the sign of
  (T_1 - T_w). This mirrors the paper's direct-flux coupling and avoids
  the alphat-inversion pathology for near-adiabatic walls at coarse y+
  (LiWallFunctions/FINDINGS_M6Tw076_energy_coupling.md).

  Sign convention: Griffin's q_w (via the B_q definition) is positive for
  heat entering the fluid, which coincides with OpenFOAM's outward-normal
  boundary-flux sign; qTargetOF = q_w with no flip.

  The patch T value under this boundary condition is the tangent-
  extrapolated face temperature, NOT the physical wall temperature. Tw is
  supplied through the wall-function coefficients shared by all Griffin
  boundary conditions on the patch; the viscosity law is likewise
  dictionary-supplied. The wall's alphat patch should be a trivial
  fixedValue 0.

  Usage (isothermal wall):
      wall
      {
          type            GriffinWallHeatFluxTemperature;
          Tw              300;
          viscosity
          {
              mode        sutherland;
              As          1.458e-06;
              Ts          110.4;
          }
          value           uniform 300;
      }
\*---------------------------------------------------------------------------*/

#ifndef GriffinWallHeatFluxTemperatureFvPatchScalarField_H
#define GriffinWallHeatFluxTemperatureFvPatchScalarField_H

#include "GriffinWallFunction.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class GriffinWallHeatFluxTemperatureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    GriffinWallFunction::Coeffs coeffs_;
    scalar gradientRelaxation_;
    scalar TFaceMin_;
    scalar TFaceMax_;
    bool diagnostics_;
    bool diagnosticLog_;
    bool diagnosticWrite_;
    label diagnosticWriteInterval_;
    label diagnosticWritePrecision_;
    label diagnosticTimeIndex_;

    void checkType() const;
    void checkEntries() const;

    void writeDiagnostics
    (
        const scalarField& alphaEff,
        const scalarField& heGradient,
        const scalarField& TFace,
        const GriffinWallFunction::PatchData& data
    );

public:
    TypeName("GriffinWallHeatFluxTemperature");

    GriffinWallHeatFluxTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    GriffinWallHeatFluxTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    GriffinWallHeatFluxTemperatureFvPatchScalarField
    (
        const GriffinWallHeatFluxTemperatureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    GriffinWallHeatFluxTemperatureFvPatchScalarField
    (
        const GriffinWallHeatFluxTemperatureFvPatchScalarField&
    );

    GriffinWallHeatFluxTemperatureFvPatchScalarField
    (
        const GriffinWallHeatFluxTemperatureFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual tmp<fvPatchField<scalar>> clone() const
    {
        return fvPatchField<scalar>::Clone(*this);
    }

    virtual tmp<fvPatchField<scalar>> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const
    {
        return fvPatchField<scalar>::Clone(*this, iF);
    }

    //- True when this boundary condition was set to adiabatic; used by
    //  GriffinWallFunction::isAdiabatic to resolve thermalMode auto.
    bool adiabaticIntent() const
    {
        return coeffs_.thermalMode == "adiabatic";
    }

    virtual void updateCoeffs();
    virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
