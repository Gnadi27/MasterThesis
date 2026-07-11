/*---------------------------------------------------------------------------*\
  Neumann (prescribed-gradient) energy coupling for the Wang et al. (2023)
  compressible Analytical Wall Function (AWF^Cmp).

  Boundary condition on T that imposes the AWF wall heat flux q_w (paper
  Eqs. WF43/WF44) directly by prescribing the wall-normal temperature
  gradient, so that the discrete energy flux through the wall face equals
  the model q_w with plain non-negative diffusivities, independent of the
  sign of (T_1 - T_w). Delivery mechanics are identical to
  GaoNeumannWallHeatFluxTemperature (see ../GaoWallFunctions).

  The patch T value under this boundary condition is the tangent-
  extrapolated face temperature, NOT the physical wall temperature. Tw and
  the wall transport properties are supplied through the wall-function
  coefficients (Tw, wallTransport) shared by all WangNeumann* boundary
  conditions on the patch. The wall's alphat patch should be a trivial
  fixedValue 0.

  Usage (isothermal wall):
      wall
      {
          type            WangNeumannWallHeatFluxTemperature;
          Tw              300;
          wallTransport
          {
              mode        sutherland;
              As          1.458e-06;
              Ts          110.4;
              Pr          0.7;    // paper value; omitted -> modified Eucken
          }
          value           uniform 300;
      }
\*---------------------------------------------------------------------------*/

#ifndef WangNeumannWallHeatFluxTemperatureFvPatchScalarField_H
#define WangNeumannWallHeatFluxTemperatureFvPatchScalarField_H

#include "WangNeumannWallFunction.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class WangNeumannWallHeatFluxTemperatureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    WangNeumannWallFunction::Coeffs coeffs_;
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
        const WangNeumannWallFunction::PatchData& data
    );

public:
    TypeName("WangNeumannWallHeatFluxTemperature");

    WangNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    WangNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    WangNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const WangNeumannWallHeatFluxTemperatureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    WangNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const WangNeumannWallHeatFluxTemperatureFvPatchScalarField&
    );

    WangNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const WangNeumannWallHeatFluxTemperatureFvPatchScalarField&,
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
    //  WangNeumannWallFunction::isAdiabatic to resolve thermalMode auto.
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
