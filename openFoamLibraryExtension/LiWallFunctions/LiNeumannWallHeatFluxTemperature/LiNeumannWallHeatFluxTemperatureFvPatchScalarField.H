/*---------------------------------------------------------------------------*\
  Neumann (prescribed-gradient) energy coupling for the Li et al.
  frozen-composition wall function.

  Boundary condition on T that imposes the Li wall heat flux directly by
  prescribing the wall-normal temperature gradient, so that the discrete
  energy flux through the wall face equals the model q_w with plain
  non-negative diffusivities, independent of the sign of (T_1 - T_w).
  This makes the T_1 < T_w regime (near-adiabatic walls at coarse y+)
  representable, which the alphat-inversion coupling of
  compressible::alphatLiWallFunction structurally cannot deliver.

  The patch T value under this boundary condition is the tangent-
  extrapolated face temperature, NOT the physical wall temperature. Tw and
  the wall transport properties are supplied through the wall-function
  coefficients (Tw, wallTransport) shared by all Li Neumann boundary
  conditions on the patch. The wall's alphat patch should be a trivial
  fixedValue 0.

  Usage (isothermal wall):
      wall
      {
          type            LiNeumannWallHeatFluxTemperature;
          Tw              300;
          wallTransport
          {
              mode        sutherland;
              As          1.458e-06;
              Ts          110.4;
              Pr          0.72;   // optional; omitted -> modified Eucken
          }
          value           uniform 300;
      }
\*---------------------------------------------------------------------------*/

#ifndef LiNeumannWallHeatFluxTemperatureFvPatchScalarField_H
#define LiNeumannWallHeatFluxTemperatureFvPatchScalarField_H

#include "LiNeumannWallFunction.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class LiNeumannWallHeatFluxTemperatureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    LiNeumannWallFunction::Coeffs coeffs_;
    word thermalFluxTarget_;
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
        const word& selectedTarget,
        const scalarField& alphaEff,
        const scalarField& heGradient,
        const scalarField& TFace,
        const LiNeumannWallFunction::PatchData& data
    );

public:
    TypeName("LiNeumannWallHeatFluxTemperature");

    LiNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    LiNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    LiNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const LiNeumannWallHeatFluxTemperatureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    LiNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const LiNeumannWallHeatFluxTemperatureFvPatchScalarField&
    );

    LiNeumannWallHeatFluxTemperatureFvPatchScalarField
    (
        const LiNeumannWallHeatFluxTemperatureFvPatchScalarField&,
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
    //  LiNeumannWallFunction::isAdiabatic to resolve thermalMode auto.
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
