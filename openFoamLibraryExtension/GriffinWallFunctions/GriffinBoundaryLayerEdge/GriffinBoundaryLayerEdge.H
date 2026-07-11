/*---------------------------------------------------------------------------*\
  Boundary-layer edge quantities for the Griffin wall function.

  Implements the local-reconstruction method of Griffin, Fu & Moin,
  "General method for determining the boundary layer thickness in
  nonequilibrium flows", PRF 6 (2021), in its compressible form
  (Appendix D): along a wall-normal profile the reconstructed inviscid
  velocity is
      U_I^2[y] = 2 (h_0,ref - h[y]) - V^2[y],
      h_0[y] = h[y] + |U|^2/2,   h_0,ref = max_y h_0[y],
  and the boundary-layer edge is the first wall-outward crossing of
  U/U_I = threshold (0.99). Stagnation enthalpy is conserved across
  shocks, so the reconstruction remains valid in shock-containing flows.

  The wall-normal profile is gathered by a cell walk from each wall face:
  starting at the face cell, repeatedly step to the neighbour cell with
  the largest inward-normal advancement. The walk cannot cross processor
  boundaries; decompose wall-tangentially, or use edgeMode prescribed.

  Results are cached per patch and shared by all Griffin boundary
  conditions; recomputation is throttled by edgeUpdateInterval.
\*---------------------------------------------------------------------------*/

#ifndef GriffinBoundaryLayerEdge_H
#define GriffinBoundaryLayerEdge_H

#include "GriffinWallFunction.H"
#include "Map.H"
#include "regIOobject.H"

namespace Foam
{
class fvMesh;

namespace GriffinWallFunction
{

// Pure-array edge finder on one wall-normal profile ordered from the wall
// outward. Returns 0 when a threshold crossing was found (delta/Ue/Te
// linearly interpolated at the crossing), 1 when no crossing exists (edge
// set to the outermost sample). If the ratio starts above the threshold
// (spurious near-wall values, e.g. strongly heated walls), the search is
// armed only once the ratio has dropped below the threshold, following
// the paper's advice to search from the wall outward.
label findEdge
(
    const UList<scalar>& dist,   // wall-normal distance
    const UList<scalar>& Utan,   // wall-tangential speed
    const UList<scalar>& Vn,     // wall-normal velocity component
    const UList<scalar>& h,      // sensible enthalpy
    const UList<scalar>& T,      // temperature
    const scalar threshold,
    scalar& delta,
    scalar& Ue,
    scalar& Te
);

} // End namespace GriffinWallFunction


class GriffinBoundaryLayerEdge
:
    public regIOobject
{
public:

    struct PatchEdgeData
    {
        scalarField delta;
        scalarField Ue;
        scalarField Te;
        labelField status;   // findEdge status; 2 = prescribed
        label lastTimeIndex = -1;
    };

private:

    const fvMesh& mesh_;
    Map<PatchEdgeData> patchData_;

    void reconstruct
    (
        const compressible::turbulenceModel& turbModel,
        const fvPatch& patch,
        const GriffinWallFunction::Coeffs& coeffs,
        PatchEdgeData& pd
    ) const;

public:

    TypeName("GriffinBoundaryLayerEdge");

    explicit GriffinBoundaryLayerEdge(const fvMesh& mesh);

    //- Demand-driven singleton per mesh (objectRegistry-stored)
    static GriffinBoundaryLayerEdge& New(const fvMesh& mesh);

    //- Edge quantities for one wall patch, cached per time index and
    //  recomputed every edgeUpdateInterval time indices
    const PatchEdgeData& edgeData
    (
        const compressible::turbulenceModel& turbModel,
        const fvPatch& patch,
        const GriffinWallFunction::Coeffs& coeffs
    );

    virtual bool writeData(Ostream&) const
    {
        return true;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
