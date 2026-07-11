#include "GriffinBoundaryLayerEdge.H"

#include "DynamicList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"

namespace Foam
{
namespace
{

inline vector tangentialPart(const vector& value, const vector& normal)
{
    return value - normal*(value & normal);
}

} // End anonymous namespace


namespace GriffinWallFunction
{

label findEdge
(
    const UList<scalar>& dist,
    const UList<scalar>& Utan,
    const UList<scalar>& Vn,
    const UList<scalar>& h,
    const UList<scalar>& T,
    const scalar threshold,
    scalar& delta,
    scalar& Ue,
    scalar& Te
)
{
    const label n = dist.size();

    delta = 0;
    Ue = 0;
    Te = 0;

    if (n == 0)
    {
        return 1;
    }

    // Reference stagnation enthalpy: the profile maximum (PRF 2021,
    // P_o,ref analog of Appendix D).
    scalar h0Ref = -GREAT;
    forAll(dist, i)
    {
        const scalar h0 = h[i] + 0.5*(sqr(Utan[i]) + sqr(Vn[i]));
        h0Ref = max(h0Ref, h0);
    }

    const scalar thresholdSqr = sqr(threshold);

    const auto ratioSqr = [&](const label i) -> scalar
    {
        const scalar UISqr =
            max(2*(h0Ref - h[i]) - sqr(Vn[i]), VSMALL);
        return sqr(Utan[i])/UISqr;
    };

    bool armed = false;
    label iPrev = -1;
    scalar rPrev = 0;

    forAll(dist, i)
    {
        const scalar r = ratioSqr(i);

        if (!armed)
        {
            if (r < thresholdSqr)
            {
                armed = true;
                iPrev = i;
                rPrev = r;
            }
            continue;
        }

        if (r >= thresholdSqr)
        {
            const scalar denom = max(r - rPrev, VSMALL);
            const scalar w = min(max((thresholdSqr - rPrev)/denom, scalar(0)), scalar(1));

            delta = dist[iPrev] + w*(dist[i] - dist[iPrev]);
            Ue = Utan[iPrev] + w*(Utan[i] - Utan[iPrev]);
            Te = T[iPrev] + w*(T[i] - T[iPrev]);
            return 0;
        }

        iPrev = i;
        rPrev = r;
    }

    // No crossing: use the outermost sample as the edge estimate.
    delta = dist[n - 1];
    Ue = Utan[n - 1];
    Te = T[n - 1];
    return 1;
}

} // End namespace GriffinWallFunction


defineTypeNameAndDebug(GriffinBoundaryLayerEdge, 0);


GriffinBoundaryLayerEdge::GriffinBoundaryLayerEdge(const fvMesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            GriffinBoundaryLayerEdge::typeName,
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    patchData_()
{}


GriffinBoundaryLayerEdge& GriffinBoundaryLayerEdge::New(const fvMesh& mesh)
{
    auto* ptr =
        mesh.thisDb().getObjectPtr<GriffinBoundaryLayerEdge>(typeName);

    if (!ptr)
    {
        ptr = new GriffinBoundaryLayerEdge(mesh);
        ptr->store();
    }

    return *ptr;
}


void GriffinBoundaryLayerEdge::reconstruct
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const GriffinWallFunction::Coeffs& coeffs,
    PatchEdgeData& pd
) const
{
    const labelUList& faceCells = patch.faceCells();

    const fluidThermo& thermo = turbModel.transport();
    const volVectorField& U = turbModel.U();
    const fvPatchVectorField& Uw = U.boundaryField()[patch.index()];
    const volScalarField& he = thermo.he();
    const volScalarField& p = thermo.p();
    const volScalarField& rho = turbModel.rho();
    const volScalarField& T = thermo.T();
    const bool internalEnergy = he.name() == "e";

    const vectorField nOut(patch.nf());
    const vectorField& Cf = patch.Cf();
    const vectorField& C = mesh_.C();
    const labelListList& cellCells = mesh_.cellCells();

    const label maxSamples = max(coeffs.maxEdgeCells, label(2));

    DynamicList<scalar> dist(maxSamples);
    DynamicList<scalar> Utan(maxSamples);
    DynamicList<scalar> Vn(maxSamples);
    DynamicList<scalar> h(maxSamples);
    DynamicList<scalar> Ts(maxSamples);

    label noCrossingCount = 0;

    forAll(faceCells, facei)
    {
        const vector normal = nOut[facei];
        const vector inward = -normal;
        const vector origin = Cf[facei];
        const vector wallVelocity = Uw[facei];

        dist.clear();
        Utan.clear();
        Vn.clear();
        h.clear();
        Ts.clear();

        label celli = faceCells[facei];

        for (label step = 0; step < maxSamples; ++step)
        {
            const scalar d = (C[celli] - origin) & inward;
            const vector Urel = U[celli] - wallVelocity;

            dist.append(d);
            Utan.append(mag(tangentialPart(Urel, normal)));
            Vn.append(Urel & inward);
            scalar hs = he[celli];
            if (internalEnergy)
            {
                hs += p[celli]/max(rho[celli], VSMALL);
            }
            h.append(hs);
            Ts.append(T[celli]);

            if (d > coeffs.maxEdgeDistance)
            {
                break;
            }

            // Advance to the neighbour with the largest inward-normal
            // progress; positive advancement makes the walk cycle-free.
            // Processor-boundary neighbours are not visible here: the
            // walk simply stops (documented limitation).
            label next = -1;
            scalar best = SMALL;

            const labelList& nbrs = cellCells[celli];
            forAll(nbrs, nbrI)
            {
                const scalar adv = (C[nbrs[nbrI]] - C[celli]) & inward;
                if (adv > best)
                {
                    best = adv;
                    next = nbrs[nbrI];
                }
            }

            if (next == -1)
            {
                break;
            }

            celli = next;
        }

        pd.status[facei] = GriffinWallFunction::findEdge
        (
            dist,
            Utan,
            Vn,
            h,
            Ts,
            coeffs.edgeThreshold,
            pd.delta[facei],
            pd.Ue[facei],
            pd.Te[facei]
        );

        if (pd.status[facei] != 0)
        {
            ++noCrossingCount;
        }
    }

    if (noCrossingCount)
    {
        static bool warned = false;
        if (!warned)
        {
            WarningInFunction
                << "Griffin boundary-layer edge reconstruction found no "
                << "U/U_I threshold crossing on " << noCrossingCount
                << " face(s) of patch " << patch.name()
                << "; the outermost profile sample was used instead. "
                << "Possible causes: maxEdgeCells/maxEdgeDistance too "
                << "small, the wall-normal cell walk hit a processor "
                << "boundary (decompose wall-tangentially), or the "
                << "boundary layer fills the domain. Consider edgeMode "
                << "prescribed." << endl;
            warned = true;
        }
    }
}


const GriffinBoundaryLayerEdge::PatchEdgeData&
GriffinBoundaryLayerEdge::edgeData
(
    const compressible::turbulenceModel& turbModel,
    const fvPatch& patch,
    const GriffinWallFunction::Coeffs& coeffs
)
{
    const label patchi = patch.index();
    const label nFaces = patch.size();

    if (!patchData_.found(patchi))
    {
        patchData_.insert(patchi, PatchEdgeData());
    }

    PatchEdgeData& pd = patchData_[patchi];

    const bool sizeChanged = pd.delta.size() != nFaces;
    if (sizeChanged)
    {
        pd.delta.resize(nFaces, scalar(0));
        pd.Ue.resize(nFaces, scalar(0));
        pd.Te.resize(nFaces, scalar(0));
        pd.status.resize(nFaces, label(0));
        pd.lastTimeIndex = -1;
    }

    if (coeffs.edgeMode == "prescribed")
    {
        if (pd.lastTimeIndex < 0)
        {
            pd.delta = coeffs.delta;
            pd.Ue = coeffs.Ue;
            pd.Te = coeffs.Te;
            pd.status = 2;
            pd.lastTimeIndex = mesh_.time().timeIndex();
        }
        return pd;
    }

    const label timeIndex = mesh_.time().timeIndex();
    const label interval = max(coeffs.edgeUpdateInterval, label(1));

    if
    (
        pd.lastTimeIndex < 0
     || (timeIndex - pd.lastTimeIndex) >= interval
    )
    {
        reconstruct(turbModel, patch, coeffs, pd);
        pd.lastTimeIndex = timeIndex;
    }

    return pd;
}

} // End namespace Foam

// ************************************************************************* //
