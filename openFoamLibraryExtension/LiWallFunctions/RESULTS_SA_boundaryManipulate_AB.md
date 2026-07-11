# SA `nuTilda1` matrix-constraint A/B test

Date: 2026-07-11

## Change

The local OpenFOAM v2412
`SpalartAllmarasBase::correct()` was changed from

```cpp
nuTildaEqn.ref().relax();
fvOptions.constrain(nuTildaEqn.ref());
solve(nuTildaEqn);
```

to

```cpp
nuTildaEqn.ref().relax();
fvOptions.constrain(nuTildaEqn.ref());
nuTildaEqn.ref().boundaryManipulate(nuTilda_.boundaryFieldRef());
solve(nuTildaEqn);
```

The owning `libcompressibleTurbulenceModels.so` was rebuilt.  Its patched
SHA-256 is
`0803584d0b89f814ad14245641486a4bf4536dab3f2a193e510d55209da54854`.
The original local library was also preserved for a controlled A/B; its
SHA-256 is
`88f20d810205997bf1535e45c71872b475ab1a7f250ba9b762bbaeefb11bf93a`.

## Test design

Two serial 5000-step calculations were started from identical clean fields and
used identical case dictionaries.  Serial execution was selected because this
WSL installation was much faster for the 24,000-cell case than MPI.  The only
runtime difference was the directory prepended to `LD_LIBRARY_PATH`:

- unpatched: `/home/mgnadlinger/OpenFOAM/mgnadlinger-v2412/backups/SA-boundaryManipulate/unpatched`;
- patched: `/home/mgnadlinger/OpenFOAM/mgnadlinger-v2412/backups/SA-boundaryManipulate/patched`.

The cases are:

- `/home/mgnadlinger/WF_validation/Li_et_al_SA_fix/alphaT_0p05_unpatchedSA`;
- `/home/mgnadlinger/WF_validation/Li_et_al_SA_fix/alphaT_0p05_patchedSA`.

The local unpatched control reproduces the earlier cluster result: its final
area averages differ by less than `6e-5%` in heat flux and `3e-5%` in shear.
This makes the local pair a clean test of the SA change despite the two
environments reporting different OpenFOAM build hashes.

## Was Eq. (39) actually imposed?

For every wall face, the written owner-cell `nuTilda` was compared with the
diagnostic `nuTilda1` used by `matrix.setValues()`.

| metric at time 5000 | unpatched | patched |
|---|---:|---:|
| faces | 150 | 150 |
| mean relative mismatch | `5.10975e-1` | `1.54806e-7` |
| maximum relative mismatch | `1.12422` | `5.21162e-6` |
| mean actual/target ratio | `0.840565` | `1.00000003` |

At `x=0.801195 m` specifically:

| run | actual first-cell `nuTilda` | Li `nuTilda1` | relative difference |
|---|---:|---:|---:|
| unpatched | `0.022551148432` | `0.0187034294941` | `+20.5723%` |
| patched | `0.0186438499328` | `0.0186438507568` | `-4.42e-6%` |

The added callback therefore activates the intended constraint.

## Converged wall outputs

Both wall outputs are stationary over the last 100 steps.  The last-1000 mean
flow residual is `2.866e-5` unpatched and `3.458e-5` patched.

| area average | resolved SA at 50000 | unpatched | patched | patch effect |
|---|---:|---:|---:|---:|
| `abs(qw)`, MW/m2 | 29.82338 | 33.914461 | 34.325711 | `+1.2126%` |
| `tauw`, Pa | 7520.763 | 8171.348 | 8257.586 | `+1.0554%` |
| Li `yPlus` | n/a | 125.8402 | 126.4105 | `+0.4532%` |

Consequently, the full-wall discrepancies versus the resolved SA solution
change as follows:

| error versus resolved SA | unpatched | patched |
|---|---:|---:|
| area-average heat flux | `+13.7177%` | `+15.0966%` |
| area-average shear | `+8.6505%` | `+9.7972%` |

The fix is required for implementation fidelity, but it does not explain or
remove the wall-output discrepancy.

## Spatial split

The area averages hide two opposite responses:

- Over `x >= 0.1 m`, the integrated heat-flux error improves from `+13.7487%`
  to `+12.8252%`, and shear improves from `+9.5127%` to `+8.4885%`.
- At `x=0.801195 m`, heat-flux error improves from `+13.2244%` to
  `+12.3667%`, and shear from `+9.7057%` to `+8.7831%`.
- Near `x=0.049196 m`, the response reverses sharply.  The heat-flux ratio to
  resolved SA moves from `0.8242` to `1.2860`; the shear ratio moves from
  `0.7929` to `1.1842`.  This narrow leading/transitional region raises the
  full-wall averages even though the developed region improves slightly.

All 150 Li faces retain status 0, no temperature limit is reached, and heat
flux and shear are delivered accurately.  The localized change is therefore a
closure/state sensitivity exposed by imposing Eq. (39), not a boundary-flux
delivery failure.

At `x=0.8 m`, the resolved SA `nuTilda` interpolated to the coarse first-cell
height is about `0.022155`.  The unpatched value is only `+1.79%` from it,
whereas the enforced Li target is about `-15.85%` low.  Thus the formerly
unconstrained SA solution happened to match the resolved turbulence variable
better there than Eq. (39) does.

## Remaining closure-error decomposition at `x=0.8 m`

Evaluating the implementation a-priori with resolved SA states at the coarse
first- and second-cell locations gives heat flux `+18.11%` and shear `+7.79%`
relative to the resolved wall values.  This separates into two effects:

1. Eq. (42), with `omegaH=0.72`, supplies the `+7.79%` shear error.
2. The Eq. (40)/(41) thermal factor `qw/tauw` is another `+9.57%` high.
   Their product gives the `+18.11%` heat-flux error.

For the thermal part, Eq. (40) infers `rg=0.7912`; resolved `qw`, `tauw`, and
the same reference state imply `rg=0.6880`, a `+15.0%` error.  The underlying
enthalpy at the first point differs from the relation by only about `2.16%`,
but the two-point inversion amplifies it because `U1/Uref=0.846` and the
denominator is `x(1-x)`.

Changing the reference from the second coarse point to progressively deeper
points does not cure this case: a-priori heat-flux error rises from `18.1%` to
about `22.4%`, while shear remains `6.6--7.8%` high.  The second cell is the
least-bad tested reference.

Replacing `omegaH=0.72` by the effective Sutherland exponent near `0.525`
reduces the a-priori shear error to about `5.47%` and heat-flux error to
`15.57%`.  Thus the reacting-air value is a secondary contributor, not the
dominant thermal error.  Setting `omegaH=0` would nearly remove the shear error
but still leaves heat flux `9.55%` high; it is a diagnostic decomposition, not
a justified model change.

## Conclusion

The OpenFOAM change works and should remain if the intended implementation is
to prescribe Li Eq. (39) in the first off-wall cells.  Its effect is not the
desired correction of the validation discrepancy:

1. downstream `qw` and `tauw` improve by roughly one percentage point;
2. the leading/transitional region develops a much larger overprediction;
3. full-wall heat flux and shear errors increase by about 1.38 and 1.15
   percentage points, respectively;
4. the dominant remaining issue is upstream target/closure behavior, not the
   Neumann or shear-flux delivery.

The comparison utility is
`/home/mgnadlinger/WF_validation/Li_et_al_SA_fix/check_nuTilda1_constraint.py`.
