# FINDINGS: Why nut-based wall shear prescription has no sign-flip problem (unlike alphat)

**Scope:** all four wall-function families — Gao, Li, Wang, Griffin.

**Context:** In `LiWallFunctions/` and `GaoWallFunctions/` we documented that the
alphat-based heat-flux prescription breaks when the first off-wall cell has
`T_1 < T_w`: OpenFOAM computes the wall heat flux as
`q = alphaEff * grad_n(h)|_wall`, and since `alphat >= 0` can only scale the
magnitude, the *sign* of the discrete flux is dictated by `sign(T_1 - T_w)` —
which can disagree with the sign the closure wants. This motivated the switch
to a Neumann temperature BC (`*NeumannWallHeatFluxTemperature`).

The natural follow-up question: the wall functions prescribe the wall *shear
stress* through exactly the same mechanism (a nonnegative `nut` on the wall
patch). Why does the analogous sign problem not occur there — not even under
separation or recirculation?

## 1. How the shear stress is prescribed

The solver's momentum equation computes the wall stress as the diffusive
wall-face flux:

```
tau_solver = (mu_w + rho * nut_w) * grad_n(U)|_wall
```

All four `nut*WallFunction` BCs back-solve this so it reproduces the closure's
`tau_w = rho_w * uTau^2`. The core formula (identical in all four families):

```cpp
nutShear = max(scalar(0), (tauw - muWall*magGradU)/(rhoStress*magGradU));
```

so that `(mu_w + rho*nut) * |grad_n(U_t)| = tau_w` exactly. Locations:

- `LiWallFunctions/LiWallFunction/LiWallFunction.C:794`
- `GaoWallFunctions/GaoWallFunction/GaoWallFunction.C:517`
- `GriffinWallFunctions/GriffinWallFunction/GriffinWallFunction.C:887`
- `WangWallFunction/WangNeumannWallFunction/WangNeumannWallFunction.C:1188`

So: same trick as alphat (nonnegative scalar multiplier on the discrete
gradient), same `max(0, ...)` clip.

## 2. Why the sign problem cannot occur for momentum

The alphat problem exists because **two different physical comparisons** set
the sign:

| | discrete flux sign set by | closure's desired sign set by |
|---|---|---|
| heat flux | `sign(T_1 - T_w)` | `sign(T_aw - T_w)` (wall vs. *recovery* temperature) |
| shear stress | direction of `(U_1 - U_wall)` | direction of `(U_1 - U_wall)` |

In hypersonic flow, `T_1 < T_w` can coexist with the closure wanting heat
*into* the wall (`T_aw > T_w`), so the discrete flux has the wrong sign and no
nonnegative alphat can fix it.

For momentum, both signs come from the **same comparison**:

- The discrete wall gradient is `Uw.snGrad() = (U_1 - U_wall)/delta` (Li takes
  its tangential projection at `LiWallFunction.C:706`; Gao takes
  `mag(Uw.snGrad())` at `GaoWallFunction.C:341`). Its tangential direction is
  *by construction* aligned with the first-cell tangential relative velocity.
- The closure solves `uTau` from `U1 = mag(tangential(U[celli] - wallVelocity))`
  (`LiWallFunction.C:561`), a nonnegative magnitude; the target
  `tau_w = rho_w * uTau^2 >= 0` is implicitly directed opposite that *same*
  relative velocity.

**There is no momentum analogue of a recovery temperature.** The "zero-stress
wall velocity" is exactly the no-slip wall velocity itself. A wall can heat a
fluid that is hotter than the first cell; it can never drag the fluid forward
relative to itself. The closure therefore never requests a stress whose
direction disagrees with the discrete gradient — prescribing shear is a pure
*magnitude* problem, whereas prescribing heat flux was a *signed* problem.

## 3. Separation and recirculation

The right place to look for trouble, but it produces accuracy loss only, never
a sign conflict:

- **Inside a recirculation bubble:** `U_1t` reverses. The discrete gradient
  direction *and* the closure's stress direction reverse together, because
  both derive from the same vector. The wall correctly drags against the
  reversed flow.
- **At separation / reattachment points:** `U1 -> 0`. Handled explicitly: Li
  bails out with a status code when `U1 <= SMALL` (`LiWallFunction.C:583`);
  when `magGradU` vanishes, all families fall back to
  `nutShear = mut1/rhoStress` instead of dividing by zero.
- **Sub-cell backflow** (incipient separation with the reversal below the
  first cell, y+ ~ 30–100): invisible to the discretization —
  `(U_1 - U_w)/delta` is always aligned with `U_1`, so no *inconsistency*
  arises within the scheme. The equilibrium closure is simply physically wrong
  there (`tau_w -> 0`, log-law assumptions invalid) — the well-known accuracy
  limitation of any equilibrium wall function near separation, not a
  BC-mechanism failure.

## 4. The one real (mild) analogue of the clipping

The `max(0, ...)` does bite occasionally: `nut >= 0` means the effective
stress can never fall **below** the laminar stress `mu_w * |grad_n(U)|`. If
the closure wants `tau_w < mu_w * |grad_n(U)|` — strongly accelerated flow,
relaminarization, transients, very low y+ where the fit slightly undershoots —
the prescription clips to the laminar value.

Unlike the alphat case, this error is:

- **bounded** (at most the laminar/target gap),
- **same-signed** (never flips the stress direction),
- **vanishing** as y+ enters the viscous sublayer (there `u+ = y+` makes the
  target and laminar stress coincide).

That asymmetry — momentum can scale *up* but not below laminar, versus heat
flux not being able to flip sign at all — is exactly why the energy side
needed the Neumann rework (`*NeumannWallHeatFluxTemperature`) and the momentum
side did not.

## Conclusion

The nut-based shear prescription is structurally immune to the sign-flip
failure documented for alphat. No Neumann-style rework of the momentum BC is
needed — including for detached/recirculating flow, where the limitation is
the equilibrium closure's physics, not the prescription mechanism.
