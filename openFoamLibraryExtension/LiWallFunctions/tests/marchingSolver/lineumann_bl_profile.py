#!/usr/bin/env python3
"""Generate a synthetic LiNeumann-style flat-plate boundary-layer profile.

This is a reduced-order profile generator, not a CFD replacement.  It marches
the compressible zero-pressure-gradient momentum-integral equation and closes
skin friction with the Li/White wall-law used by the local LiNeumann wall
function.  The thermal profile is a Crocco/Walz relation written in the same
quadratic enthalpy form used by the Li wall-function implementation.

The output CSV intentionally includes common ParaView/OpenFOAM-like columns:
Points_0, Points_1, Points_2, T, U_0, U_1, U_2, U_Magnitude, nut, nuTilda,
the LiNeumann-implied wall heat flux qw, heat-transfer coefficient Ch, and
friction Reynolds number Re_tau.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from typing import Iterable


ROOTVSMALL = 1.0e-12


@dataclass(frozen=True)
class Gas:
    gamma: float = 1.4
    r: float = 287.052874
    pr: float = 0.72
    sutherland_as: float = 1.458e-6
    sutherland_ts: float = 110.4

    @property
    def cp(self) -> float:
        return self.gamma * self.r / (self.gamma - 1.0)

    def mu(self, temperature: float) -> float:
        t = max(temperature, ROOTVSMALL)
        return self.sutherland_as * t * math.sqrt(t) / (t + self.sutherland_ts)

    def speed_of_sound(self, temperature: float) -> float:
        return math.sqrt(self.gamma * self.r * temperature)


@dataclass(frozen=True)
class LiCoeffs:
    kappa: float = 0.41
    b: float = 6.0
    omega_h: float = 0.72
    cv1: float = 7.1


@dataclass(frozen=True)
class EdgeState:
    mach: float
    pressure: float
    temperature: float
    wall_temperature: float
    gas: Gas

    @property
    def velocity(self) -> float:
        return self.mach * self.gas.speed_of_sound(self.temperature)

    @property
    def density(self) -> float:
        return self.pressure / (self.gas.r * self.temperature)

    @property
    def wall_density(self) -> float:
        return self.pressure / (self.gas.r * self.wall_temperature)

    @property
    def wall_viscosity(self) -> float:
        return self.gas.mu(self.wall_temperature)

    @property
    def recovery_factor(self) -> float:
        return self.gas.pr ** (1.0 / 3.0)

    @property
    def adiabatic_wall_temperature(self) -> float:
        return self.temperature * (
            1.0
            + self.recovery_factor
            * (self.gas.gamma - 1.0)
            * 0.5
            * self.mach
            * self.mach
        )


@dataclass(frozen=True)
class BoundaryLayerState:
    x: float
    theta: float
    delta: float
    displacement_thickness: float
    u_tau: float
    cf: float
    re_x: float
    re_theta: float
    re_tau: float
    beta: float
    gamma_law: float
    t_aw: float
    qw: float
    ch: float
    delta_99: float


def clamp_unit(value: float) -> float:
    return max(-1.0, min(1.0, value))


def y_white_plus(
    u_plus: float,
    beta: float,
    gamma_law: float,
    coeffs: LiCoeffs,
) -> float:
    discriminant = beta * beta + 4.0 * gamma_law
    if gamma_law <= ROOTVSMALL or discriminant <= ROOTVSMALL:
        return math.exp(-coeffs.kappa * coeffs.b + min(50.0, coeffs.kappa * u_plus))

    q = math.sqrt(discriminant)
    phi = math.asin(clamp_unit(-beta / q))
    arg = clamp_unit((2.0 * gamma_law * u_plus - beta) / q)
    exponent = coeffs.kappa / math.sqrt(gamma_law) * (math.asin(arg) - phi)
    exponent = min(50.0, max(-50.0, exponent))
    return math.exp(exponent - coeffs.kappa * coeffs.b)


def d_y_white_du_plus(
    u_plus: float,
    beta: float,
    gamma_law: float,
    coeffs: LiCoeffs,
) -> float:
    discriminant = beta * beta + 4.0 * gamma_law
    if gamma_law <= ROOTVSMALL or discriminant <= ROOTVSMALL:
        return coeffs.kappa * math.exp(
            -coeffs.kappa * coeffs.b + min(50.0, coeffs.kappa * u_plus)
        )

    q = math.sqrt(discriminant)
    arg = clamp_unit((2.0 * gamma_law * u_plus - beta) / q)
    bracket = max(1.0e-14, 1.0 - arg * arg)
    return (
        2.0
        * y_white_plus(u_plus, beta, gamma_law, coeffs)
        * coeffs.kappa
        * math.sqrt(gamma_law)
        / q
        / math.sqrt(bracket)
    )


def li_law_y_plus(
    u_plus: float,
    beta: float,
    gamma_law: float,
    coeffs: LiCoeffs,
) -> float:
    y0 = math.exp(-coeffs.kappa * coeffs.b)
    ku = coeffs.kappa * u_plus

    if gamma_law <= ROOTVSMALL and abs(beta) <= ROOTVSMALL:
        return u_plus + y0 * (
            math.exp(min(50.0, ku))
            - (1.0 + ku + 0.5 * ku * ku + (ku**3) / 6.0)
        )

    u_gm = u_plus * (
        1.0
        + 0.5 * coeffs.omega_h * beta * u_plus
        - (coeffs.omega_h * gamma_law * u_plus * u_plus) / 3.0
    )
    return u_gm + y_white_plus(u_plus, beta, gamma_law, coeffs) - y0 * (
        1.0 + ku + 0.5 * ku * ku + (ku**3) / 6.0
    )


def d_li_law_y_plus_du_plus(
    u_plus: float,
    beta: float,
    gamma_law: float,
    coeffs: LiCoeffs,
) -> float:
    y0 = math.exp(-coeffs.kappa * coeffs.b)
    ku = coeffs.kappa * u_plus
    derivative = (
        1.0
        + coeffs.omega_h * beta * u_plus
        - coeffs.omega_h * gamma_law * u_plus * u_plus
        + d_y_white_du_plus(u_plus, beta, gamma_law, coeffs)
        - coeffs.kappa * y0 * (1.0 + ku + 0.5 * ku * ku)
    )
    return max(derivative, ROOTVSMALL)


def solve_nu_tilda(
    nut_target: float,
    nu: float,
    cv1: float,
    max_iter: int = 80,
    tolerance: float = 1.0e-8,
) -> float:
    if nut_target <= 0.0:
        return 0.0

    nu_eff = max(nu, ROOTVSMALL)

    def residual(nu_tilda: float) -> float:
        chi = nu_tilda / nu_eff
        chi3 = chi**3
        fv1 = chi3 / (chi3 + cv1**3)
        return nu_tilda * fv1 - nut_target

    lo = 0.0
    hi = max(nut_target, nu_eff)
    for _ in range(80):
        if residual(hi) >= 0.0:
            break
        hi *= 2.0

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if residual(mid) < 0.0:
            lo = mid
        else:
            hi = mid
        if abs(hi - lo) / max(mid, ROOTVSMALL) < tolerance:
            break

    return 0.5 * (lo + hi)


def wall_law_parameters(edge: EdgeState, u_tau: float) -> tuple[float, float, float]:
    gas = edge.gas
    h_w = gas.cp * edge.wall_temperature
    h_e = gas.cp * edge.temperature
    t_aw = edge.adiabatic_wall_temperature
    h_rg = gas.cp * t_aw
    u_e = edge.velocity
    rg = 2.0 * (h_rg - h_e) / max(u_e * u_e, ROOTVSMALL)
    beta = u_tau * (h_rg - h_w) / max(u_e * h_w, ROOTVSMALL)
    gamma_law = rg * u_tau * u_tau / max(2.0 * h_w, ROOTVSMALL)
    return beta, gamma_law, t_aw


def crocco_walz_temperature(edge: EdgeState, speed_ratio: float) -> float:
    s = max(0.0, min(1.0, speed_ratio))
    return (
        edge.wall_temperature
        + (edge.adiabatic_wall_temperature - edge.wall_temperature) * s
        + (edge.temperature - edge.adiabatic_wall_temperature) * s * s
    )


def wall_heat_flux(edge: EdgeState, u_tau: float) -> float:
    """Return paper-sign q_w [W/m^2], positive from the fluid into the wall."""
    tau_w = edge.wall_density * u_tau * u_tau
    h_rg = edge.gas.cp * edge.adiabatic_wall_temperature
    h_w = edge.gas.cp * edge.wall_temperature
    return (h_rg - h_w) * tau_w / max(
        edge.gas.pr * edge.velocity,
        ROOTVSMALL,
    )


def heat_transfer_coefficient(edge: EdgeState, u_tau: float) -> float:
    """Return Ch = qw/[rho_inf*U_inf*cp*(Tr - Tw)]."""
    tau_w = edge.wall_density * u_tau * u_tau
    denominator = edge.gas.pr * edge.density * edge.velocity * edge.velocity
    return tau_w / max(denominator, ROOTVSMALL)


def delta_at_speed_ratio(
    edge: EdgeState,
    coeffs: LiCoeffs,
    u_tau: float,
    beta: float,
    gamma_law: float,
    speed_ratio: float,
) -> float:
    """Return y where the synthetic profile reaches U/U_inf = speed_ratio."""
    ratio = max(0.0, min(1.0, speed_ratio))
    u_plus = ratio * edge.velocity / max(u_tau, ROOTVSMALL)
    wall_length = edge.wall_viscosity / max(
        edge.wall_density * u_tau,
        ROOTVSMALL,
    )
    return wall_length * li_law_y_plus(u_plus, beta, gamma_law, coeffs)


def integrate_profile_moments(
    edge: EdgeState,
    coeffs: LiCoeffs,
    u_tau: float,
    n_quad: int,
) -> tuple[float, float, float]:
    beta, gamma_law, _ = wall_law_parameters(edge, u_tau)
    u_e = edge.velocity
    u_plus_edge = u_e / max(u_tau, ROOTVSMALL)
    scale = edge.wall_viscosity / max(edge.wall_density * u_tau, ROOTVSMALL)
    delta = scale * li_law_y_plus(u_plus_edge, beta, gamma_law, coeffs)

    theta_integral = 0.0
    disp_integral = 0.0
    prev_u = 0.0
    prev_theta = 0.0
    prev_disp = scale * d_li_law_y_plus_du_plus(0.0, beta, gamma_law, coeffs)

    for i in range(1, n_quad + 1):
        u_plus = u_plus_edge * i / n_quad
        speed_ratio = u_plus / max(u_plus_edge, ROOTVSMALL)
        temperature = crocco_walz_temperature(edge, speed_ratio)
        density_ratio = edge.temperature / max(temperature, ROOTVSMALL)
        dy_du = scale * d_li_law_y_plus_du_plus(u_plus, beta, gamma_law, coeffs)
        theta_value = density_ratio * speed_ratio * (1.0 - speed_ratio) * dy_du
        disp_value = (1.0 - density_ratio * speed_ratio) * dy_du
        theta_integral += 0.5 * (theta_value + prev_theta) * (u_plus - prev_u)
        disp_integral += 0.5 * (disp_value + prev_disp) * (u_plus - prev_u)
        prev_u = u_plus
        prev_theta = theta_value
        prev_disp = disp_value

    return theta_integral, disp_integral, delta


def state_from_theta(
    x: float,
    theta_target: float,
    edge: EdgeState,
    coeffs: LiCoeffs,
    n_quad: int,
) -> BoundaryLayerState:
    u_e = edge.velocity

    def residual(u_tau: float) -> tuple[float, float, float, float]:
        theta, disp, delta = integrate_profile_moments(edge, coeffs, u_tau, n_quad)
        return theta - theta_target, theta, disp, delta

    lo = u_e / 5000.0
    hi = u_e / 3.0
    f_lo, _, _, _ = residual(lo)
    f_hi, _, _, _ = residual(hi)

    for _ in range(80):
        if f_lo > 0.0:
            break
        lo *= 0.5
        f_lo, _, _, _ = residual(lo)

    for _ in range(80):
        if f_hi < 0.0:
            break
        hi *= 1.5
        f_hi, _, _, _ = residual(hi)

    if f_lo <= 0.0 or f_hi >= 0.0:
        raise RuntimeError(
            "Could not bracket u_tau for theta="
            f"{theta_target:.6e}; residuals are {f_lo:.6e}, {f_hi:.6e}"
        )

    theta = disp = delta = float("nan")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid, theta, disp, delta = residual(mid)
        if f_mid > 0.0:
            lo = mid
        else:
            hi = mid
        if abs(hi - lo) / max(mid, ROOTVSMALL) < 1.0e-9:
            break

    u_tau = 0.5 * (lo + hi)
    beta, gamma_law, t_aw = wall_law_parameters(edge, u_tau)
    tau_w = edge.wall_density * u_tau * u_tau
    cf = 2.0 * tau_w / max(edge.density * edge.velocity * edge.velocity, ROOTVSMALL)
    mu_e = edge.gas.mu(edge.temperature)
    re_x = edge.density * edge.velocity * x / mu_e
    re_theta = edge.density * edge.velocity * theta_target / mu_e
    delta_99 = delta_at_speed_ratio(
        edge,
        coeffs,
        u_tau,
        beta,
        gamma_law,
        0.99,
    )
    re_tau = (
        edge.wall_density * u_tau * delta_99
        / max(edge.wall_viscosity, ROOTVSMALL)
    )
    return BoundaryLayerState(
        x=x,
        theta=theta_target,
        delta=delta,
        displacement_thickness=disp,
        u_tau=u_tau,
        cf=cf,
        re_x=re_x,
        re_theta=re_theta,
        re_tau=re_tau,
        beta=beta,
        gamma_law=gamma_law,
        t_aw=t_aw,
        qw=wall_heat_flux(edge, u_tau),
        ch=heat_transfer_coefficient(edge, u_tau),
        delta_99=delta_99,
    )


def initial_theta(edge: EdgeState, x_start: float) -> float:
    mu_e = edge.gas.mu(edge.temperature)
    re_x = max(edge.density * edge.velocity * x_start / mu_e, 1.0)
    return 0.036 * x_start * re_x ** (-0.2)


def march_to_x(
    x_target: float,
    edge: EdgeState,
    coeffs: LiCoeffs,
    x_start: float,
    n_steps: int,
    n_quad: int,
) -> BoundaryLayerState:
    if x_target <= x_start:
        raise ValueError("x_target must be larger than x_start")

    theta = initial_theta(edge, x_start)
    x = x_start

    for step in range(n_steps):
        fraction0 = step / n_steps
        fraction1 = (step + 1) / n_steps
        x0 = x_start * (x_target / x_start) ** fraction0
        x1 = x_start * (x_target / x_start) ** fraction1
        if step == 0:
            x = x0

        state0 = state_from_theta(x, theta, edge, coeffs, n_quad)
        dtheta_dx0 = 0.5 * state0.cf
        theta_mid = theta + 0.5 * (x1 - x0) * dtheta_dx0
        x_mid = 0.5 * (x0 + x1)
        state_mid = state_from_theta(x_mid, theta_mid, edge, coeffs, n_quad)
        theta += (x1 - x0) * 0.5 * state_mid.cf
        x = x1

    return state_from_theta(x_target, theta, edge, coeffs, n_quad)


def solve_u_plus_for_y_plus(
    y_plus_target: float,
    u_plus_edge: float,
    beta: float,
    gamma_law: float,
    coeffs: LiCoeffs,
) -> float:
    if y_plus_target <= 0.0:
        return 0.0
    edge_y_plus = li_law_y_plus(u_plus_edge, beta, gamma_law, coeffs)
    if y_plus_target >= edge_y_plus:
        return u_plus_edge

    lo = 0.0
    hi = u_plus_edge
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if li_law_y_plus(mid, beta, gamma_law, coeffs) < y_plus_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def profile_rows(
    state: BoundaryLayerState,
    edge: EdgeState,
    coeffs: LiCoeffs,
    y_values: Iterable[float],
    mu_t_ratio_inf: float,
) -> list[dict[str, float]]:
    beta = state.beta
    gamma_law = state.gamma_law
    u_plus_edge = edge.velocity / max(state.u_tau, ROOTVSMALL)
    y_plus_scale = edge.wall_density * state.u_tau / edge.wall_viscosity
    rows = []
    edge_y_plus = li_law_y_plus(u_plus_edge, beta, gamma_law, coeffs)

    for y in y_values:
        y_plus = y * y_plus_scale
        outside_layer = y_plus >= edge_y_plus
        if outside_layer:
            u_plus = u_plus_edge
            speed = edge.velocity
            speed_ratio = 1.0
            temperature = edge.temperature
        else:
            u_plus = solve_u_plus_for_y_plus(
                y_plus,
                u_plus_edge,
                beta,
                gamma_law,
                coeffs,
            )
            speed = min(edge.velocity, state.u_tau * u_plus)
            speed_ratio = speed / edge.velocity
            temperature = crocco_walz_temperature(edge, speed_ratio)
        density = edge.pressure / (edge.gas.r * temperature)
        mu = edge.gas.mu(temperature)
        nu = mu / max(density, ROOTVSMALL)
        if outside_layer:
            nut = mu_t_ratio_inf * nu
            mut = density * nut
        else:
            law_derivative = d_li_law_y_plus_du_plus(
                u_plus,
                beta,
                gamma_law,
                coeffs,
            )
            mut = max(0.0, edge.wall_viscosity * law_derivative - mu)
            nut = mut / max(density, ROOTVSMALL)
        nu_tilda = solve_nu_tilda(nut, nu, coeffs.cv1)
        mach = speed / edge.gas.speed_of_sound(temperature)
        rows.append(
            {
                "x": state.x,
                "y": y,
                "yPlus": y_plus,
                "U": speed,
                "T": temperature,
                "rho": density,
                "mu": mu,
                "nu": nu,
                "mut": mut,
                "nut": nut,
                "nuTilda": nu_tilda,
                "qw": state.qw,
                "Ch": state.ch,
                "Re_tau": state.re_tau,
                "Ma": mach,
                "Points_0": state.x,
                "Points_1": y,
                "Points_2": 0.0,
                "U_0": speed,
                "U_1": 0.0,
                "U_2": 0.0,
                "U_Magnitude": speed,
            }
        )
    return rows


def stretched_y_values(y_max: float, n: int, stretch: float) -> list[float]:
    if n < 2:
        return [0.0]
    if abs(stretch - 1.0) < 1.0e-12:
        return [y_max * i / (n - 1) for i in range(n)]
    return [y_max * (math.exp(stretch * i / (n - 1)) - 1.0) / (math.exp(stretch) - 1.0) for i in range(n)]


def write_csv(path: str, rows: list[dict[str, float]], metadata: dict[str, float]) -> None:
    fields = [
        "x",
        "y",
        "yPlus",
        "Ma",
        "rho",
        "mu",
        "nu",
        "T",
        "U",
        "qw",
        "Ch",
        "Re_tau",
        "mut",
        "nut",
        "nuTilda",
        "Points_0",
        "Points_1",
        "Points_2",
        "U_0",
        "U_1",
        "U_2",
        "U_Magnitude",
    ]
    with open(path, "w", newline="") as handle:
        handle.write("# LiNeumann synthetic boundary-layer profile\n")
        for key, value in metadata.items():
            handle.write(f"# {key}: {value:.12g}\n")
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a LiNeumann-style synthetic U(y), T(y) profile."
    )
    parser.add_argument("--x", type=float, default=0.045, help="Target x location [m].")
    parser.add_argument("--mach", type=float, default=5.84, help="Freestream Mach number.")
    parser.add_argument("--p-inf", type=float, default=697.0, help="Freestream pressure [Pa].")
    parser.add_argument("--t-inf", type=float, default=55.2, help="Freestream static temperature [K].")
    parser.add_argument("--tw", type=float, default=97.5, help="Wall temperature [K].")
    parser.add_argument("--out", default="lineumann_profile_x0p045.csv", help="Output CSV path.")
    parser.add_argument("--n-y", type=int, default=200, help="Number of profile points.")
    parser.add_argument("--y-max", type=float, default=None, help="Maximum y [m]; default is y-max-factor*delta.")
    parser.add_argument("--y-max-factor", type=float, default=1.25, help="y_max/delta if --y-max is omitted.")
    parser.add_argument("--stretch", type=float, default=4.0, help="Exponential y-grid stretch; 1 gives nearly uniform.")
    parser.add_argument("--x-start", type=float, default=1.0e-5, help="Virtual turbulent start location [m].")
    parser.add_argument("--steps", type=int, default=180, help="Momentum-integral marching steps.")
    parser.add_argument("--quad", type=int, default=320, help="Quadrature points for profile moments.")
    parser.add_argument("--gamma", type=float, default=1.4)
    parser.add_argument("--r", type=float, default=287.052874, help="Specific gas constant [J/(kg K)].")
    parser.add_argument("--pr", type=float, default=0.72, help="Molecular Prandtl number for recovery factor.")
    parser.add_argument("--kappa", type=float, default=0.41)
    parser.add_argument("--B", type=float, default=6.0)
    parser.add_argument("--omega-h", type=float, default=0.72)
    parser.add_argument("--cv1", type=float, default=7.1)
    parser.add_argument("--mu-t-ratio", type=float, default=0.001, help="Freestream nut/nu ratio outside the synthetic BL.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gas = Gas(gamma=args.gamma, r=args.r, pr=args.pr)
    coeffs = LiCoeffs(kappa=args.kappa, b=args.B, omega_h=args.omega_h, cv1=args.cv1)
    edge = EdgeState(
        mach=args.mach,
        pressure=args.p_inf,
        temperature=args.t_inf,
        wall_temperature=args.tw,
        gas=gas,
    )
    state = march_to_x(
        x_target=args.x,
        edge=edge,
        coeffs=coeffs,
        x_start=args.x_start,
        n_steps=args.steps,
        n_quad=args.quad,
    )
    y_max = args.y_max if args.y_max is not None else args.y_max_factor * state.delta
    y_values = stretched_y_values(y_max, args.n_y, args.stretch)
    rows = profile_rows(state, edge, coeffs, y_values, args.mu_t_ratio)
    nu_inf = gas.mu(args.t_inf) / edge.density
    nut_inf = args.mu_t_ratio * nu_inf
    metadata = {
        "x": state.x,
        "M_inf": args.mach,
        "p_inf": args.p_inf,
        "T_inf": args.t_inf,
        "Tw": args.tw,
        "U_inf": edge.velocity,
        "rho_inf": edge.density,
        "mu_inf": gas.mu(args.t_inf),
        "nut_inf": nut_inf,
        "nuTilda_inf": solve_nu_tilda(nut_inf, nu_inf, coeffs.cv1),
        "T_aw": state.t_aw,
        "qw": state.qw,
        "Ch": state.ch,
        "theta": state.theta,
        "delta_star": state.displacement_thickness,
        "delta": state.delta,
        "delta_99": state.delta_99,
        "u_tau": state.u_tau,
        "Cf": state.cf,
        "Re_x": state.re_x,
        "Re_theta": state.re_theta,
        "Re_tau": state.re_tau,
        "beta": state.beta,
        "Gamma": state.gamma_law,
        "cv1": coeffs.cv1,
        "mu_t_ratio": args.mu_t_ratio,
    }
    write_csv(args.out, rows, metadata)
    print(f"Wrote {len(rows)} points to {args.out}")
    print(
        "x={x:.6g} m, delta={delta:.6g} m, theta={theta:.6g} m, "
        "Cf={cf:.6g}, u_tau={u_tau:.6g} m/s, qw={qw:.6g} W/m^2, "
        "Ch={ch:.6g}, delta_99={delta_99:.6g} m, Re_tau={re_tau:.6g}, "
        "Re_x={re_x:.6g}".format(
            x=state.x,
            delta=state.delta,
            theta=state.theta,
            cf=state.cf,
            u_tau=state.u_tau,
            qw=state.qw,
            ch=state.ch,
            delta_99=state.delta_99,
            re_tau=state.re_tau,
            re_x=state.re_x,
        )
    )


if __name__ == "__main__":
    main()
