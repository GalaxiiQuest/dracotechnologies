#!/usr/bin/env python3
"""
Dragon Engineering: flight-scaling model
=========================================

Companion code for:
    Choueiri, A. G. (2026). Dragon Engineering: A Speculative Biology
    Thought Experiment. Preprint.

What it does
------------
Scales a flying animal isometrically (same shape, different size) from the
giant pterosaur Quetzalcoatlus and compares two quantities at each body mass:

  * P_req : mechanical power needed for level flapping flight at the
            minimum-power speed (Pennycuick-style model).
  * P_av  : power the flight muscles can deliver
            (muscle mass x work per beat x wingbeat frequency).

The ratio P_av / P_req tells whether flapping flight is possible (>= 1).
The script reproduces Tables 2, 3 and 6 and Figure 3 of the paper.

This is an order-of-magnitude model. Every assumption is listed in
PARAMETERS below and in the paper's Methods section.

Usage
-----
    python3 dragon_flight_model.py            # print all tables
    python3 dragon_flight_model.py --plot     # also save fig_power.pdf

Requires Python 3.8+. Plotting needs matplotlib (optional).
"""

import math
import sys

# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------
G = 9.81            # gravitational acceleration, m/s^2
RHO = 1.225         # air density at sea level, kg/m^3

# Reference animal: Quetzalcoatlus planform (Witton 2008, as tabulated in
# Witton & Habib 2010, Table 1).
M_REF = 259.06      # mass, kg
B_REF = 9.64        # wingspan, m
S_REF = 11.36       # wing area, m^2

K_INDUCED = 1.2     # induced-drag factor
CD_BODY = 0.1       # body drag coefficient
PROFILE_RATIO = 1.2 # profile power as a multiple of (induced + parasite)
CL_MAX = 2.2        # maximum lift coefficient for stall speed (membrane wing)

SUSTAINED_ANCHOR_KG = 16.0   # 17% muscle animal's sustained limit (calibration)
ANCHOR_MUSCLE_FRACTION = 0.17
BURST_FACTOR = 2.5           # burst (anaerobic) / sustained output


# ---------------------------------------------------------------------------
# MODEL
# ---------------------------------------------------------------------------
def geometry(m):
    """Isometric scaling from the reference animal: span ~ m^(1/3), area ~ m^(2/3)."""
    s = (m / M_REF) ** (1.0 / 3.0)
    return B_REF * s, S_REF * s * s


def body_frontal_area(m):
    """Body frontal area (m^2), Pennycuick's empirical bird relation."""
    return 0.00813 * m ** 0.666


def power_required(m):
    """Minimum mechanical power for level flight (W) and the speed (m/s) at which it occurs."""
    b, _ = geometry(m)
    disc = math.pi * b * b / 4.0           # actuator-disc area
    weight = m * G
    s_body = body_frontal_area(m)
    best_p, best_v = float("inf"), 0.0
    for i in range(10, 1200):              # speeds 1.0 to 119.9 m/s
        v = i / 10.0
        p_induced = K_INDUCED * weight ** 2 / (2.0 * RHO * v * disc)
        p_parasite = 0.5 * RHO * v ** 3 * s_body * CD_BODY
        p = p_induced + p_parasite
        if p < best_p:
            best_p, best_v = p, v
    return best_p * (1.0 + PROFILE_RATIO), best_v


def wingbeat_frequency(m):
    """Pennycuick's empirical wingbeat-frequency formula (Hz)."""
    b, s = geometry(m)
    return (m ** (3 / 8) * G ** 0.5 * b ** (-23 / 24)
            * s ** (-1 / 3) * RHO ** (-3 / 8))


def power_available(m, muscle_fraction, q):
    """Muscle power (W) = muscle mass x work per kg per beat (J/kg) x frequency (Hz)."""
    return muscle_fraction * m * q * wingbeat_frequency(m)


def mass_limit(muscle_fraction, q):
    """Largest mass (kg) at which available power still meets required power (bisection on log scale)."""
    lo, hi = 0.1, 1e5
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        p_req, _ = power_required(mid)
        if power_available(mid, muscle_fraction, q) > p_req:
            lo = mid
        else:
            hi = mid
    return lo


def calibrate_q(anchor_kg=SUSTAINED_ANCHOR_KG, fraction=ANCHOR_MUSCLE_FRACTION):
    """Find work per kg per beat (J/kg) that puts the sustained limit at anchor_kg."""
    lo, hi = 1.0, 500.0
    for _ in range(100):
        q = (lo + hi) / 2.0
        if mass_limit(fraction, q) < anchor_kg:
            lo = q
        else:
            hi = q
    return q


def stall_speed(m):
    _, s = geometry(m)
    return math.sqrt(2.0 * m * G / (RHO * s * CL_MAX))


# ---------------------------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------------------------
def main():
    q_sus = calibrate_q()
    q_burst = BURST_FACTOR * q_sus
    print(f"Calibrated sustained work per beat Q = {q_sus:.1f} J/kg muscle\n")

    print("Table 2. Mass limits for flapping flight")
    print(f"{'Muscle':>8} {'Sustained (kg)':>16} {'Burst (kg)':>12}")
    for f in (0.17, 0.25, 0.30):
        print(f"{f:>8.0%} {mass_limit(f, q_sus):>16.0f} {mass_limit(f, q_burst):>12.0f}")

    print("\nTable 3. Geometry at three sizes")
    print(f"{'Mass':>7} {'Span m':>7} {'Area m2':>8} {'Load N/m2':>10} {'Stall m/s':>10} {'Beats/s':>8}")
    for m in (500, 1000, 5000):
        b, s = geometry(m)
        print(f"{m:>7} {b:>7.1f} {s:>8.1f} {m * G / s:>10.0f} {stall_speed(m):>10.1f} {wingbeat_frequency(m):>8.2f}")

    print("\nFigure 3 data (25% flight muscle)")
    print(f"{'Mass kg':>8} {'P_req W':>10} {'Sus ratio':>10} {'Burst ratio':>12}")
    masses = (10, 70, 259, 500, 1000, 5000)
    fig_rows = []
    for m in masses:
        p_req, _ = power_required(m)
        rs = power_available(m, 0.25, q_sus) / p_req
        rb = power_available(m, 0.25, q_burst) / p_req
        fig_rows.append((m, rs, rb))
        print(f"{m:>8} {p_req:>10.0f} {rs:>10.2f} {rb:>12.2f}")

    print("\nTable 6. Sensitivity of the burst mass limit (kg) at 25% flight muscle")
    anchors = (12.0, 16.0, 22.0)
    factors = (2.0, 2.5, 3.0)
    print(f"{'Anchor kg':>10} " + " ".join(f"{'burst x' + str(x):>11}" for x in factors))
    for a in anchors:
        q = calibrate_q(anchor_kg=a)
        cells = " ".join(f"{mass_limit(0.25, x * q):>11.0f}" for x in factors)
        print(f"{a:>10.0f} {cells}")

    print("\nTable 6b. Burst mass limit (kg) by muscle fraction and burst factor (anchor 16 kg)")
    print(f"{'Muscle':>8} " + " ".join(f"{'burst x' + str(x):>11}" for x in factors))
    for f in (0.17, 0.25, 0.30, 0.35):
        cells = " ".join(f"{mass_limit(f, x * q_sus):>11.0f}" for x in factors)
        print(f"{f:>8.0%} {cells}")

    if "--plot" in sys.argv:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.0, 3.0))
        ax.axhline(1, color="0.35", ls="--", lw=0.8)
        ax.plot([r[0] for r in fig_rows], [r[1] for r in fig_rows], "o-", color="0.6", label="Sustained (aerobic)")
        ax.plot([r[0] for r in fig_rows], [r[2] for r in fig_rows], "o-", color="#1f5fae", label="Burst (~1 min, anaerobic)")
        ax.set_xscale("log")
        ax.set_xlabel("Body mass (kg, log scale)")
        ax.set_ylabel("Power available / needed")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig("fig_power.pdf")
        print("\nSaved fig_power.pdf")


if __name__ == "__main__":
    main()
