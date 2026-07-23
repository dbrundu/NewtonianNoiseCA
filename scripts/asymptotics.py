"""
asymptotics.py -- analytic asymptotic behaviours from Brundu et al. (2022), used
as cross-checks overlaid on the numerical curves (as the 2022 paper did).

All functions return sqrt(S_h)-shaped curves anchored to one numerical point
(anchor = (x_a, y_a)); only the *shape* is physical here, since the absolute
prefactors involve E, E_T, alpha, ... The anchoring makes the slope / envelope
directly comparable to the data.

Regimes (for sqrt(S_h), with omega = 2 pi f):

* Wind-advection dominated (HI, frozen limit, omega r0 / U >> 1), Eq. (frozen ld):
      S_g   ~ omega^{-25/6} exp(-2 omega r0 / U)
      => sqrt(S_h) ~ omega^{-49/12} exp(-omega r0 / U)          (slope ~ -4.08)

* Eddy-decay dominated power-law tail (HI weak wind large r0, and PBL large
  omega), Eqs. (weak wind), (large omega):  S_g ~ omega^{-8}
      => sqrt(S_h) ~ f^{-6}

* Depth scaling at large r0 (model independent), Eq. general scaling:
      S_g ~ 1/r0^2   =>   sqrt(S_h) ~ 1/r0
"""

import math

TAIL_SLOPE = -6.0          # sqrt(S_h) power-law tail in frequency
FROZEN_SLOPE = -49.0 / 12  # sqrt(S_h) power in the frozen wind-dominated envelope


def powerlaw(f, slope, anchor):
    fa, ya = anchor
    return [ya * (fi / fa) ** slope for fi in f]


def frozen_envelope(f, U, r0, anchor):
    """HI frozen wind-dominated shape: f^{-49/12} exp(-omega r0 / U)."""
    fa, ya = anchor
    out = []
    for fi in f:
        shape = (fi / fa) ** FROZEN_SLOPE * math.exp(-2 * math.pi * (fi - fa) * r0 / U)
        out.append(ya * shape)
    return out


def inv_r0(r0, anchor):
    """Model-independent large-depth scaling sqrt(S_h) ~ 1/r0."""
    ra, ya = anchor
    return [ya * (ra / ri) for ri in r0]
