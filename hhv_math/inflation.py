"""
inflation.py

Slow-roll inflation helper functions.

All formulas here are standard in the literature for single-field,
canonical slow-roll inflation in Planck units (M_pl = 1 unless stated).

Key relations:

  P_t = 2 H_I^2 / (pi^2 M_pl^2)
  P_s = A_s
  r   = P_t / P_s = 2 H_I^2 / (pi**2 A_s M_pl**2)

So:

  H_I = (pi / sqrt(2)) * sqrt(r * A_s) * M_pl

We also provide the approximate inflationary energy scale:

  V^(1/4) ~ (3^(1/4)) * (H_I)^(1/2) * M_pl^(1/2)

assuming V ~ 3 H_I^2 M_pl^2.
"""

import numpy as np


def H_from_r_As(r, As, Mpl=1.0):
    """
    Compute Hubble scale H_I from tensor-to-scalar ratio r and scalar amplitude A_s.

    Parameters
    ----------
    r : float or array
        Tensor-to-scalar ratio at pivot scale.
    As : float
        Scalar amplitude A_s at pivot (e.g. ~2.1e-9 from Planck).
    Mpl : float
        Reduced Planck mass (default 1.0 for Planck units).

    Returns
    -------
    H_I : float or ndarray
        Inflationary Hubble parameter H_I in the same units as Mpl.

    Formula:
      H_I = (pi / sqrt(2)) * sqrt(r * As) * Mpl
    """
    r = np.asarray(r, dtype=float)
    return (np.pi / np.sqrt(2.0)) * np.sqrt(r * As) * Mpl


def r_from_H_As(HI, As, Mpl=1.0):
    """
    Compute r from H_I and A_s by inverting the same relation.

    Parameters
    ----------
    HI : float or array
        Hubble scale H_I.
    As : float
        Scalar amplitude A_s.
    Mpl : float
        Reduced Planck mass.

    Returns
    -------
    r : float or ndarray
        Tensor-to-scalar ratio.

    Formula:
      r = 2 H_I^2 / (pi^2 A_s M_pl^2)
    """
    HI = np.asarray(HI, dtype=float)
    return 2.0 * (HI**2) / (np.pi**2 * As * Mpl**2)


def V_quarter_from_H(HI, Mpl=1.0):
    """
    Estimate the inflationary energy scale V^(1/4) from H_I.

    Assumes V ~ 3 H_I^2 M_pl^2, so:

      V^(1/4) = (3^(1/4)) * H_I^(1/2) * M_pl^(1/2).

    Parameters
    ----------
    HI : float or array
        Hubble scale H_I.
    Mpl : float
        Reduced Planck mass.

    Returns
    -------
    V_quarter : float or ndarray
        V^(1/4) in the same units as Mpl^(1/2) * HI^(1/2).
    """
    HI = np.asarray(HI, dtype=float)
    return (3.0 ** 0.25) * np.sqrt(HI * Mpl)


def slow_roll_epsilon_from_r(r):
    """
    For canonical single-field inflation, slow-roll parameter epsilon is:

      epsilon = r / 16

    Parameters
    ----------
    r : float or array

    Returns
    -------
    epsilon : float or ndarray
    """
    r = np.asarray(r, dtype=float)
    return r / 16.0
