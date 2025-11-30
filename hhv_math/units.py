"""
units.py

Basic physical constants and unit conversions for use with the HHV–TRP
engine and inflation helpers.

Values here are standard and can be cross-checked against PDG or
textbook references. We keep only what we actually need.

We work mostly in reduced Planck units where Mpl = 1, but provide
conversion factors for GeV as well.
"""

import numpy as np

# Reduced Planck mass (Mpl = (8 pi G)^(-1/2)) in GeV
# Commonly quoted ~ 2.435e18 GeV
MPL_REDUCED_GEV = 2.435e18

# Conversion from H in (units of Mpl) to GeV
def H_in_GeV(H_in_Mpl):
    """
    Convert H in units of reduced Planck mass (Mpl=1) to GeV.

    Parameters
    ----------
    H_in_Mpl : float or array
        Hubble scale H in Planck units.

    Returns
    -------
    H_GeV : float or ndarray
        H in GeV.
    """
    H_in_Mpl = np.asarray(H_in_Mpl, dtype=float)
    return H_in_Mpl * MPL_REDUCED_GEV


def V_quarter_in_GeV(V_quarter_in_Mpl_half):
    """
    Convert V^(1/4) from "Mpl^(1/2) * H^(1/2)" units to GeV, assuming
    Mpl is the reduced Planck mass.

    For practical use, if you compute V^(1/4) in Planck units and want GeV:

      V_quarter[GeV] ~ V_quarter[Planck] * Mpl_reduced_GeV

    Parameters
    ----------
    V_quarter_in_Mpl : float or array
        V^(1/4) in Planck units.

    Returns
    -------
    V_quarter_GeV : float or ndarray
        V^(1/4) in GeV.
    """
    V_quarter_in_Mpl = np.asarray(V_quarter_in_Mpl_half, dtype=float)
    return V_quarter_in_Mpl * MPL_REDUCED_GEV
