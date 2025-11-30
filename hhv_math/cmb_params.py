"""
cmb_params.py

Numerical constants from CMB observations (Planck 2018 and combined
Planck + BICEP/Keck analyses), used as reference values in the HHV–TRP
engine.

These are not "fitted" by this repo; they are imported as external
inputs. Values are approximate and should be updated if newer
constraints are adopted.

References (for documentation level; not enforced in code):

- Planck 2018 results. VI. Cosmological parameters.
- Planck 2018 results. X. Constraints on inflation.
- BICEP/Keck Array XII / XIII etc. combined limits on r.
"""

# Scalar amplitude A_s at k=0.05 Mpc^-1 (Planck 2018 best fit)
A_s_planck_2018 = 2.1e-9

# Example upper bounds on r at CMB scales (approximate, 95% CL).
# These should be treated as representative, not exact.
r_upper_approx_0p005 = 0.03  # r at k ~ 0.005 Mpc^-1

# Pivot Ne for our use (typical value for CMB scales)
N_e_pivot = 60.0

# For convenience, store a "canonical" pivot tuple
PIVOT = {
    "As": A_s_planck_2018,
    "r_upper": r_upper_approx_0p005,
    "Ne": N_e_pivot,
}
