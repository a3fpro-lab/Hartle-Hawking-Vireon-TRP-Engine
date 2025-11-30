"""
hhv_math: core mathematical utilities for the HHV–TRP Engine.

Everything in this package is:
  - Based on standard, referenceable physics/maths; or
  - An explicit, clearly labeled model/ansatz.

Modules:
  inflation : slow-roll inflation relations (H_I, r, A_s, energy scale).
  cmb_params: CMB parameter constants (Planck results) used for calibration.
  units     : basic Planck/GeV conversions and physical constants.
"""
from . import inflation
from . import cmb_params
from . import units

__all__ = ["inflation", "cmb_params", "units"]
