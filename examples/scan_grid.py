"""
examples/scan_grid.py

Scan a grid in (Ne, H_I) and compute the maximal allowed |eps|
according to the TRP bound, using Planck-like pivot calibration.

Output:
  - CSV file 'trp_eps_max_grid.csv' with columns:
      Ne, HI, eps_max

Requirements:
  - numpy
"""

import csv
import numpy as np

from hhv_math import inflation, cmb_params
from trp_engine import TRPEngine


def main():
    # --- 1. Planck-like pivot parameters ---

    As = cmb_params.A_s_planck_2018
    r_star = cmb_params.r_upper_approx_0p005
    Ne_star = cmb_params.N_e_pivot

    # H_I at pivot from slow-roll relation
    HI_star = inflation.H_from_r_As(r_star, As, Mpl=1.0)

    # --- 2. Instantiate TRPEngine calibrated at pivot ---

    eps_star = 0.02  # from Planck-like statistical isotropy bound

    engine = TRPEngine(
        S0=1e5,
        T_min=3.0,
        k_A=1e3,
        sigma=1.0,
        Ne_star=Ne_star,
        HI_star=HI_star,
        eps_star=eps_star,
    )

    print(f"Calibrated mu ≈ {engine.mu:.3e}")
    print(f"Pivot: Ne*={Ne_star:.1f}, H_I*={HI_star:.3e}, eps*={eps_star:.3f}")

    # --- 3. Define scan grid in (Ne, H_I) ---

    # Ne from 50 to 70 e-folds
    Ne_vals = np.linspace(50.0, 70.0, 41)  # step 0.5

    # H_I from ~0.3 * H_I* to ~3 * H_I* (log spaced)
    HI_factors = np.logspace(-0.5, 0.5, 31)  # ~0.316..3.162
    HI_vals = HI_star * HI_factors

    # --- 4. Compute eps_max on the grid ---

    rows = []
    for Ne in Ne_vals:
        for HI in HI_vals:
            eps_max = float(engine.epsilon_max(Ne, HI))
            rows.append((Ne, HI, eps_max))

    # --- 5. Write CSV file ---

    out_path = "trp_eps_max_grid.csv"
    with open(out_path, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Ne", "H_I_PlanckUnits", "eps_max"])
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {out_path}")
    print("Columns: Ne, H_I (Planck units), eps_max")


if __name__ == "__main__":
    main()
