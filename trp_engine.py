

```python
import numpy as np


class TRPEngine:
    """
    Hartle–Hawking–Vireon TRP engine with explicit formulas and data-backed calibration.

    Definitions (Planck units, 4 G ħ = 1):

      A(Ne, HI)    = k_A * exp(2 Ne) / HI^2
      S_geom       = A / 4
      R(Ne, HI)    = log(S_geom / S0)
      C(eps)       = eps^2 / (2 sigma^2)
      T(Ne,HI,eps) = R * exp(-mu * C)

    TRP viability:
      T >= T_min  <=>  C <= C_max(Ne,HI)

      C_max(Ne,HI) = (1/mu) * log(R / T_min),
                     R = log(S_geom / S0),
                     defined only for R > T_min.

      |eps|_max(Ne,HI) = sigma * sqrt(2 * C_max(Ne,HI)).

    Calibration:
      mu is fixed by requiring
          T(Ne_star, HI_star, eps_star) = T_min
      at a pivot (Ne_star, HI_star) consistent with CMB (e.g. r and A_s),
      and eps_star consistent with Planck's statistical isotropy bounds.
    """

    def __init__(self,
                 S0=1e5,
                 T_min=3.0,
                 k_A=1e3,
                 sigma=1.0,
                 Ne_star=60.0,
                 HI_star=1.8e-5,
                 eps_star=0.02):
        """
        Parameters
        ----------
        S0 : float
            Reference entropy scale for the minimal recursion engine.
        T_min : float
            Minimum TRP threshold for viability.
        k_A : float
            Geometric prefactor in A(Ne, HI).
        sigma : float
            Normalization scale in C(eps) = eps^2 / (2 sigma^2).
        Ne_star : float
            Pivot e-fold number for calibration (e.g. 60).
        HI_star : float
            Pivot Hubble scale for calibration (Planck units).
        eps_star : float
            Pivot anisotropy amplitude used to calibrate mu.
        """
        self.S0 = float(S0)
        self.T_min = float(T_min)
        self.k_A = float(k_A)
        self.sigma = float(sigma)

        self.Ne_star = float(Ne_star)
        self.HI_star = float(HI_star)
        self.eps_star = float(eps_star)

        # Calibrate mu so that T(Ne_star, HI_star, eps_star) = T_min
        self.mu = self._calibrate_mu()

    # ---------- Geometry and entropy ----------

    def area(self, Ne, HI):
        """
        Exit surface area A(Ne, HI) = k_A * exp(2 Ne) / HI^2.
        Ne, HI can be scalars or arrays.
        """
        Ne = np.array(Ne, dtype=float)
        HI = np.array(HI, dtype=float)
        return self.k_A * np.exp(2.0 * Ne) / (HI ** 2)

    def S_geom(self, Ne, HI):
        """
        Geometric entropy S_geom = A / 4.
        """
        return 0.25 * self.area(Ne, HI)

    def R(self, Ne, HI):
        """
        Reality bandwidth R = log(S_geom / S0).
        """
        Sg = self.S_geom(Ne, HI)
        return np.log(Sg / self.S0)

    # ---------- Calibration of mu ----------

    def _calibrate_mu(self):
        """
        Solve for mu from the condition
            T(Ne_star, HI_star, eps_star) = T_min,
        i.e.
            R_star * exp(-mu * C_star) = T_min,
        with C_star = eps_star^2 / (2 sigma^2).
        """
        R_star = float(self.R(self.Ne_star, self.HI_star))
        C_star = (self.eps_star ** 2) / (2.0 * self.sigma ** 2)

        if R_star <= self.T_min:
            raise ValueError(
                f"Calibration invalid: R_star={R_star:.3g} <= T_min={self.T_min:.3g}."
            )
        if C_star <= 0.0:
            raise ValueError(
                f"Calibration invalid: eps_star^2/(2 sigma^2) must be > 0."
            )

        mu = (1.0 / C_star) * np.log(R_star / self.T_min)
        return float(mu)

    # ---------- TRP bound and epsilon_max ----------

    def C_max(self, Ne, HI):
        """
        TRP complexity bound:

          C_max(Ne,HI) = (1/mu) * log(R / T_min),
          R = log(S_geom / S0).

        Returns
        -------
        Cmax : ndarray or float
            Maximal allowed curvature complexity. 0 if R <= T_min.
        """
        R_val = np.array(self.R(Ne, HI), dtype=float)
        Cmax = np.zeros_like(R_val)

        mask = R_val > self.T_min
        Cmax[mask] = (1.0 / self.mu) * np.log(R_val[mask] / self.T_min)
        Cmax = np.where(Cmax > 0.0, Cmax, 0.0)
        return Cmax

    def epsilon_max(self, Ne, HI):
        """
        Maximal allowed |eps| at given (Ne, HI):
          |eps|_max = sigma * sqrt(2 C_max(Ne,HI)).

        Returns
        -------
        eps_max : ndarray or float
        """
        Cmax = np.array(self.C_max(Ne, HI), dtype=float)
        with np.errstate(invalid="ignore"):
            eps_max = self.sigma * np.sqrt(2.0 * Cmax)
            eps_max[~np.isfinite(eps_max)] = 0.0
        return eps_max

    # ---------- Full TRP and viability ----------

    def T(self, Ne, HI, eps):
        """
        Total Recursive Processing:

          T(Ne,HI,eps) = R(Ne,HI) * exp(-mu * C(eps)),
          C(eps) = eps^2 / (2 sigma^2).
        """
        R_val = self.R(Ne, HI)
        eps_arr = np.array(eps, dtype=float)
        C_val = (eps_arr ** 2) / (2.0 * self.sigma ** 2)
        return R_val * np.exp(-self.mu * C_val)

    def is_viable(self, Ne, HI, eps):
        """
        Check TRP viability: True if T(Ne,HI,eps) >= T_min.
        """
        return self.T(Ne, HI, eps) >= self.T_min


if __name__ == "__main__":
    # Example: instantiate with default Planck-like calibration
    engine = TRPEngine(
        S0=1e5,
        T_min=3.0,
        k_A=1e3,
        sigma=1.0,
        Ne_star=60.0,
        HI_star=1.8e-5,
        eps_star=0.02,
    )

    print(f"Calibrated mu ≈ {engine.mu:.3g}")

    # Pivot check
    eps_max_star = engine.epsilon_max(engine.Ne_star, engine.HI_star)
    print(f"|eps|_max(Ne*={engine.Ne_star}, HI*={engine.HI_star:.2e}) = {eps_max_star:.5f}")

    # Some nearby sample points
    for (Ne, HI) in [(60.0, 1.0e-5), (55.0, 1.0e-5), (60.0, 5.0e-6), (65.0, 1.0e-5)]:
        eps_max = engine.epsilon_max(Ne, HI)
        print(f"Ne={Ne:5.1f}, HI={HI:7.2e} -> |eps|_max = {eps_max:.5f}")
