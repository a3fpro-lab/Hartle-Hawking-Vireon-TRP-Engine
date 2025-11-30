# Impact of the TRP Anisotropy Bound on Early-Universe Models

This document explains how the Hartle–Hawking–Vireon TRP engine in this
repository constrains early-universe scenarios, using only standard
slow-roll + CMB inputs and the explicit TRP ansatz implemented in
`trp_engine.py`.

We keep all assumptions and formulas explicit so that every prediction is
reproducible from the code in this repo.

---

## 1. TRP structure

The core TRP law is multiplicative:
\[
T = R \times P.
\]

In the engine, we define:

- **Reality bandwidth** (from exit-geometry entropy):
  \[
  A(N_e, H_I) = k_A\,\frac{e^{2N_e}}{H_I^2},
  \qquad
  S_{\rm geom}(N_e,H_I) = \frac{A(N_e,H_I)}{4},
  \]
  \[
  R(N_e,H_I) = \log\!\left(\frac{S_{\rm geom}(N_e,H_I)}{S_0}\right),
  \]
  where \(S_0\) is a reference entropy scale and \(k_A\) is a geometric
  prefactor.

- **Curvature complexity** (quadratic ansatz):
  \[
  C(\epsilon) = \frac{\epsilon^2}{2\sigma^2},
  \]
  where \(\epsilon\) parametrizes exit anisotropy and \(\sigma\) is a
  scale parameter.

- **Perception gain**:
  \[
  P(\epsilon) = \exp\big(-\mu\,C(\epsilon)\big).
  \]

The total TRP capacity is:
\[
T(N_e,H_I,\epsilon) = R(N_e,H_I)\,\exp\big(-\mu\,C(\epsilon)\big).
\]

The viability condition is:
\[
T(N_e,H_I,\epsilon) \;\ge\; T_{\min}.
\]

All of these are implemented directly in `trp_engine.py` as:

- `R(Ne, HI)`
- `C(eps)`
- `P(eps)`
- `T(Ne, HI, eps) = R * P`

---

## 2. Calibration: fixing \(\mu\) from CMB anisotropy

There is a one-parameter freedom in how much “difficulty” we assign to
geometry vs perception: any rescaling
\[
R \rightarrow a\,R, 
\qquad
P \rightarrow \frac{1}{a} P
\]
leaves \(T = R P\) invariant. To make the framework predictive we fix
that freedom with a **pivot calibration** at a Planck-like point:
\[
(N_e^\*, H_I^\*, \epsilon^\*).
\]

In code, the pivot is defined in `hhv_math` via:

- `cmb_params.A_s_planck_2018` (scalar amplitude),
- `cmb_params.r_upper_approx_0p005` (an upper bound on \(r\)),
- `cmb_params.N_e_pivot` (a representative pivot \(N_e^\*\)),

and we obtain \(H_I^\*\) using:
```python
HI_star = inflation.H_from_r_As(r_star, As, Mpl=1.0)

We choose a reference anisotropy level (\epsilon^*) at that pivot
(e.g. (\epsilon^* = 0.02)) and define (\mu) by requiring:
[
T(N_e^*,H_I^*,\epsilon^*) = T_{\min}.
]

This condition is implemented in _calibrate_mu() as:
[
R_* ,\exp(-\mu C_*) = T_{\min}
\quad\Rightarrow\quad
\mu = \frac{1}{C_*},\log!\left(\frac{R_*}{T_{\min}}\right),
]
where
[
R_* = R(N_e^*,H_I^*), \qquad
C_* = C(\epsilon^*).
]

Once (\mu) is fixed, there is no further freedom: every
(|\epsilon|_{\max}(N_e,H_I)) that the engine reports is a direct
consequence of this calibration plus the TRP ansatz.

⸻

3. Complexity bound and anisotropy bound

From the condition (T \ge T_{\min}) we obtain a bound on the allowed
complexity:
[
R(N_e,H_I),\exp\big(-\mu,C(\epsilon)\big) \ge T_{\min}.
]

Rearranging:
[
\exp\big(-\mu,C(\epsilon)\big) \ge \frac{T_{\min}}{R(N_e,H_I)}
\quad\Rightarrow\quad
-\mu,C(\epsilon) \ge \log\frac{T_{\min}}{R(N_e,H_I)}.
]

For regions where (R(N_e,H_I) > T_{\min}), this is equivalent to:
[
C(\epsilon) \le \frac{1}{\mu},\log!\left(\frac{R(N_e,H_I)}{T_{\min}}\right)
;:=; C_{\max}(N_e,H_I).
]

With the quadratic ansatz (C(\epsilon) = \epsilon^2/(2\sigma^2)), we
get:
[
\frac{\epsilon^2}{2\sigma^2} \le C_{\max}(N_e,H_I)
\quad\Rightarrow\quad
|\epsilon| \le \sigma,\sqrt{2,C_{\max}(N_e,H_I)}.
]

In code, this is exactly what epsilon_max(Ne, HI) computes:
	•	C_max(Ne, HI) implements the formula above,
	•	epsilon_max(Ne, HI) = sigma * sqrt(2 * C_max).

Thus, for any early-universe model that predicts an exit anisotropy
(\epsilon_{\rm model}(N_e,H_I)), HHV-TRP declares it viable iff:
[
|\epsilon_{\rm model}(N_e,H_I)| \le \epsilon_{\max}(N_e,H_I),
]
otherwise the branch is non-viable (TRP too low).

⸻

4. Impact on early-universe models

Because (R(N_e,H_I)) grows with exit area (and hence with (N_e) and
lower (H_I)), the TRP bound has clear qualitative consequences:

4.1 High-scale, short inflation is tightly constrained
	•	Large (H_I) (high-scale inflation) and modest (N_e) imply smaller
exit area:
[
A(N_e,H_I) \propto \frac{e^{2N_e}}{H_I^2}.
]
	•	This yields smaller (S_{\rm geom}) and therefore smaller (R).
	•	For smaller (R) at fixed (T_{\min}), the allowed complexity:
[
C_{\max}(N_e,H_I)
= \frac{1}{\mu},\log!\left(\frac{R(N_e,H_I)}{T_{\min}}\right)
]
decreases, which tightens the bound (|\epsilon|_{\max}).

Interpretation:
High-scale, minimally long inflation (large (r), near-minimal
(N_e)) is only TRP-viable if it yields very isotropic exit
geometries. Models with significant primordial anisotropy at high
(H_I) are disfavored/ruled out by HHV-TRP.

4.2 Low-scale, long inflation is more tolerant
	•	For lower (H_I) or larger (N_e), the area grows and so does
(S_{\rm geom}), hence (R(N_e,H_I)) increases.
	•	Larger (R) at fixed (T_{\min}) yields a larger (C_{\max}), and
therefore a larger (|\epsilon|_{\max}).

Interpretation:
Low-scale, long-duration inflation can tolerate more pre-inflation
geometric “mess” while still having enough TRP to support recursive
observers.

4.3 Anisotropic / Bianchi models get a sharp filter

If a specific model (e.g. Bianchi IX seeds, vector hair, etc.) provides
a prediction (\epsilon_{\rm model}(N_e,H_I)), you can test it against

HHV-TRP in a few lines:
from trp_engine import TRPEngine
from hhv_math import inflation, cmb_params

# 1. Build engine from a Planck-like pivot
As = cmb_params.A_s_planck_2018
r_star = cmb_params.r_upper_approx_0p005
Ne_star = cmb_params.N_e_pivot
HI_star = inflation.H_from_r_As(r_star, As, Mpl=1.0)

engine = TRPEngine(
    S0=1e5,
    T_min=3.0,
    k_A=1e3,
    sigma=1.0,
    Ne_star=Ne_star,
    HI_star=HI_star,
    eps_star=0.02,
)

# 2. Model prediction at some (Ne_model, HI_model)
Ne_model = 60.0
HI_model = HI_star
eps_model = 0.01  # example: model's predicted exit anisotropy

# 3. TRP anisotropy bound
eps_max = engine.epsilon_max(Ne_model, HI_model)
viable = abs(eps_model) <= eps_max

print("eps_model =", eps_model)
print("eps_max   =", float(eps_max))
print("TRP-viable:", bool(viable))

•	If viable is False, that branch is excluded by HHV-TRP, even
if it is allowed by the underlying inflationary dynamics.
	•	If viable is True, the branch survives the TRP selection factor.

⸻

5. What is and is not claimed

What this repo does claim:
	•	A fully explicit, test-backed implementation of a TRP selection
factor on top of the Hartle–Hawking no-boundary framework.
	•	A clear mapping from standard slow-roll + CMB inputs to an anisotropy
bound (|\epsilon|_{\max}(N_e,H_I)).
	•	A reproducible way to filter early-universe models based on their
predicted exit anisotropy.

What this repo does not claim (yet):
	•	A derivation of (\mu) or the quadratic complexity ansatz from a
microscopic quantum-gravity calculation.
	•	A stronger bound on inflationary parameters than existing CMB data
alone; instead, TRP adds an observer-centric selection layer that
refines which branches are viable within a given theory.

In other words, HHV-TRP is designed as a bridge:
given any anisotropic Hartle–Hawking (or tunneling) calculation, you can
plug its exit geometry into this engine and obtain a transparent, fully
documented bound on the maximal allowed anisotropy consistent with a
minimal TRP requirement.


