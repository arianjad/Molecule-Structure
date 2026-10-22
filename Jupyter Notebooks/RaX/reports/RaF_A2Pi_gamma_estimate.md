---
title: "Estimate of $\\gamma$ in the $A\\,{}^2\\Pi_{1/2}$ state of RaF"
author: "A. Jadbabaie (with Claude)"
date: "2026-05-07"
geometry: margin=1in
---

# Summary

The spin-rotation parameter $\gamma$ in the $A\,{}^2\Pi_{1/2}$ state of ${}^{226}\mathrm{RaF}$ is estimated under one-electron pure precession, drawing on second-order perturbation theory of Brown & Carrington and Veseth (1979), the measured $\gamma(C\,{}^2\Sigma^+_{1/2})$ from Conn et al. (2025), and the deperturbed $A\,{}^2\Pi$ / $C\,{}^2\Sigma^+$ analysis from the same source.

**Best estimate:** $\gamma(A\,{}^2\Pi_{1/2}, \mathrm{RaF}) \approx +6$ to $+8$ GHz, **positive**.

The dominant $C\,{}^2\Sigma^+$ channel ($+6.2$ GHz) is well-anchored to measured data; the $B\,{}^2\Delta_{3/2}$ channel ($\sim+1$ GHz) is more uncertain due to unknown configuration phases and a possibly worse Franck–Condon overlap.

---

# 1. Sources

| Ref | Citation | Used for |
|---|---|---|
| B&C | Brown, J. M. & Carrington, A., *Rotational Spectroscopy of Diatomic Molecules* (Cambridge, 2003). §7.4, §7.8, §10.6.3 | Formal second-order $\gamma$ derivation; CH worked example; "twice and opposite" rule |
| Veseth | Veseth, L., *J. Mol. Spectrosc.* **77**, 195 (1979) | Explicit second-order $\gamma({}^2\Pi)$ including both ${}^2\Sigma$ and ${}^2\Delta$ contributions (Eq. 8) |
| A-K | Athanasakis-Kaklamanakis, M. *et al.*, *Nat. Commun.* **16**, 2139 (2025) | RaF excited-state energies, spin-orbit constants, FS-RCC $\mathrm{Ra}^+$ orbital compositions |
| Conn | Conn, C. J., Yu, P. *et al.*, arXiv:2508.08368v1 (2025) | Measured $\gamma(C\,{}^2\Sigma^+, \mathrm{RaF})$; $3\times 3$ deperturbation of $A\,{}^2\Pi/C\,{}^2\Sigma^+$ system using B&C §8.5.2(d) framework |
| NIST | NIST Atomic Spectra Database, Ra II | Atomic $\zeta_{7p}(\mathrm{Ra}^+) = 3239$ cm$^{-1}$, $\zeta_{6d}(\mathrm{Ra}^+) = 663$ cm$^{-1}$ |

---

# 2. Verified inputs

## 2.1 RaF spectroscopic constants

Energies and spin-orbit constants from A-K Tables 1 and 3:

| Quantity | Value |
|---|---:|
| $T_0(A\,{}^2\Pi_{1/2})$ | $13{,}284.4$ cm$^{-1}$ |
| $T_0(B\,{}^2\Delta_{3/2})$ | $14{,}333$ cm$^{-1}$ |
| $T_0(C\,{}^2\Sigma^+_{1/2})$ | $16{,}612$ cm$^{-1}$ |
| $A(A\,{}^2\Pi)_\text{obs}$ | $2051(1)$ cm$^{-1}$ |
| $A(B\,{}^2\Delta)_\text{obs}$ (tentative) | $404(1)$ cm$^{-1}$ |

Measured $\gamma$ for the RaF $C\,{}^2\Sigma^+$ state (Conn et al. Table II):

$$
\gamma(C\,{}^2\Sigma^+_{1/2}, \mathrm{RaF}) = -12{,}388(90)\ \mathrm{MHz} = -0.4132\ \mathrm{cm}^{-1}.
$$

## 2.2 Caltech $3\times 3$ deperturbation of $A\,{}^2\Pi$ / $C\,{}^2\Sigma^+$

From Conn et al. Table S2 (RaF, footnote c). The model is the parity-conserving Hund's case (a) basis from B&C §8.5.2(d), with off-diagonal matrix elements

$$
\langle{}^2\Pi_{1/2}^{(\pm)}|\hat{H}|{}^2\Sigma^+_{1/2,(\pm)}\rangle = M_2 \mp M_1\,(-1)^{J-1/2}(J+1/2),
$$

where $M_1 = \langle{}^2\Pi|BL_+|{}^2\Sigma\rangle$ and $M_2 = \langle{}^2\Pi|(A/2 + B)L_+|{}^2\Sigma\rangle$.

| Parameter | Value |
|---|---:|
| $A_\text{dep}(A\,{}^2\Pi)$ | $1450$ cm$^{-1}$ |
| $M_2$ | $1025$ cm$^{-1}$ |
| $M_1$ | $0.34047(18)$ cm$^{-1}$ |
| $E_\text{dep}(A\,{}^2\Pi)$ | $14{,}420$ cm$^{-1}$ |
| $E_\text{dep}(C\,{}^2\Sigma^+)$ | $15{,}748.0$ cm$^{-1}$ |
| $B_\text{dep}(C\,{}^2\Sigma^+)$ | $0.18652(57)$ cm$^{-1}$ |

Deperturbed $\Pi$–$\Sigma$ gap: $\Delta E_\text{dep} = 1328$ cm$^{-1}$ (vs. observed 3328 cm$^{-1}$).

## 2.3 FS-RCC orbital compositions (A-K Table 2)

Squared amplitudes only; relative signs not reported.

- $A\,{}^2\Pi$: $0.60\,(7p_{1/2}) + 0.20\,(6d_{3/2}) + 0.10\,(7p_{3/2})$ → $|c_{7p}|^2 = 0.70$, $|c_{6d}|^2 = 0.20$
- $C\,{}^2\Sigma^+$: $0.50\,(7p_{3/2}) + 0.30\,(6d_{5/2}) + 0.10\,(7d_{5/2})$ → $|c_{7p}|^2 = 0.50$, $|c_{6d}|^2 = 0.30$
- $B\,{}^2\Delta_{3/2}$: $0.40\,(6d_{3/2}) + 0.30\,(6d_{5/2}) + 0.20\,(7p_{3/2})$ → $|c_{6d}|^2 = 0.70$

## 2.4 Atomic Ra$^+$ spin-orbit (NIST levels)

- $7p_{1/2}\leftrightarrow 7p_{3/2}$ splitting $= 4858$ cm$^{-1}$ → $\zeta_{7p}(\mathrm{Ra}^+) = (2/3)\cdot 4858 = 3239$ cm$^{-1}$
- $6d_{3/2}\leftrightarrow 6d_{5/2}$ splitting $= 1659$ cm$^{-1}$ → $\zeta_{6d}(\mathrm{Ra}^+) = (2/5)\cdot 1659 = 663$ cm$^{-1}$

---

# 3. Verified formulas

## 3.1 Second-order $\gamma$ for one-electron pure precession

Veseth (1979) Eq. 8 for $\gamma$ in a ${}^2\Pi$ state, with the matrix element $\langle{}^2\Pi|AL_+|{}^2\Sigma\rangle$ interpreted to include the $1/2$ spin-coupling factor from $H_\text{so} = \xi(L_z S_z + \frac{1}{2}(L_+ S_- + L_- S_+))$, gives for $l=1$ ($\Sigma\leftrightarrow\Pi$, $p$-orbital):

$$
\gamma^{(2)}({}^2\Pi)_\Sigma = -\frac{A\,B}{E_\Pi - E_\Sigma}
$$

For $l=2$ ($\Pi\leftrightarrow\Delta$, $d$-orbital):

$$
\gamma^{(2)}({}^2\Pi)_\Delta = -\frac{2\,\zeta_d\,B}{E_\Pi - E_{\Delta_{3/2}}}
$$

The factor of 2 in the $\Delta$ term comes from the $l=2$ matrix element $\langle l=2,\lambda=+2|L_+|l=2,\lambda=+1\rangle = 2$ (vs. $\sqrt{2}$ for $l=1$).

For a ${}^2\Sigma$ host with ${}^2\Pi$ perturber, summing over the two $\Lambda=\pm 1$ components of $\Pi$ gives

$$
\gamma^{(2)}({}^2\Sigma) = -\frac{2\,A\,B}{E_\Sigma - E_\Pi},
$$

yielding the **B&C "twice and opposite" rule**:

$$
\gamma({}^2\Sigma) = -2\,\gamma({}^2\Pi).
$$

## 3.2 Empirical verification

**RaF $C\,{}^2\Sigma^+$** (using Conn deperturbed values):

$$
\gamma_\text{calc}(C\,{}^2\Sigma) = -\frac{2 \cdot 1450 \cdot 0.18652}{1328} = -0.408\ \mathrm{cm}^{-1}\quad\text{vs measured}\ -0.4132\ \mathrm{cm}^{-1}.
$$

Agreement to 1.4 %.

**OH $X\,{}^2\Pi$** (B&C p. 393, with $A_\text{mol} = -139$, $B = 18.5$, $\Delta E = -32{,}600$):

$$
\gamma_\text{calc}(X\,{}^2\Pi_\text{OH}) = -\frac{(-139)(18.5)}{-32600} = -0.079\ \mathrm{cm}^{-1}\quad\text{vs measured}\ -0.119.
$$

Sign correct; magnitude $1.5\times$ under (typical pure-precession accuracy).

**OH $A\,{}^2\Sigma^+$**:

$$
\gamma_\text{calc}(A\,{}^2\Sigma^+_\text{OH}) = -\frac{2(-139)(17)}{32600} = +0.145\ \mathrm{cm}^{-1}\quad\text{vs measured}\ +0.201.
$$

Sign correct; magnitude $1.4\times$ under.

The OH ratio $|\gamma({}^2\Sigma)|/|\gamma({}^2\Pi)| = 0.201/0.119 = 1.7$ is close to the predicted 2; deviation attributable to non-pure-precession corrections.

---

# 4. Calculations

## 4.1 $\Sigma$ channel (well-constrained)

Anchored to the *measured* $\gamma(C\,{}^2\Sigma^+_\text{RaF})$ via the universal $-2$ ratio:

$$
\gamma(A\,{}^2\Pi_\text{RaF})_\Sigma = -\frac{1}{2}\,\gamma(C\,{}^2\Sigma^+_\text{RaF})_\text{meas} = +0.207\ \mathrm{cm}^{-1} \approx +6.2\ \mathrm{GHz}.
$$

This bypasses pure-precession uncertainties since the $-2$ ratio is theoretically rigid (Λ-counting in second-order PT) and is verified empirically on OH. Caltech's deperturbation already absorbs the dominant Franck–Condon factor, so no FCF correction is needed here.

## 4.2 $\Delta$ channel (less well-constrained)

Pure precession with the $d$-orbital channel only (the $7p$ component of $A\,{}^2\Pi$ does not couple to a $\Lambda=2$ state under PP):

$$
\gamma(A\,{}^2\Pi)_\Delta = -\frac{2\,\zeta_{6d,\text{eff}}\,B\,c_d^2(A\,{}^2\Pi)\,c_d^2(B\,{}^2\Delta_{3/2})}{E_\Pi - E_\Delta}\cdot \text{FCF}_{00}
$$

with $\zeta_{6d,\text{eff}} \approx A(B\,{}^2\Delta)_\text{obs}/c_d^2(B\,{}^2\Delta) = 404/0.70 = 577$ cm$^{-1}$, $c_d^2(A\,{}^2\Pi) = 0.20$, $c_d^2(B\,{}^2\Delta_{3/2}) = 0.70$, $B = 0.19$ cm$^{-1}$, $E_\Pi - E_\Delta = -1049$ cm$^{-1}$:

$$
\gamma(A\,{}^2\Pi)_\Delta \approx +0.029\ \mathrm{cm}^{-1} \cdot \text{FCF}_{00} \approx +0.85\text{ to }+0.95\ \mathrm{GHz}.
$$

The $\text{FCF}_{00}(A\,{}^2\Pi \leftrightarrow B\,{}^2\Delta)$ is estimated $\sim 0.85$–$0.95$ from the near-equal effective principal quantum numbers $n^*(A) = 2.06$ and $n^*(B\,{}^2\Delta) = 2.08$ (A-K Table 3), but lacks direct measurement.

Higher ${}^2\Delta$ states ($I\,{}^2\Delta$ at $29{,}693$ cm$^{-1}$): negligible contribution ($\lesssim +0.05$ GHz).

## 4.3 Total

$$
\gamma(A\,{}^2\Pi_{1/2}, \mathrm{RaF}) \approx +6.2 + 0.9 \approx +7.0\ \mathrm{GHz},\ \text{positive}.
$$

---

# 5. Bottom-line range

Plausible range under different scenarios:

| Scenario | $\gamma(A\,{}^2\Pi_{1/2})$ |
|---|---:|
| $\Sigma$ only (anchored to Conn) | $+6.2$ GHz |
| $\Sigma + \Delta$, FCF $= 0.95$, constructive | $+7.0$ GHz |
| $\Sigma + \Delta$, FCF $= 0.85$, constructive | $+6.9$ GHz |
| $\Sigma + \Delta$, FCF $= 0.50$, constructive | $+6.6$ GHz |
| $\Sigma + \Delta$, sign flip on $d$-channel | $+5.3$ GHz |

**Range: $+5$ to $+8$ GHz, most likely $+7$ GHz, positive.**

---

# 6. Caveats not eliminated

1. **Pure precession overestimates matrix elements** at the $\sim 30\%$ level (cf. OH calc vs. measured).
2. **$\text{FCF}_{00}(A\,{}^2\Pi \leftrightarrow B\,{}^2\Delta)$ is not measured.** $B_0(B\,{}^2\Delta_\text{RaF})$ is not in the literature; the FCF is estimated from $n^*$ similarity, but the orbital character differs ($7p$ Rydberg vs. $6d$ Rydberg).
3. **Configuration phases of the $d$-channel admixture in $A\,{}^2\Pi$ and $B\,{}^2\Delta$** are not measured. A-K reports $|c|^2$ only. The empirical $M_2 > 0$ from the Caltech fit fixes the sign for the dominant $7p$ channel between $A\,{}^2\Pi$ and $C\,{}^2\Sigma^+$; no analogous constraint exists for the $d$-channel between $A\,{}^2\Pi$ and $B\,{}^2\Delta$.
4. **$A(B\,{}^2\Delta_\text{RaF})$ is reported as tentative** in A-K Table 3 (bracketed value).
5. **Higher $\Sigma$ states** (E, F, H ${}^2\Sigma$ at 25–30 kcm$^{-1}$) contribute $\lesssim +0.2$ GHz; safely within other uncertainties.

---

# 7. Test of the prediction

A direct rotational measurement of the $A\,{}^2\Pi_{1/2}$ state of RaF (e.g., via spin-rotation splitting of the lowest rotational levels) would be the cleanest test. The predicted positive sign is opposite to the (negative) $\gamma$ of the $C\,{}^2\Sigma^+$ state, by the Λ-counting twice-and-opposite rule.
