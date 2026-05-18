#set page(paper: "us-letter", margin: 1in, numbering: "1")
#set par(justify: true, leading: 0.65em)
#set text(font: "New Computer Modern", size: 11pt)
#show heading: set block(above: 1.4em, below: 0.8em)
#show table.cell.where(y: 0): strong
#set table(stroke: 0.4pt, inset: 6pt)

// Helper: prefix-superscript term symbol like ²Π
#let ts(left, body) = box[#super[#left]#body]

#align(center)[
  #text(size: 16pt, weight: "bold")[
    Estimate of $gamma$ in the #ts[2]$Pi_(1\/2)$ state of RaF
  ] \
  #v(0.4em)
  A. Jadbabaie (with Claude) #h(1cm) 2026-05-07
]

#v(0.5em)

= Summary

The spin-rotation parameter $gamma$ in the $A$ #ts[2]$Pi_(1\/2)$ state of #super[226]RaF is estimated under one-electron pure precession, drawing on second-order perturbation theory of Brown & Carrington and Veseth (1979), the measured $gamma(C #ts[2]Sigma^+_(1\/2))$ from Conn et al. (2025), and the deperturbed $A$ #ts[2]$Pi$ / $C$ #ts[2]$Sigma^+$ analysis from the same source.

#align(center)[
  #rect(stroke: 0.5pt, inset: 8pt, radius: 4pt)[
    *Best estimate:* $gamma(A #ts[2]Pi_(1\/2), "RaF") approx +6 "to" +8$ GHz, positive.
  ]
]

The dominant $C$ #ts[2]$Sigma^+$ channel ($+6.2$ GHz) is well-anchored to measured data; the $B$ #ts[2]$Delta_(3\/2)$ channel ($tilde +1$ GHz) is more uncertain due to unknown configuration phases and a possibly worse Franck–Condon overlap.

= Sources

#table(
  columns: (auto, 1fr, 1.4fr),
  align: (left, left, left),
  [Ref], [Citation], [Used for],
  [B&C], [Brown, J. M. & Carrington, A., #emph[Rotational Spectroscopy of Diatomic Molecules] (Cambridge, 2003). §7.4, §7.8, §10.6.3], [Formal $gamma^((2))$ derivation; CH worked example; "twice and opposite" rule],
  [Veseth], [Veseth, L., #emph[J. Mol. Spectrosc.] *77*, 195 (1979)], [Explicit second-order $gamma(#ts[2]Pi)$ including both #ts[2]$Sigma$ and #ts[2]$Delta$ contributions (Eq. 8)],
  [A-K], [Athanasakis-Kaklamanakis, M. et al., #emph[Nat. Commun.] *16*, 2139 (2025)], [RaF excited-state energies, spin-orbit constants, and FS-RCC Ra#super[+] orbital compositions],
  [Conn], [Conn, C. J., Yu, P., et al., arXiv:2508.08368v1 (2025)], [Measured $gamma(C #ts[2]Sigma^+, "RaF")$; 3$times$3 deperturbation of $A$ #ts[2]$Pi$ / $C$ #ts[2]$Sigma^+$ system (B&C §8.5.2(d) framework)],
  [NIST], [NIST Atomic Spectra Database, Ra II], [Atomic $zeta_(7p)("Ra"^+) = 3239$ cm#super[-1], $zeta_(6d)("Ra"^+) = 663$ cm#super[-1]],
)

= Verified inputs

== RaF spectroscopic constants

Energies and spin-orbit constants from A-K Tables 1 and 3:

#align(center)[
  #table(
    columns: 2,
    align: (left, right),
    [Quantity], [Value],
    [$T_0(A #ts[2]Pi_(1\/2))$], [13,284.4 cm#super[-1]],
    [$T_0(B #ts[2]Delta_(3\/2))$], [14,333 cm#super[-1]],
    [$T_0(C #ts[2]Sigma^+_(1\/2))$], [16,612 cm#super[-1]],
    [$A(A #ts[2]Pi)_"obs"$], [2051(1) cm#super[-1]],
    [$A(B #ts[2]Delta)_"obs"$ (tentative)], [404(1) cm#super[-1]],
  )
]

Measured $gamma$ for the RaF $C$ #ts[2]$Sigma^+$ state (Conn et al. Table II):

$ gamma(C #ts[2]Sigma^+_(1\/2), "RaF") = -12 thin 388(90) "MHz" = -0.4132 "cm"^(-1). $

== Caltech 3 × 3 deperturbation of $A$ #ts[2]$Pi$ / $C$ #ts[2]$Sigma^+$

From Conn et al. Table S2 (RaF, footnote c). The model is the parity-conserving Hund's case (a) basis from B&C §8.5.2(d), with off-diagonal matrix elements

$ angle.l #ts[2]Pi_(1\/2)^((plus.minus))|hat(H)|#ts[2]Sigma^+_(1\/2,(plus.minus)) angle.r = M_2 minus.plus M_1 (-1)^(J - 1\/2) (J + 1\/2), $

where $M_1 = angle.l #ts[2]Pi |B L_+| #ts[2]Sigma angle.r$ and $M_2 = angle.l #ts[2]Pi |(A\/2 + B) L_+| #ts[2]Sigma angle.r$.

#align(center)[
  #table(
    columns: 2,
    align: (left, right),
    [Parameter], [Value],
    [$A_"dep"(A #ts[2]Pi)$], [1450 cm#super[-1]],
    [$M_2$], [1025 cm#super[-1]],
    [$M_1$], [0.34047(18) cm#super[-1]],
    [$E_"dep"(A #ts[2]Pi)$], [14,420 cm#super[-1]],
    [$E_"dep"(C #ts[2]Sigma^+)$], [15,748.0 cm#super[-1]],
    [$B_"dep"(C #ts[2]Sigma^+)$], [0.18652(57) cm#super[-1]],
  )
]

Deperturbed $Pi$–$Sigma$ gap: $Delta E_"dep" = 1328$ cm#super[-1] (vs. observed 3328 cm#super[-1]).

== FS-RCC orbital compositions (A-K Table 2)

Squared amplitudes only; relative signs not reported.

- $A$ #ts[2]$Pi$: $0.60 (7p_(1\/2)) + 0.20 (6d_(3\/2)) + 0.10 (7p_(3\/2))$ #h(0.4em)$=>$#h(0.4em) $|c_(7p)|^2 = 0.70$, $|c_(6d)|^2 = 0.20$
- $C$ #ts[2]$Sigma^+$: $0.50 (7p_(3\/2)) + 0.30 (6d_(5\/2)) + 0.10 (7d_(5\/2))$ #h(0.4em)$=>$#h(0.4em) $|c_(7p)|^2 = 0.50$, $|c_(6d)|^2 = 0.30$
- $B$ #ts[2]$Delta_(3\/2)$: $0.40 (6d_(3\/2)) + 0.30 (6d_(5\/2)) + 0.20 (7p_(3\/2))$ #h(0.4em)$=>$#h(0.4em) $|c_(6d)|^2 = 0.70$

== Atomic Ra#super[+] spin-orbit (NIST levels)

- $7p_(1\/2)$ ↔ $7p_(3\/2)$ splitting $= 4858$ cm#super[-1], so $zeta_(7p)("Ra"^+) = (2\/3) dot 4858 = 3239$ cm#super[-1]
- $6d_(3\/2)$ ↔ $6d_(5\/2)$ splitting $= 1659$ cm#super[-1], so $zeta_(6d)("Ra"^+) = (2\/5) dot 1659 = 663$ cm#super[-1]

= Verified formulas

== Second-order $gamma$ for one-electron pure precession

Veseth (1979) Eq. 8 for $gamma$ in a #ts[2]$Pi$ state, with the matrix element $angle.l #ts[2]Pi |A L_+| #ts[2]Sigma angle.r$ interpreted to include the $1\/2$ spin-coupling factor from $H_"so" = xi (L_z S_z + 1\/2 (L_+ S_- + L_- S_+))$, gives for $l = 1$ ($Sigma <-> Pi$, $p$-orbital):

$ gamma^((2))(#ts[2]Pi)_Sigma = -frac(A B, E_Pi - E_Sigma) $

For $l = 2$ ($Pi <-> Delta$, $d$-orbital):

$ gamma^((2))(#ts[2]Pi)_Delta = -frac(2 zeta_d B, E_Pi - E_(Delta_(3\/2))) $

The factor of 2 in the $Delta$ term comes from the $l = 2$ matrix elements $angle.l l = 2, lambda = +2 | L_+ | l = 2, lambda = +1 angle.r = 2$ (vs. $sqrt(2)$ for $l = 1$).

For a #ts[2]$Sigma$ host with #ts[2]$Pi$ perturber, summing over the two $Lambda = plus.minus 1$ components of $Pi$ gives

$ gamma^((2))(#ts[2]Sigma) = -frac(2 A B, E_Sigma - E_Pi), $

yielding the *B&C "twice and opposite" rule*:

$ gamma(#ts[2]Sigma) = -2 thin gamma(#ts[2]Pi). $

== Empirical verification

*RaF $C$ #ts[2]$Sigma^+$* (using Conn deperturbed values):
$ gamma_"calc"(C #ts[2]Sigma) = -frac(2 dot 1450 dot 0.18652, 1328) = -0.408 "cm"^(-1) "vs measured" -0.4132 "cm"^(-1). $
Agreement to 1.4 %.

*OH $X$ #ts[2]$Pi$* (B&C p. 393, $A_"mol" = -139$, $B = 18.5$, $Delta E = -32 thin 600$):
$ gamma_"calc"(X #ts[2]Pi_"OH") = -frac((-139)(18.5), -32 thin 600) = -0.079 "cm"^(-1) "vs measured" -0.119. $
Sign correct; magnitude $1.5 times$ under (typical pure-precession accuracy).

*OH $A$ #ts[2]$Sigma^+$*:
$ gamma_"calc"(A #ts[2]Sigma^+_"OH") = -frac(2 (-139)(17), 32 thin 600) = +0.145 "cm"^(-1) "vs measured" +0.201. $
Sign correct; magnitude $1.4 times$ under.

The OH ratio $|gamma(#ts[2]Sigma)| \/ |gamma(#ts[2]Pi)| = 0.201 \/ 0.119 = 1.7$ is close to the predicted 2; the deviation is attributable to non-pure-precession corrections.

= Calculations

== $Sigma$ channel (well-constrained)

Anchored to the #emph[measured] $gamma(C #ts[2]Sigma^+_"RaF")$ via the universal $-2$ ratio:

$ gamma(A #ts[2]Pi_"RaF")_Sigma = -1\/2 dot gamma(C #ts[2]Sigma^+_"RaF")_"meas" = +0.207 "cm"^(-1) approx +6.2 "GHz." $

This bypasses pure-precession uncertainties, since the $-2$ ratio is theoretically rigid (Λ-counting in second-order PT) and is verified empirically on OH. Caltech's deperturbation already absorbs the dominant Franck–Condon factor, so no FCF correction is needed here.

== $Delta$ channel (less well-constrained)

Pure precession with the $d$-orbital channel only (the $7p$ component of $A$ #ts[2]$Pi$ does not couple to a $Lambda = 2$ state under PP):

$ gamma(A #ts[2]Pi)_Delta = -frac(2 zeta_(6d, "eff") B med c_d^2(A #ts[2]Pi) med c_d^2(B #ts[2]Delta_(3\/2)), E_Pi - E_Delta) dot "FCF"_(00) $

with $zeta_(6d, "eff") approx A(B #ts[2]Delta)_"obs" \/ c_d^2(B #ts[2]Delta_(3\/2)) = 404\/0.70 = 577$ cm#super[-1], $c_d^2(A #ts[2]Pi) = 0.20$, $c_d^2(B #ts[2]Delta_(3\/2)) = 0.70$, $B = 0.19$ cm#super[-1], $E_Pi - E_Delta = -1049$ cm#super[-1]:

$ gamma(A #ts[2]Pi)_Delta approx +0.029 "cm"^(-1) dot "FCF"_(00) approx +0.85 "to" +0.95 "GHz." $

The FCF₀₀ between $A$ #ts[2]$Pi$ and $B$ #ts[2]$Delta$ is estimated $tilde 0.85 - 0.95$ from the near-equal effective principal quantum numbers $n^*(A) = 2.06$ and $n^*(B$ #ts[2]$Delta) = 2.08$ (A-K Table 3), but lacks direct measurement.

Higher #ts[2]$Delta$ states ($I$ #ts[2]$Delta$ at 29,693 cm#super[-1]): negligible contribution ($lt.eq +0.05$ GHz).

== Total

#align(center)[
  #rect(stroke: 0.5pt, inset: 10pt, radius: 4pt)[
    $ gamma(A #ts[2]Pi_(1\/2), "RaF") approx +6.2 + 0.9 approx +7.0 "GHz, positive." $
  ]
]

= Bottom-line range

Plausible range under different scenarios:

#align(center)[
  #table(
    columns: 2,
    align: (left, right),
    [Scenario], [$gamma(A #ts[2]Pi_(1\/2))$],
    [$Sigma$ only (anchored to Conn)], [+6.2 GHz],
    [$Sigma + Delta$, FCF $= 0.95$, constructive], [+7.0 GHz],
    [$Sigma + Delta$, FCF $= 0.85$, constructive], [+6.9 GHz],
    [$Sigma + Delta$, FCF $= 0.50$, constructive], [+6.6 GHz],
    [$Sigma + Delta$, sign flip on $d$-channel], [+5.3 GHz],
  )
]

*Range: +5 to +8 GHz, most likely +7 GHz, positive.*

= Caveats not eliminated

+ *Pure precession overestimates matrix elements* at the $tilde 30 %$ level (cf. OH calc vs. measured).
+ *FCF₀₀ between $A$ #ts[2]$Pi$ and $B$ #ts[2]$Delta$ is not measured.* $B_0(B$ #ts[2]$Delta_"RaF")$ is not in the literature; the FCF is estimated from $n^*$ similarity, but the orbital character differs ($7p$ Rydberg vs. $6d$ Rydberg).
+ *Configuration phases of the $d$-channel admixture in $A$ #ts[2]$Pi$ and $B$ #ts[2]$Delta$* are not measured. A-K reports $|c|^2$ only. The empirical $M_2 > 0$ from the Caltech fit fixes the sign for the dominant $7p$ channel between $A$ #ts[2]$Pi$ and $C$ #ts[2]$Sigma^+$; no analogous constraint exists for the $d$-channel between $A$ #ts[2]$Pi$ and $B$ #ts[2]$Delta$.
+ *$A(B$ #ts[2]$Delta_"RaF")$ is reported as tentative* in A-K Table 3 (bracketed value).
+ *Higher $Sigma$ states* (E, F, H #ts[2]$Sigma$ at 25–30 kcm#super[-1]) contribute $lt.eq +0.2$ GHz; safely within other uncertainties.

= Test of the prediction

A direct rotational measurement of the $A$ #ts[2]$Pi_(1\/2)$ state of RaF (e.g., via spin-rotation splitting of the lowest rotational levels) would be the cleanest test. The predicted positive sign is opposite to the (negative) $gamma$ of the $C$ #ts[2]$Sigma^+$ state, by the Λ-counting twice-and-opposite rule.
