"""External validation: migrated RaF A0 (boson, ²²⁶Ra¹⁹F) vs the published
Nat. Phys. 2024 N² constants. Numbers transcribed verbatim (with citations) in
docs/superpowers/plans/2026-05-16-n2-r2-formalism-converter-natphys-oracle.md
from Udrescu et al., Nature Physics 2024, DOI 10.1038/s41567-023-02296-w,
Table I (p.7), ²²⁶Ra¹⁹F, PGopher/N² fit.

Acceptance #3: the paper tabulates constants (spectra are figures, no single
N=1 line individually tabulated), so validation is the constants + physical
band-origin comparison; the corrected Origin must equal T_Π1/2,0.

If a comparison FAILs, do NOT loosen TOL to force green — a real disagreement
is a finding to investigate/report (most likely a unit or N²/R² back-out
error), per the plan."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from molecule_parameters import get_molecule_params, c

# --- transcribed from the oracle markdown (Table I, p.7; A²Π₁/₂, ν=0) ---
PAPER_A_Be_N2_CM    = 0.191015      # B′₀(N²)  cm⁻¹   (5)[15]
PAPER_A_p2q_N2_CM   = -0.41071      # p₀(N²)   cm⁻¹   (3)[7]
PAPER_BAND_ORIGIN_CM = 13284.427    # T_Π1/2,0 cm⁻¹   (1)[20]
# Tolerance in cm⁻¹: far inside the paper's 1σ (B′₀ ~1.6e-5, p₀ ~8e-5) — the
# values match exactly post-correction (entered as the paper value; the
# converter + inverse round-trips), so this is a tightness, not a fudge.
TOL_CM = 1e-6

LAM2 = 1                                            # A²Π ⇒ Λ²=1
b = get_molecule_params('RaF', 'A', '0', 'boson')   # stored R² (post-converter)

# invert the converter (formalism.py, B&C Table 7.2) R² → N²:
be_n2  = b['Be'] + 2 * LAM2 * b['D']                       # B row
p2q_n2 = b['p+2q'] - LAM2 * b['p2q_D']                      # generic-X row
origin_phys = b['Origin'] - LAM2 * be_n2 / c + (LAM2 ** 2) * (b['D'] / c)  # G row

be_n2_cm, p2q_n2_cm = be_n2 / c, p2q_n2 / c

print(f"  Be(N²)   converter={be_n2_cm:.6f} cm⁻¹  paper={PAPER_A_Be_N2_CM:.6f}  "
      f"Δ={abs(be_n2_cm-PAPER_A_Be_N2_CM):.2e}")
print(f"  p+2q(N²) converter={p2q_n2_cm:.6f} cm⁻¹  paper={PAPER_A_p2q_N2_CM:.6f}  "
      f"Δ={abs(p2q_n2_cm-PAPER_A_p2q_N2_CM):.2e}")
print(f"  Origin   converter={origin_phys:.6f} cm⁻¹  paper={PAPER_BAND_ORIGIN_CM:.6f}  "
      f"Δ={abs(origin_phys-PAPER_BAND_ORIGIN_CM):.2e}")

ok = (abs(be_n2_cm - PAPER_A_Be_N2_CM) < TOL_CM
      and abs(p2q_n2_cm - PAPER_A_p2q_N2_CM) < TOL_CM
      and abs(origin_phys - PAPER_BAND_ORIGIN_CM) < TOL_CM)
print("OK natphys" if ok else "FAIL natphys")
sys.exit(0 if ok else 1)
