# RaF Effective-Hamiltonian Constants — Reference Document (Design / Spec)

**Date:** 2026-05-19
**Author:** Arian Jadbabaie (with Claude)
**Status:** Approved design, pending implementation plan
**Deliverable location:** `~/obsidian-vault/04-resources/` (LaTeX + HTML)
**This spec (scaffolding):** `Code/Molecule-Structure/docs/superpowers/specs/2026-05-19-raf-hamiltonian-constants-reference-design.md`

---

## 1. Goal

A single curated reference for the effective-Hamiltonian constants used to model **²²⁶Ra¹⁹F** in the Molecule-Structure engine, in two renderings:

- `RaF-effective-hamiltonian-constants.tex` — self-contained, Overleaf-compilable (`pdflatex`).
- `RaF-effective-hamiltonian-constants.html` — self-contained, MathJax-3 (SVG); authored to be browser-rendered then copy-pasted into OneNote with equations as crisp images.

For each constant the document gives, in order: (1) the Hamiltonian form, (2) the Brown & Carrington pages/equations that define and reduce it, (3) the important effective-Hamiltonian contributions that fold into the fitted constant, (4) conventions around phase/prefactor/interpretation as a three-column **B&C / N²-paper / PGopher** comparison with explicit transformation maps.

## 2. Hard invariants (non-negotiable)

These come directly from the user and override convenience:

1. **Every claim is inline-cited.** No uncited sentence ships. Citation taxonomy:

   | Tag | Source |
   |---|---|
   | `[B&C eq. X.XX, p.NNN]` | Brown & Carrington, *Rotational Spectroscopy of Diatomic Molecules*, transcribed verbatim via pdf-mcp from the Zotero PDF (`storage/CKZKCGXY`) |
   | `[Udrescu 2024, Tbl N]` | Udrescu et al., *Nat. Phys.* **20**, 202 (2024), DOI 10.1038/s41567-023-02296-w (Zotero `storage/MSQLIKLJ`) |
   | `[PGopher: Linear §…]` | PGopher online linear-molecule Hamiltonian documentation (read directly as primary source) |
   | `[Thesis App A.3.1]` | A. Jadbabaie, Springer Thesis 2025 (`~/Downloads/AJ_Springer_Thesis.pdf`) |
   | `[molecule_parameters.py:LL]` | Molecule-Structure codebase, file:line |
   | `[n2r2-spec §…]` | `docs/superpowers/specs/2026-05-15-n2-r2-formalism-converter-design.md` |
   | `[per Arian, 2026-05-19]` | The user, as the ultimate domain authority |
   | `⚠[OPEN]` | Unresolved — no source available, pending the user's ruling |

2. **No assumptions.** Anything not citable is *not written as a claim*. It becomes either a question to the user or a visible `⚠[OPEN]` callout. Guessed signs, assumed values, and "the standard convention is…" without a citation are forbidden (matches `brown-carrington` skill red-flags and CLAUDE.md epistemic rules).

3. **The user is a citable authority, consulted at every assumption point.** Sign ambiguities, value divergences, and interpretations not stated by any source are surfaced to the user; the answer is recorded as `[per Arian, 2026-05-19]`. The implementation plan must contain explicit user-consultation checkpoints; it is not autonomous-to-completion.

4. **B&C formulas are transcribed verbatim, never paraphrased.** Equation number + PDF page cited every time. Sourced via `mcp__pdf-mcp__pdf_search` / `pdf_read_pages`, not memory (per `tool-selection.md`: primary-source fidelity overrides any context-mode routing hint on `mcp__pdf-mcp__*`).

## 3. Scope

**In scope:** ²²⁶Ra¹⁹F only.

- **X²Σ⁺** (Hund's case (b)): ground electronic state.
- **A²Π₁/₂** (Hund's case (a)): the laser-cooling excited state.
- **¹⁹F hyperfine only** (I(¹⁹F)=½). ²²⁶Ra has I=0 → no Ra hyperfine, no nuclear electric quadrupole.

**Constant inventory** (from the `molecules['RaF']['boson']` X0/A0 dicts, `molecule_parameters.py:43-79` — values to re-verify at build, not asserted here):

| Group | Constants | State(s) |
|---|---|---|
| Rotational | `B` (Be), `D` | X²Σ⁺ & A²Π |
| Fine structure | `γ` (Gamma_SR) | X²Σ⁺ |
| | `A_SO` (ASO) | A²Π |
| Λ-doubling | `p`, `q`, `p+2q`, `p₂q_D` (p2q_D) | A²Π |
| Hyperfine (¹⁹F, I=½) | `b_F` (bF), `c` | X²Σ⁺ |
| | `a`, `d`, `h₁/₂` (= a−(b_F+2c/3) = A‖/2) | A²Π |
| Field-interaction* | `g_S`, `g_l`, `g_l′` (g_lp) ; `μ_E` (muE) | A²Π ; X & A |
| Energy | `T₀` / `Origin` | A relative to X |

\*g-factors and μ_E are Stark/Zeeman parameters, **not** field-free effective-Hamiltonian terms. They are grouped separately and the distinction is stated explicitly in the document, not lumped silently.

**Out of scope (non-goals):** ²²⁵RaF / fermion branch; Ra hyperfine; nuclear electric-quadrupole machinery (eQq); other isotopologues (²²³, etc.); polyatomic structure (Renner-Teller, K-resonance, ℓ); spectral simulation / fitting workflow; any change to the codebase or its constants.

## 4. Document architecture

### 4.1 Front matter
- Scope statement (²²⁶Ra¹⁹F; X²Σ⁺ case (b) / A²Π₁/₂ case (a)).
- **Provenance & Citation Policy box** — states invariant #1–#3 so any reader knows every claim is sourced and what `[per Arian, …]` / `⚠[OPEN]` mean.
- Master constant table (symbol · value · units · isotope · primary source), all cited.
- **Conventions Preamble**: the N²↔R² formalism relation (`[Thesis App A.3.1]`, `[n2r2-spec §…]`, codebase `formalism:'N2'` → R² converter); Condon–Shortley + Brown/Carrington/Hirota parity E\* = σ_xz R_y(π) `[B&C …]`; the PGopher↔B&C sign-flip ledger (`d`, `p+2q` minus-flips noted in `molecule_parameters.py`).

### 4.2 Per-constant entry (identical six-part schema for every constant)

1. **Header** — symbol · name · electronic state · Hund's case · codebase key & stored value (units, ²²⁶Ra¹⁹F) `[molecule_parameters.py:LL]` · primary-source provenance.
2. **Hamiltonian form** — the operator term in the appropriate case basis (X²Σ⁺ → case (b); A²Π₁/₂ → case (a)), transcribed verbatim `[B&C eq. X.XX, p.NNN]`.
3. **Brown & Carrington pages** — operator definition (Ch 4), case reduction / matrix element (Ch 5), fitted-constant meaning (Ch 7), Σ/Π worked example (Ch 9–11). Multiple eq+page citations.
4. **Effective-Hamiltonian contributions** — microscopic interactions folded in via Van Vleck / contact transformation (e.g. γ = γ⁽¹⁾+γ⁽²⁾, γ⁽²⁾ from 2nd-order spin-orbit×rotation, sign-significant in heavy RaF; q ~ pB/A curl relation; microscopic vs effective A_SO). Each statement cited; quantitative cross-checks against `[Udrescu 2024, Tbl N]`.
5. **Conventions — three-column table** `B&C │ N²-paper (Udrescu/codebase-stored) │ PGopher` with explicit transformation/sign map (PGopher minus-flips on `d`, `p+2q`; N²→R² converter row; phase, prefactor, interpretation notes). Every cell cited.
6. **→ Appendix N.x** — cross-reference to the numbered appendix entry holding the background/derivation for this constant.

### 4.3 Appendices (the "layered" depth target)
- **A. Formalism** — full N²↔R² transformation derivation and the engine's converter rows `[Thesis App A.3.1]`, `[n2r2-spec]`.
- **B. Contact-transformation sketches** — γ⁽²⁾ and Λ-doubling p,q origins `[B&C Ch 7]`.
- **C. PGopher LinearMolecule Hamiltonian reference** — relevant operator definitions and sign conventions, from the online doc `[PGopher: Linear §…]`.
- **D. Known divergences & their resolution** — every point where sources disagreed and how it was resolved (`[per Arian, …]` or `⚠[OPEN]`).
- **E. Bibliography** — full source list with Zotero keys / paths.

### 4.4 Known divergences to resolve with the user (seed list; expanded at build)
Surfaced, never silently reconciled:
- **A_SO**: `molecule_parameters.py` `ASO = 1350 cm⁻¹` ("Fixed") vs `225_RaF.pgo` `A = 2067.6 cm⁻¹`.
- **Band origin**: codebase `Origin = 13284.427` vs `225_RaF.pgo` `Origin = 14318.33`.
- **B (A²Π)**: codebase `0.191015 cm⁻¹` vs `225_RaF.pgo` `0.19110 cm⁻¹`.
- Any sign that cannot be pinned in B&C verbatim.

## 5. Build pipeline

- **Single internal source** (pandoc-Markdown + LaTeX math) → two standalone deliverables. One transcription of every equation → no LaTeX/HTML sign-drift (the correctness-critical reason for single-source).
- **LaTeX target**: standalone, CTAN-standard packages only (`amsmath`, `amssymb`, `booktabs`, `mdframed` or `tcolorbox`, `hyperref`, `longtable`); compiles on Overleaf with `pdflatex` and no local filters/shell-escape. Three-column convention tables via `booktabs`/`longtable`; appendix layering via `\appendix`.
- **HTML target**: single `.html`, MathJax-3 from CDN configured for **SVG output** (so browser-rendered equations copy into OneNote as images), self-contained inline CSS only (no JS beyond MathJax, no `<details>` — neither OneNote nor the workflow supports it), bordered/tinted blocks for set-apart content, plain `<table>` for the three-column comparisons.
- Build tooling (pandoc, Lua filter, CSS) lives with the scaffolding; only the two rendered standalone files are delivered to the vault. Exact tooling chosen in the implementation plan.

## 6. Sources & verification protocol

| Source | Access method | Use |
|---|---|---|
| Brown & Carrington | pdf-mcp on Zotero `storage/CKZKCGXY`; `pdf_search`/`pdf_read_pages`, verbatim, read ±2–4 pp around hits, triangulate Ch 4/5/7/9–11 | Hamiltonian forms, matrix elements, sign conventions |
| Udrescu et al. Nat. Phys. 2024 | Zotero `storage/MSQLIKLJ`, read directly | Fitted constant values, N²-paper convention |
| PGopher online doc | Read directly (primary source for sign conventions) | PGopher column, sign-flip ledger |
| Springer thesis App A.3.1 | `~/Downloads/AJ_Springer_Thesis.pdf` (pdf-mcp) | N²↔R² formalism |
| n2-r2 converter spec | `docs/superpowers/specs/2026-05-15-n2-r2-formalism-converter-design.md` | What the engine actually applies |
| Codebase | `molecule_parameters.py`, `matrix_elements.py`, `hamiltonian_builders.py` | Stored values, the form the engine builds |
| **Arian** | Direct consultation at every assumption point | Authority on unresolved signs/divergences/interpretation |

Per `physics-verification.md`: state isotope (²²⁶Ra¹⁹F) and units (MHz; cm⁻¹ where the codebase multiplies by `c`) per constant; flag ab-initio vs experimental; cross-reference values across sources and report divergences.

## 7. Acceptance criteria

- [ ] Every constant in §3's inventory has a complete six-part entry.
- [ ] Zero uncited sentences (citation-tag audit passes).
- [ ] Every B&C citation is eq#+PDF page, verbatim, verified against the PDF this build.
- [ ] Three-column conventions table present and cited for every constant.
- [ ] All §4.4 divergences resolved as `[per Arian, …]` or rendered as visible `⚠[OPEN]`.
- [ ] `RaF-effective-hamiltonian-constants.tex` compiles clean on Overleaf (`pdflatex`, standard packages).
- [ ] `RaF-effective-hamiltonian-constants.html` renders in a browser with MathJax SVG; equations survive browser→OneNote copy as images.
- [ ] Both files in `~/obsidian-vault/04-resources/`; vault committed.
- [ ] No change to any codebase file or constant.

## 8. Non-goals / explicitly deferred

²²⁵RaF, Ra hyperfine, eQq quadrupole, other isotopologues, polyatomic effective Hamiltonians, spectral simulation, fitting, codebase edits. Detailed task breakdown deferred to the writing-plans phase (this spec is scaffolding; per code-style rule, durable physics conclusions that emerge land in the vault, not here).
