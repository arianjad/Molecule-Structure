# RaF Effective-Hamiltonian Constants Reference — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a citation-bound LaTeX + HTML reference for the ²²⁶Ra¹⁹F effective-Hamiltonian constants used in Molecule-Structure, per the approved spec.

**Architecture:** One pandoc-Markdown source (LaTeX math) → standalone Overleaf-ready `.tex` + standalone MathJax-SVG `.html`. A citation-lint script enforces "zero uncited claims" as a build gate (the test-first harness). Physics content is produced by *verified primary-source lookup* (B&C verbatim via pdf-mcp, Udrescu 2024, PGopher doc, codebase) with explicit Arian-consultation gates for anything sources do not settle.

**Tech Stack:** pandoc 3.9, MathJax-3 (SVG, CDN), `tectonic` (LaTeX compile gate, conda `claude-code`), Python 3 (citation-lint), pdf-mcp (B&C/Udrescu/thesis), Playwright via `document-skills:webapp-testing` (HTML render check).

**Spec:** `docs/superpowers/specs/2026-05-19-raf-hamiltonian-constants-reference-design.md` (invariants §2 are non-negotiable).

**Adaptation (honest):** This is a research deliverable, not code. TDD spirit is preserved: Task 0 builds the citation-lint + build + render harness *first* (it must pass on a valid empty skeleton before any content). Every content task must pass citation-lint + build + a B&C-page spot-recheck before commit. The plan does **not** pre-write physics claims or B&C equations/pages — fabricating uncited content would violate spec invariants #2/#4. Content tasks specify the exact verified-lookup procedure and acceptance gate; the physics text is produced by executing those cited lookups during execution.

**Paths:**
- Build dir (source + tooling + working files): `~/obsidian-vault/04-resources/references/raf-effective-hamiltonian/`
- Deliverables (same dir): `RaF-effective-hamiltonian-constants.tex`, `RaF-effective-hamiltonian-constants.html`
- `git` on macOS: `/Library/Developer/CommandLineTools/usr/bin/git` (REPL rule). Molecule-Structure repo for plan/spec; obsidian-vault repo for deliverables. **Commit freely (reversible); do NOT push either repo without explicit authorization.**

---

## File Structure

| File | Responsibility |
|---|---|
| `…/raf-effective-hamiltonian/source.md` | Single source of truth: pandoc-Markdown + LaTeX math; front matter, per-constant entries, appendices A–F |
| `…/raf-effective-hamiltonian/overleaf-template.tex` | pandoc LaTeX template — full Overleaf-safe preamble (CTAN-standard packages only) |
| `…/raf-effective-hamiltonian/style.css` | Standalone HTML CSS (set-apart blocks, 3-col tables); embedded into the .html at build |
| `…/raf-effective-hamiltonian/cite-lint.py` | Citation auditor: fails if any content paragraph lacks a §2-taxonomy citation token |
| `…/raf-effective-hamiltonian/build.sh` | pandoc source.md → .tex (+template) and → .html (+MathJax SVG, embedded CSS); runs cite-lint first |
| `…/raf-effective-hamiltonian/sources-index.md` | Working file (not delivered): per-constant located B&C eq/pages, Udrescu table, PGopher §, codebase line |
| `…/raf-effective-hamiltonian/pgopher-linear.html` | Saved PGopher online linear-molecule Hamiltonian doc (primary source, archived) |
| `RaF-effective-hamiltonian-constants.{tex,html}` | The two delivered standalone renderings |

Dependency order: **0 → 1 → 2 → {3,4,5,6,7,8} → 9 → 10**. Tasks 3–8 are mutually independent (all depend on 1 for sources and 2 for the convention spine).

---

### Task 0: Verify harness — citation-lint + build + render gate on an empty-but-valid skeleton

**Goal:** Build dir, source skeleton (front matter + section/appendix stubs incl. reserved Appendix F + Provenance box), Overleaf template, HTML CSS+MathJax, `cite-lint.py`, `build.sh`; prove the gate passes on a valid skeleton before any physics content (test-first).

**Files:**
- Create: `~/obsidian-vault/04-resources/references/raf-effective-hamiltonian/{source.md,overleaf-template.tex,style.css,cite-lint.py,build.sh}`

**Acceptance Criteria:**
- [ ] `cite-lint.py source.md` exits 0 on the skeleton (skeleton has no uncited *claims*; `⚠[OPEN]`/`[per Arian…]` tokens and headings/code/tables are exempt by rule).
- [ ] `build.sh` produces `RaF-effective-hamiltonian-constants.tex` and `.html`, both standalone.
- [ ] `.tex` compiles: `tectonic` exits 0 (or, if `tectonic` uninstallable offline, the structural package-allowlist check passes and the audit records "Overleaf-pending").
- [ ] `.html` renders in headless Chromium with MathJax typesetting a probe equation (no raw `$…$` visible) — screenshot inspected.

**Verify:** `cd <build dir> && python3 cite-lint.py source.md && bash build.sh && conda run -n claude-code tectonic RaF-effective-hamiltonian-constants.tex` → all exit 0; then render check (Step 6).

**Steps:**

- [ ] **Step 1: Create build dir + skeleton `source.md`.** Skeleton contains, with NO physics claims: title; **Provenance & Citation Policy** box (verbatim restatement of spec §2 invariants + the citation-token table); scope line; `## Conventions Preamble` (stub: single line `Convention spine — populated in Task 2. ⚠[OPEN — Task 2]`); `## Constants` with one `### <name>` stub per spec-§3 inventory entry, each containing the literal line `⚠[OPEN — content pending]`; `## Appendix A … F` headings, with `### F. ²²⁵RaF / fermion extension — reserved` containing `⚠[OPEN — future pass]` and the enumerated reserved constant list from spec §3. Use Write tool (not heredoc).

- [ ] **Step 2: Write `cite-lint.py`.** Exact behavior: read the markdown; split into blocks separated by blank lines; skip blocks that are (a) ATX headings, (b) fenced code, (c) Markdown tables (line starts with `|`), (d) the Provenance box (delimited by HTML comments `<!--PROVENANCE-START-->`/`<!--PROVENANCE-END-->`), (e) blocks whose stripped content is only a display-math fence `$$…$$`. For every remaining block, require ≥1 regex match of:
  `r"\[B&C eq\.|\[B&C §|\[Udrescu 2024|\[PGopher:|\[Thesis App|\[molecule_parameters\.py:|\[matrix_elements\.py:|\[hamiltonian_builders\.py:|\[n2r2-doc|\[per Arian|⚠\[OPEN"`.
  Print `FAIL <line-no>: <first 80 chars>` for each uncited block; exit 1 if any, else print `cite-lint: PASS (<n> content blocks, all cited)` exit 0. Complete script written in this step (no placeholder).

- [ ] **Step 3: Write `overleaf-template.tex`.** pandoc LaTeX template; preamble uses ONLY: `amsmath, amssymb, booktabs, longtable, array, mdframed, xcolor, hyperref, geometry, parskip` (all in Overleaf's default TeX Live). No `\write18`, no shell-escape, no local `\input`. Include `$body$`, `$title$` pandoc vars. Provenance box → `mdframed`.

- [ ] **Step 4: Write `style.css`** — readable serif body, `.provenance`/`.background` set-apart bordered tinted blocks, `table{border-collapse}` + `th,td` borders for the 3-col tables, max-width for OneNote-friendly paste.

- [ ] **Step 5: Write `build.sh`.** Exactly:
  ```bash
  set -euo pipefail
  cd "$(dirname "$0")"
  python3 cite-lint.py source.md
  pandoc source.md --standalone --template=overleaf-template.tex \
    -o RaF-effective-hamiltonian-constants.tex
  pandoc source.md --standalone --embed-resources --css=style.css \
    --mathjax='https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js' \
    -o RaF-effective-hamiltonian-constants.html
  echo "build: OK"
  ```

- [ ] **Step 6: Install compile gate + run full verify.**
  Run: `conda run -n claude-code bash -lc 'command -v tectonic || conda install -n claude-code -c conda-forge tectonic -y'`
  Expected: tectonic resolves. If offline/unavailable: add `pkg-allowlist-check.py` (assert `overleaf-template.tex` `\usepackage` set ⊆ allowlist; grep absent `write18`/local `\input`) and record "Overleaf-pending" — do NOT claim compiled.
  Run: `bash build.sh && conda run -n claude-code tectonic RaF-effective-hamiltonian-constants.tex`
  Expected: `cite-lint: PASS`, `build: OK`, tectonic exit 0, PDF produced.
  Render check: use `document-skills:webapp-testing` (Playwright) to open the `.html` `file://` URL, wait for MathJax, screenshot; `await Read()` the screenshot in REPL; confirm the probe equation is typeset, not raw `$…$`.

- [ ] **Step 7: Commit (Molecule-Structure repo for plan/spec is already committed; this commits the build harness to the vault).**
  ```bash
  cd ~/obsidian-vault && /Library/Developer/CommandLineTools/usr/bin/git add 04-resources/references/raf-effective-hamiltonian
  /Library/Developer/CommandLineTools/usr/bin/git -c user.name="Arian Jadbabaie" commit -m "raf-ref: verify harness (cite-lint + build + render) on valid skeleton"
  ```
  (No push.)

---

### Task 1: Source resolution & B&C chapter mapping (read-only — zero claims written)

**Goal:** Resolve every primary source and pre-locate (not transcribe) per-constant B&C equations/pages, Udrescu tables, PGopher sections, codebase lines, into `sources-index.md`. Front-loads lookup so content tasks are verify+transcribe, not search.

**Files:**
- Create: `…/raf-effective-hamiltonian/sources-index.md`, `…/raf-effective-hamiltonian/pgopher-linear.html`

**Acceptance Criteria:**
- [ ] B&C PDF opens via `mcp__pdf-mcp__pdf_info` (Zotero `storage/CKZKCGXY`); Udrescu 2024 (`storage/MSQLIKLJ`) and thesis (`~/Downloads/AJ_Springer_Thesis.pdf`) open.
- [ ] PGopher online linear-molecule Hamiltonian doc fetched and saved to `pgopher-linear.html`.
- [ ] `docs/n2_r2_formalism_converter.md` + `Source Code/formalism.py` docstring read (the **authoritative** N²↔R²; the 2026-05-15 spec is SUPERSEDED — use only for the still-accurate §4 Table 7.2, labeled historical).
- [ ] Every spec-§3 inventory constant has, in `sources-index.md`: ≥1 candidate B&C eq+page (from `pdf_search`, not yet transcribed), Udrescu table ref, PGopher § ref, and `molecule_parameters.py` line.

**Verify:** `python3 -c "import re,sys; t=open('sources-index.md').read(); [sys.exit('MISSING '+k) for k in ['B','D','gamma','A_SO','p+2q','p2q_D','b_F','c','a','d','h1/2','g_S','g_l','muE','Origin'] if k not in t]; print('index OK')"`

**Steps:**
- [ ] **Step 1: Open/confirm PDFs** — `mcp__pdf-mcp__pdf_info` on B&C, Udrescu, thesis. Record page counts.
- [ ] **Step 2: Fetch PGopher doc.** Per `tool-selection.md`, the PGopher convention text is citation-grade → read directly (not context-mode digest): `document-skills` Playwright or WebFetch the PGopher manual "Linear molecule"/"Hamiltonian" pages; save raw HTML to `pgopher-linear.html`; note the canonical URL for citation `[PGopher: Linear §…]`.
- [ ] **Step 3: Read authoritative N²↔R².** `cat docs/n2_r2_formalism_converter.md`; `sed -n` the `formalism.py` module docstring + `CENTRIFUGAL_PARTNERS`. Note exact converter rows for `Be, ASO, p+2q, Origin, D`.
- [ ] **Step 4: Per-constant B&C location pass.** For each spec-§3 constant, `mcp__pdf-mcp__pdf_search` with the precise operator term (use Unicode Λ; e.g. "spin–rotation", "Λ-doubling", "Fermi contact", "spin–orbit coupling", "centrifugal distortion"). Record top ~5 hits (eq# + page) and tag which is Ch4 operator / Ch5 case-reduction / Ch7 effective-constant / Ch9–11 example. **Do not transcribe** — locate only.
- [ ] **Step 5: Codebase lines.** `grep -n` `molecule_parameters.py` (RaF boson X0/A0), `matrix_elements.py`, `hamiltonian_builders.py` for each constant's engine term; record file:line.
- [ ] **Step 6: Write `sources-index.md`** (working file, no deliverable claims) and commit to vault (`raf-ref: source index + PGopher doc archived`).

---

### Task 2: Convention spine — Conventions Preamble + Appendix A (N²↔R²) + Appendix C (PGopher) + parity/sign ledger

**Goal:** Author the convention backbone every per-constant 3-col table references. **First Arian-consultation gate.**

**Files:** Modify `source.md` (Conventions Preamble; Appendix A; Appendix C; PGopher↔B&C sign-flip ledger).

**Acceptance Criteria:**
- [ ] N²↔R² relation stated with the engine's actual converter rows, cited `[n2r2-doc …]` + `[formalism.py:LL]` + `[B&C eq. 9.138 …]` (engine R² = N² − Λ²·I), verbatim where a formula.
- [ ] Parity convention (Condon–Shortley; Brown/Carrington/Hirota E\*=σ_xz R_y(π)) stated, B&C verbatim + page.
- [ ] PGopher↔B&C sign-flip ledger: `d`, `p+2q` (and any other) with the PGopher-doc statement quoted `[PGopher: …]` and the codebase comment cited `[molecule_parameters.py:LL]`.
- [ ] `cite-lint.py` PASS; `build.sh` OK; tectonic exit 0.
- [ ] Every point sources did NOT settle raised to Arian; answers written `[per Arian, 2026-05-19]`; unresolved → visible `⚠[OPEN]`.

**Verify:** `bash build.sh && conda run -n claude-code tectonic RaF-effective-hamiltonian-constants.tex` → exit 0; manual: open 2 cited B&C pages, confirm transcription matches.

**Steps:**
- [ ] **Step 1:** Draft Conventions Preamble + Appendix A from Task-1 authoritative N²↔R² sources; B&C eq. 9.138 transcribed verbatim via `pdf_read_pages` (read ±3pp).
- [ ] **Step 2:** Parity convention — `pdf_search` "E\* " / "Condon" / parity in Ch 3/6; transcribe verbatim + page.
- [ ] **Step 3:** PGopher Appendix C + sign ledger from `pgopher-linear.html` (quote the sign-defining sentence) cross-referenced to codebase comments.
- [ ] **Step 4: Arian-consultation gate.** Compile the list of anything unsettled (e.g. exact PGopher `d` sign sentence ambiguity, any N²↔R² edge). Present to Arian; record answers `[per Arian, 2026-05-19]`; leave residue as `⚠[OPEN]`. Do not guess.
- [ ] **Step 5:** `bash build.sh`; tectonic; fix until green. Commit vault (`raf-ref: convention spine (preamble + App A/C + sign ledger)`).

---

### Task 3: Rotational — B, D (X²Σ⁺ & A²Π)

**Goal:** Full six-part entries for `B` and `D` in both states.

**Files:** Modify `source.md` (`### B`, `### D`).

**Acceptance Criteria:**
- [ ] Six-part schema complete for B and D (Header w/ codebase value+line; Hamiltonian form verbatim B&C; B&C Ch4/5/7/9–11 pages; eff-H contributions cited; 3-col B&C/N²/PGopher table cited; reserved `⚠[OPEN — ²²⁵RaF]` block; → Appendix xref).
- [ ] Case basis correct: X²Σ⁺ → case (b); A²Π₁/₂ → case (a). B&C forms transcribed verbatim, page-cited.
- [ ] `cite-lint.py` PASS; `build.sh` OK; tectonic exit 0.
- [ ] B&C-page spot-recheck: reopen ≥1 cited page per constant, confirm equation text matches transcription.

**Verify:** `bash build.sh && conda run -n claude-code tectonic …tex` exit 0 + spot-recheck logged.

**Steps:**
- [ ] **Step 1:** From `sources-index.md`, `pdf_read_pages` the located B (rotational) pages (Ch 7 effective rotational; Ch 9.138 region) ±2–4pp; transcribe the operator term verbatim per state/case.
- [ ] **Step 2:** Eff-H contributions for B (vibrational averaging Bᵥ; R² vs N² centrifugal-partner coupling per the converter) — each sentence cited.
- [ ] **Step 3:** 3-col table: B&C form │ N²-paper (Udrescu Tbl, codebase stored) │ PGopher `B` — with the N²→R² converter row cited.
- [ ] **Step 4:** Repeat Steps 1–3 for D (centrifugal distortion). Add reserved ²²⁵RaF block + Appendix B/D xref.
- [ ] **Step 5:** Any divergence (e.g. A-state B 0.191015 vs `225_RaF.pgo` 0.19110) → record in Appendix D and flag for the Task-9 Arian gate; do not reconcile here. `build.sh`; tectonic; commit (`raf-ref: B, D entries`).

---

### Task 4: Fine structure — γ (X²Σ⁺ spin-rotation), A_SO (A²Π spin-orbit) — divergence-bearing

**Goal:** Six-part entries for γ and A_SO. A_SO carries the 1350 vs 2067.6 cm⁻¹ divergence.

**Files:** Modify `source.md` (`### gamma`, `### A_SO`).

**Acceptance Criteria:**
- [ ] γ: case-(b) `γ N·S` form verbatim B&C + page; eff-H contributions incl. γ = γ⁽¹⁾+γ⁽²⁾ with γ⁽²⁾ (2nd-order spin-orbit×rotation) heavy-molecule dominance/sign — every clause cited (B&C Ch7 + Udrescu + `about-arian-physics` only as orientation, not as the citation).
- [ ] A_SO: case-(a) `A Lz·Sz` form verbatim B&C + page; microscopic-vs-effective distinction cited.
- [ ] **A_SO divergence (1350 vs 2067.6) NOT silently reconciled** — both values with provenance, flagged to Task-9 Arian gate, rendered `⚠[OPEN]` until ruled.
- [ ] 3-col tables cited; reserved ²²⁵RaF block; `cite-lint` PASS; build+tectonic exit 0; spot-recheck logged.

**Verify:** `bash build.sh && conda run -n claude-code tectonic …tex` exit 0; grep shows A_SO has provenance for both values + `⚠[OPEN]`.

**Steps:**
- [ ] **Step 1:** `pdf_read_pages` γ operator (Ch4) + case-(b) reduction (Ch5) + effective γ (Ch7); transcribe verbatim.
- [ ] **Step 2:** γ eff-H contributions, cited (B&C Ch7 second-order treatment; Udrescu fitted value cross-check).
- [ ] **Step 3:** A_SO operator (Ch4) + case-(a) (Ch5) verbatim; record codebase 1350 `[molecule_parameters.py:60]` vs `225_RaF.pgo` 2067.6 — both, with source, `⚠[OPEN]`, Appendix D entry.
- [ ] **Step 4:** 3-col tables; reserved blocks; xrefs. `build.sh`; tectonic; commit (`raf-ref: γ, A_SO entries (A_SO divergence flagged)`).

---

### Task 5: Λ-doubling — p, q, p+2q, p₂q_D (A²Π)

**Goal:** Six-part entries for the Λ-doubling parameters; pin the `p+2q` combination and PGopher sign.

**Files:** Modify `source.md` (`### p`, `### q`, `### p+2q`, `### p2q_D`).

**Acceptance Criteria:**
- [ ] Λ-doubling operator verbatim B&C (Ch7/Ch9) + page; the `o,p,q` (Brown) ↔ `p+2q` effective combination for ²Π₁/₂ stated and cited; q=0 codebase choice cited `[molecule_parameters.py:LL]`.
- [ ] Curl relation q ~ pB/A noted as eff-H contribution, B&C-cited (not from memory).
- [ ] PGopher `p+2q` sign vs B&C: ledger reference + the codebase minus-flip comment cited.
- [ ] 3-col tables cited; reserved ²²⁵RaF block; `cite-lint` PASS; build+tectonic exit 0; spot-recheck logged.

**Verify:** `bash build.sh && conda run -n claude-code tectonic …tex` exit 0.

**Steps:**
- [ ] **Step 1:** `pdf_search`/`pdf_read_pages` "Λ-doubling" Ch7 + ²Π worked example Ch9–11; transcribe o/p/q verbatim.
- [ ] **Step 2:** Effective combination + curl relation, cited; p₂q_D centrifugal correction form cited.
- [ ] **Step 3:** 3-col tables incl. PGopher sign-flip (ref Appendix C ledger). Any unsettled sign → Task-9 gate.
- [ ] **Step 4:** `build.sh`; tectonic; commit (`raf-ref: Λ-doubling entries`).

---

### Task 6: Hyperfine X²Σ⁺ — b_F, c (¹⁹F, I=½)

**Goal:** Six-part entries for Fermi-contact `b_F` and dipolar `c` in X²Σ⁺.

**Files:** Modify `source.md` (`### b_F`, `### c`).

**Acceptance Criteria:**
- [ ] Frosch–Foley relation b_F = b + c/3 stated verbatim B&C + page (the BaF comment cites "B&C Eq.8.511"; verify the RaF-relevant eq independently, do not copy the BaF line uncritically).
- [ ] Case-(b) hyperfine matrix-element form verbatim B&C + page; ¹⁹F I=½ stated.
- [ ] 3-col tables cited; reserved ²²⁵RaF block names BOTH heavy-Ra and ¹⁹F slots; `cite-lint` PASS; build+tectonic exit 0; spot-recheck logged.

**Verify:** `bash build.sh && conda run -n claude-code tectonic …tex` exit 0.

**Steps:**
- [ ] **Step 1:** `pdf_read_pages` Ch4 hyperfine operators (Fermi contact, dipolar) + Ch8/9 b_F=b+c/3; transcribe verbatim.
- [ ] **Step 2:** Eff-H contributions (spin density at ¹⁹F; b vs b_F convention), cited.
- [ ] **Step 3:** 3-col tables; reserved 2-nucleus block. `build.sh`; tectonic; commit (`raf-ref: X²Σ⁺ hyperfine b_F, c`).

---

### Task 7: Hyperfine A²Π — a, d, h₁/₂ (¹⁹F) — PGopher-sign focus

**Goal:** Six-part entries for orbital `a`, the Π-state `d`, and the case-(a) `h₁/₂ = a − (b_F + 2c/3) = A‖/2` combination.

**Files:** Modify `source.md` (`### a`, `### d`, `### h1/2`).

**Acceptance Criteria:**
- [ ] a (orbital hyperfine) and d (²Π parity/Λ-doubling hyperfine) operator forms verbatim B&C + page; the case-(a) ²Π₁/₂ combination h₁/₂ = a−(b_F+2c/3) = A‖/2 stated and cited (codebase formula `[molecule_parameters.py:55]`).
- [ ] **`d` PGopher minus-flip**: PGopher-doc sign sentence quoted `[PGopher:…]`, codebase comment cited `[molecule_parameters.py:LL]`, B&C sign verbatim — the three reconciled in the 3-col table (or `⚠[OPEN]` + Task-9 gate if not settleable from sources).
- [ ] 3-col tables cited; reserved ²²⁵RaF block; `cite-lint` PASS; build+tectonic exit 0; spot-recheck logged.

**Verify:** `bash build.sh && conda run -n claude-code tectonic …tex` exit 0; grep shows `d` entry has all three convention columns + sign reconciliation or `⚠[OPEN]`.

**Steps:**
- [ ] **Step 1:** `pdf_read_pages` Ch4 a,d hyperfine + Ch9–11 ²Π hyperfine example; transcribe verbatim; derive/cite h₁/₂ identity from B&C (not memory).
- [ ] **Step 2:** PGopher `d` sign: quote `pgopher-linear.html`; reconcile vs B&C vs codebase comment; unsettled → Task-9 gate.
- [ ] **Step 3:** 3-col tables; reserved blocks. `build.sh`; tectonic; commit (`raf-ref: A²Π hyperfine a, d, h₁/₂`).

---

### Task 8: Field-interaction (g_S, g_l, g_l′, μ_E) + Energy origin T₀ — divergence-bearing

**Goal:** Six-part entries for the Zeeman g-factors and electric dipole μ_E (explicitly labelled NOT field-free effective-H) and the band origin T₀ (carries the 13284 vs 14318 divergence).

**Files:** Modify `source.md` (`### g_S`, `### g_l`, `### g_l'`, `### muE`, `### T0`).

**Acceptance Criteria:**
- [ ] g_S/g_l/g_l′ Zeeman operator forms (Ch6) verbatim B&C + page; μ_E Stark operator (Ch6) verbatim; an explicit boxed note that these are field-interaction, not field-free effective-H, cited to B&C Ch6 framing.
- [ ] T₀/Origin: definition + the engine's R² G-row offset cited `[n2r2-doc]`/`[formalism.py:LL]`; **divergence 13284.427 (codebase) vs 14318.33 (`225_RaF.pgo`) NOT reconciled** — both with provenance, `⚠[OPEN]`, Task-9 gate, Appendix D.
- [ ] 3-col tables cited; reserved ²²⁵RaF block; `cite-lint` PASS; build+tectonic exit 0; spot-recheck logged.

**Verify:** `bash build.sh && conda run -n claude-code tectonic …tex` exit 0; grep shows T₀ has both values + `⚠[OPEN]`.

**Steps:**
- [ ] **Step 1:** `pdf_read_pages` Ch6 Zeeman (g_S, g_l, parity-dependent g_l′) + Stark (μ_E); transcribe verbatim.
- [ ] **Step 2:** Field-interaction caveat box, B&C Ch6-cited. g_lp/g_l codebase mapping cited `[molecule_parameters.py:LL]` — state which codebase key maps to which B&C g; if ambiguous, `⚠[OPEN]` + Task-9 gate (do not assume the mapping).
- [ ] **Step 3:** T₀ definition + R² G-row; record both origin values w/ provenance, `⚠[OPEN]`, Appendix D.
- [ ] **Step 4:** 3-col tables; reserved blocks. `build.sh`; tectonic; commit (`raf-ref: g-factors, μ_E, T₀ (origin divergence flagged)`).

---

### Task 9: Appendix D consolidation + single consolidated Arian-consultation gate

**Goal:** Collect EVERY `⚠[OPEN]` and source disagreement from Tasks 2–8 into Appendix D; present the full list to Arian in ONE consultation; record rulings.

**Files:** Modify `source.md` (Appendix D; resolved entries throughout).

**Acceptance Criteria:**
- [ ] Appendix D lists every divergence/open item with provenance of each candidate value/sign.
- [ ] All items presented to Arian in one consolidated message; each answer written into the relevant entry AND Appendix D as `[per Arian, 2026-05-19]`.
- [ ] Items Arian leaves open remain visible `⚠[OPEN]` (never silently resolved or guessed).
- [ ] `cite-lint` PASS; build+tectonic exit 0.

**Verify:** `grep -c "⚠\[OPEN" source.md` equals the count Arian explicitly left open (logged in Appendix D); `bash build.sh && conda run -n claude-code tectonic …tex` exit 0.

**Steps:**
- [ ] **Step 1:** Assemble Appendix D from all flagged items (A_SO 1350/2067.6, Origin 13284/14318, A-B 0.191015/0.19110, any g-mapping/sign opens).
- [ ] **Step 2: Arian-consultation gate (single).** Present the consolidated list with each option's provenance; ask for rulings. Do not proceed past unresolved physics by assumption.
- [ ] **Step 3:** Write rulings as `[per Arian, 2026-05-19]` into entries + Appendix D; keep residue `⚠[OPEN]`. `build.sh`; tectonic; commit (`raf-ref: Appendix D + Arian rulings`).

---

### Task 10: Front matter finalize, master table, bibliography, full acceptance audit, deliver to vault

**Goal:** Complete front matter (master constant table, all cited) + Appendix E bibliography; run the full spec §7 acceptance audit; finalize the two deliverables in the vault.

**Files:** Modify `source.md`; produce final `RaF-effective-hamiltonian-constants.{tex,html}`.

**Acceptance Criteria (= spec §7):**
- [ ] Every spec-§3 constant has a complete six-part entry.
- [ ] **Zero uncited sentences** — `cite-lint.py` PASS *and* a manual sentence-level audit of a 15% random sample finds no uncited claim.
- [ ] Every B&C citation is eq#+PDF page, re-verified this build on a ≥20% random sample (reopen page, match text).
- [ ] 3-col conventions table present+cited for every constant.
- [ ] All Task-9 divergences resolved `[per Arian,…]` or visible `⚠[OPEN]`.
- [ ] `.tex` tectonic exit 0 (or Overleaf-pending explicitly recorded); `.html` MathJax-SVG renders in headless Chromium (screenshot inspected) and a select-all copy yields typeset equations (simulated; OneNote paste is the user's manual final step).
- [ ] Reserved ²²⁵RaF structure present (per-entry block + Appendix F skeleton); adding it later needs no reflow.
- [ ] Both files in `~/obsidian-vault/04-resources/references/raf-effective-hamiltonian/`; vault committed.
- [ ] No codebase file/constant changed (`git -C ~/Code/Molecule-Structure status` clean except docs/).

**Verify:** `cd <build dir> && python3 cite-lint.py source.md && bash build.sh && conda run -n claude-code tectonic RaF-effective-hamiltonian-constants.tex && echo AUDIT-OK`; render screenshot inspected; sample re-verification logged in `sources-index.md`.

**Steps:**
- [ ] **Step 1:** Master constant table (symbol·value·units·isotope·source) — every row cited to `[molecule_parameters.py:LL]` + paper.
- [ ] **Step 2:** Appendix E bibliography (B&C edition+Zotero key, Udrescu DOI, thesis, PGopher URL, n2r2-doc, codebase commit).
- [ ] **Step 3:** Full acceptance audit per criteria above; fix any miss inline; re-run gate.
- [ ] **Step 4:** Final `build.sh`; render check; commit vault (`raf-ref: final — front matter, bibliography, acceptance audit`). **Ask Arian before any `git push`** of the vault (push is the careful outward action).

---

## Self-Review

**1. Spec coverage:** §1 goal → Tasks 0/10. §2 invariants → cite-lint (T0) enforces #1; T2/T9 gates enforce #2/#3; "B&C verbatim" in every content task enforces #4. §3 inventory → every constant mapped to T3–T8; ²²⁵RaF reserved → T0 skeleton + every content task's reserved block + Appendix F. §4 architecture → T0 (skeleton/appendices), T2 (preamble/spine), T3–T8 (entries), T10 (front matter/bib). §4.4 divergences → T4/T8 flag, T9 consolidates. §5 build → T0. §6 sources → T1. §7 acceptance → T10 audit. §8 non-goals → respected (no codebase edits; ²²⁵RaF reserved only). No gaps.

**2. Placeholder scan:** Scaffolding tasks (0,1) contain complete scripts/commands. Content tasks (2–9) intentionally specify the *verified-lookup procedure* not pre-written physics — this is required by spec invariants #2/#4 (pre-writing B&C equations/pages here would be fabricated uncited content) and is stated explicitly in the header Adaptation note. Not a placeholder failure: the procedure, files, acceptance, and verify command are concrete in every task.

**3. Type/name consistency:** `source.md`, `cite-lint.py`, `build.sh`, `overleaf-template.tex`, `style.css`, `sources-index.md`, `pgopher-linear.html` and the two deliverable filenames are used identically across all tasks. Citation tokens identical to spec §2 and to the cite-lint regex. Task dependency order stated once and consistent.

**Status:** Plan self-review complete; no inline fixes required.
