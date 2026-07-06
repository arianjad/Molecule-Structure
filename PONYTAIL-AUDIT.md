# Ponytail audit — advisory (2026-07-06)

Over-engineering suggestions for THIS repo from a multi-agent audit (per-repo auditor + adversarial verifier + cross-repo duplication pass; load-bearing items spot-checked by the orchestrating session). **Advisory only — nothing here is a decision.** A session working in this repo, with local context, decides per item: execute, adapt, or reject. Scope is complexity only (no correctness/security/performance).

If executing an item: re-grep the symbol first (INCLUDING *.ipynb — notebooks are callers), make the change, run Verification below, commit (do not push). One coherent batch per commit.

Full machine-readable evidence: `~/.claude/handoffs/2026-07-06-ponytail-audit-findings.json`.

## Do NOT touch

- `config_path.py` copies in 9 notebook dirs — refuted: documented deliberate (`Jupyter Notebooks/CLAUDE.md:27` — keep identical, in sync).

## Suggestions (largest cut first, ~422 LOC total)

### 1. `delete` — ~100 LOC — `Source Code/case_bBJ_MEs.py:1`

- **Suggestion:** case_bBJ_MEs.py — 100-line file of bBJ matrix elements (Rot_bBJ, SR_bBJ, IS_bBJ, own kronecker) that duplicates a subset of matrix_elements.py and is imported nowhere
- **Replacement:** nothing (git history keeps it)
- **Auditor evidence:** grep -rl 'case_bBJ_MEs' over 'Source Code' *.py AND 'Jupyter Notebooks' *.ipynb returned zero hits; git log for the file shows only ancient 'Reorganizing files'/'Cleaning up' commits; read L1-30 — functions are name-identical to matrix_elements.py versions. README lists it in the Source Code map, which is the only reason confidence isn't higher.
- **Skeptic verdict:** confirmed — Whole-repo grep (py/ipynb/md) for 'case_bBJ_MEs' hits only README.md:108 (a table row, not an import). wc -l = exactly 100. Defs Rot_bBJ/SR_bBJ/IS_bBJ are signature-identical to matrix_elements.py:108/115/122. git log shows only ancient reorg commits (af7735f, d5fcb9c). No getattr/eval dynamic dispatch exists in Source Code. Deletion also needs the one README table row removed. (auditor confidence: medium)

### 2. `delete` — ~95 LOC — `Source Code/quantum_numbers.py:5`

- **Suggestion:** quantum_numbers.py 'DEVELOPMENT BEGINNING…END' block (L5-115): count_momenta, recursive_add_J, base_add_J, list_add_J, q_numbers_bBJ_new — a closed dead subgraph; q_numbers_bBJ_new has no callers and is the only consumer of the add_J family
- **Replacement:** delete outright
- **Auditor evidence:** grep for q_numbers_bBJ_new/count_momenta/recursive_add_J across Source Code *.py and all *.ipynb: only intra-block calls (quantum_numbers.py:11-56) plus notebooks (YbOH_Branching.ipynb, Generic/Untitled.ipynb) that define their OWN local copies of count_momenta rather than importing it. molecule_library_class dispatches to q_numbers_even_bBJ/odd_bBJ/decoupled_mJ, never _new. git -S: last touched in old 'Cleaning up file directory' commit. Explicit DEVELOPMENT marker = parked WIP, hence medium not high.
- **Skeptic verdict:** downgraded (main-session spot-check) — DEVELOPMENT markers at quantum_numbers.py:4 and :118 (~115-line block). q_numbers_bBJ_new: zero callers repo-wide incl. ipynb (only the def at :47). molecule_library_class.py dispatch table uses q_numbers_even_bBJ/odd_bBJ/bBS/decoupled_mJ via partial(); even the commented legacy entry references qn.q_numbers_bBJ, not _new. count_momenta/add_J family called only intra-block; YbOH_Branching.ipynb and Generic/Untitled.ipynb define their own local copies (verified def lines in the ipynb JSON). grep for getattr/globals()/eval in the dispatch path: empty. (auditor confidence: medium)
- **⚠ CORRECTION: `count_momenta` is LIVE — 5 uses in `Jupyter Notebooks/YbOH/Even Isotope/YbOH_Branching.ipynb`. Keep `count_momenta` (move above the block); the add_J trio only appears in scratch `Generic/Untitled.ipynb`. ~95 LOC, not 111.**

### 3. `delete` — ~81 LOC — `Source Code/Energy_Levels.py:1057`

- **Suggestion:** MoleculeLevels.gen_state_str_old — 81-line superseded version of gen_state_str kept in the class
- **Replacement:** delete outright
- **Auditor evidence:** grep 'gen_state_str_old' over all *.ipynb (zero hits, verified twice incl. per-file count grep) and over Source Code excluding the def line (zero internal callers). git log -S shows it was orphaned by 'Bug fixes for gen_state_str'.
- **Skeptic verdict:** confirmed — Only occurrence repo-wide (py + ipynb) is the def at Energy_Levels.py:1057; next def (select_q) at :1139, so the method is 82 lines. No dynamic dispatch (getattr/eval grep empty). The live gen_state_str sits directly above it. (auditor confidence: high)

### 4. `delete` — ~48 LOC — `Source Code/Energy_Levels.py:250`

- **Suggestion:** ~45 lines of commented-out mp.Pool/perf_counter timing scaffolding repeated verbatim in AngleMap, ZeemanMap, and StarkMap (superseded by the live diagonalize_batch path), plus the three imports now used ONLY inside those comments: 'import sympy as sy', 'from functools import partial', 'from time import perf_counter'
- **Replacement:** delete outright
- **Auditor evidence:** Read L243-370: each of the three sweep methods carries the same ~15-line commented pool/perf_counter block. grep confirms perf_counter( and partial( appear only on '#' lines (L250-359), and '\bsy\b'/'sy.' appears nowhere outside line 3. npLA/sciLA/torch were also checked and ARE live (diagonalize/diagonalize_batch L1363-1395) — deliberately NOT flagged.
- **Skeptic verdict:** confirmed — grep confirms perf_counter(/partial(/mp.Pool appear ONLY on '#'-commented lines (L250-359, three identical blocks in AngleMap/ZeemanMap/StarkMap, superseded by live diagonalize_batch calls at each site). 'sy.' appears nowhere in Energy_Levels.py (import sympy as sy at L3 is dead), and no file anywhere does 'from Energy_Levels import *', so removing the three imports (L3, L9, L10) breaks no re-export. (auditor confidence: high)

### 5. `shrink` — ~40 LOC — `Source Code/CPT_Sims.py:5`

- **Suggestion:** CPT_Sims.py: four near-copy prep_sim variants (prep_sim, prep_sim_1ph, prep_sim_1d_weighted, prep_sim2d, each ~30 lines sharing the decoherence/decay/time defaults) plus parallel mp_cpt_sim* wrappers
- **Replacement:** one prep_sim(...) with the scan axis as a parameter; keep thin aliases so archival notebooks still import
- **Auditor evidence:** smart_outline shows the 4x prep + 4x mp + 3x sim family with identical default signatures; diff of prep_sim vs prep_sim_1ph bodies shows substantial but not total overlap (27 diff lines). Callers are 5 archival YbOH thesis notebooks (grep *.ipynb) — churn risk on archival deliverables is why confidence is low.
- **Skeptic verdict:** downgraded — The 4x prep_sim / 4x mp_cpt_sim / 3-4x cpt_sim family exists as described (defs at L5/43/90/147 etc., 285-line file) and the signatures share the decoherence/decay/time defaults. But behavior-identity of a merged prep_sim is unproven: the variants differ structurally (prep_sim_1ph takes one_photon_array as leading arg, _1d_weighted and 2d change the sweep pipeline, 27 diff lines per the finding's own diff), callers are 6 (not 5) frozen archival YbOH thesis notebooks, and the repo's verification gate (RaF tutorial notebook) never exercises CPT_Sims, so a regression would ship silently. Plausible consolidation, negative risk/benefit on archival-only code. (auditor confidence: low)

### 6. `shrink` — ~25 LOC — `Source Code/Energy_Levels.py:268`

- **Suggestion:** the ~14-line state-ordering + phase-fix loop (state_ordering, reorder, np.sign(np.diag(...)) phase correction) is copy-pasted three times in AngleMap, ZeemanMap, StarkMap
- **Replacement:** one module-level helper, e.g. track_states(evals_t, evecs_t, initial_evecs, round) -> (evals_t, evecs_t), called from all three sweeps
- **Auditor evidence:** Read L243-394: AngleMap (L268-281), ZeemanMap (L318-333), StarkMap (~L370-385) contain the identical loop body; awk over EB_grid (L1161-1243) confirmed no fourth copy (it doesn't inline the phase-fix), so the extraction touches exactly three verbatim sites.
- **Skeptic verdict:** confirmed — Read all three sweeps: the ordering+phase-fix loop body (state_ordering -> reorder -> np.sign(np.diag(evecs_ordered@evecs_old.T)) -> evecs_old update, incl. the 'i==0 and initial_evecs is None: continue' guard) is verbatim-identical in AngleMap (L268-281), ZeemanMap (L322-335, inside 'if order:'), StarkMap (L371-384, inside 'if order:') modulo array names. A helper taking (evals_t, evecs_t, initial_evecs, round) is behavior-identical; each site keeps its own 'if order:' gate. 3x14 lines -> ~16-line helper + 3 calls, net ~23-25. (auditor confidence: medium)

### 7. `native` — ~9 LOC — `Source Code/hamiltonian_builders.py:387`

- **Suggestion:** hand-rolled pure-Python matrix algebra in hamiltonian_builders.py: matmult (zip-comprehension matrix product) and scalarmult are never called (matmult only in a commented line and in notebooks' own local copies); matadd is called 5 times on operands that are ALREADY numpy arrays being .tolist()'d just to feed it
- **Replacement:** delete matmult/scalarmult; replace matadd(H0, (X).tolist()) with numpy '+' on the object-dtype array (np is already imported and computes every operand via '@')
- **Auditor evidence:** Read L385-399 (defs) and grepped call sites: matadd used at L77,79,81,165,167 in the form matadd(H0,(-params['D']*N0@N0).tolist()); matmult/scalarmult appear elsewhere only in commented L82 and in *.ipynb cells that define their own local copies (YbOH_X010_AOM_OBEs_N1.ipynb:3027).
- **Skeptic verdict:** confirmed — matmult/scalarmult: zero live callers — matmult only in commented L82 and its own def; notebook hits (e.g. YbOH_X010_AOM_OBEs_N1.ipynb:3027) are notebook-local 'def matmult' copies, and the only hamiltonian_builders import in any notebook is build_operator. matadd call sites verified: H0 is np.zeros(...).tolist() filled with float(...) entries, operands are float ndarrays, so (np.asarray(H0)+X).tolist() is behavior-identical (keeping .tolist() sidesteps the file's 'Sympy does not like numpy arrays' comment). One evidence inaccuracy: L165/167 pass the ndarray WITHOUT .tolist() (matadd zips it fine) — doesn't change the conclusion. Defs span L387-396 ≈ 9 lines. (auditor confidence: medium)

### 8. `delete` — ~9 LOC — `Source Code/Energy_Levels.py:693`

- **Suggestion:** MoleculeLevels.filter_evecs — 9-line eigenvector filter with zero callers anywhere
- **Replacement:** delete outright
- **Auditor evidence:** grep 'filter_evecs' over all *.ipynb: zero hits; over Energy_Levels.py: only the def at L693 (no self.filter_evecs anywhere). Neighboring select_q has 59 notebook callers, so the zero is signal not grep failure.
- **Skeptic verdict:** confirmed — filter_evecs is exactly 9 lines (Energy_Levels.py:693-701). Repo-wide grep (py + ipynb): the def line is the sole occurrence — no self.filter_evecs, no notebook call. Neighboring methods do get notebook use, so the zero is signal. (auditor confidence: high)

### 9. `delete` — ~8 LOC — `Source Code/matrix_elements.py:458`

- **Suggestion:** MQM_bBS_old — superseded 'combined tensor form' variant of MQM_bBS kept beside the active one
- **Replacement:** delete outright
- **Auditor evidence:** grep 'MQM_bBS_old' over all *.ipynb and Source Code *.py excluding the def: zero callers (hamiltonian_builders imports MQM_bBS, not _old). Kept-alternate-form cross-check is plausible for a physicist (comment '#Using combined tensor form'), hence low confidence per the alternate-paths rule.
- **Skeptic verdict:** confirmed — grep 'MQM_bBS_old' repo-wide: only the two defs (matrix_elements.py:458 and a twin at matrix_elements_sym.py:280 the finding under-scopes) — zero callers. hamiltonian_builders.py:5 imports the active MQM_bBS and uses it at L289/293. The _old def spans 458-466 = 9 lines. 'Kept alternate tensor form' is a preference, not a usage; git history retains it. (auditor confidence: low)

### 10. `stdlib` — ~6 LOC — `Source Code/matrix_elements.py:24`

- **Suggestion:** hand-rolled 5-line kronecker delta (if a==b: return 1 else: return 0), duplicated in matrix_elements.py and matrix_elements_sym.py
- **Replacement:** def kronecker(a,b): return int(a==b) — or inline (a==b) since callers only multiply it
- **Auditor evidence:** Read both defs (matrix_elements.py:24, matrix_elements_sym.py:17); call sites are guard products like kronecker(K0,K1)*kronecker(N0,N1)*... throughout both files, where a bool/int works identically.
- **Skeptic verdict:** confirmed — Both defs read (matrix_elements.py:24, matrix_elements_sym.py:17; a third copy lives in case_bBJ_MEs.py, moot if finding 0 lands). 130 call sites in matrix_elements.py alone, in product-guard, 'if kronecker(...)', and 'return kronecker(...)*...' forms — int(a==b) returns the identical 0/1 with identical truthiness and arithmetic in every context; args are scalar quantum numbers, no array inputs. Fully behavior-identical. (auditor confidence: high)

### 11. `delete` — ~1 LOC — `Source Code/molecule_parameters.py:3`

- **Suggestion:** unused 'from sympy.physics.wigner import wigner_3j,wigner_6j,wigner_9j' at the top of the constants file
- **Replacement:** delete outright
- **Auditor evidence:** grep -n 'wigner|kronecker' in molecule_parameters.py: only the import line itself matches; the file is a constants registry with no matrix-element math.
- **Skeptic verdict:** confirmed — grep in molecule_parameters.py: the wigner import at L3 is the only wigner match in the file. Checked re-export risk: no 'from molecule_parameters import *' anywhere; all consumers import get_molecule_params/params_general/c by name; the one 'import molecule_parameters as _mp' consumer (RaX/zeeman_slower_2.ipynb) accesses only _mp.molecules. Safe 1-line deletion. (auditor confidence: high)

## Rejected by adversarial verification (recorded so they are not re-flagged)

- `Jupyter Notebooks/config_path.py` — config_path.py sys.path-walking bootstrap copy-pasted into 9 notebook directories → **refuted:** Direct counter-evidence: 'Jupyter Notebooks/CLAUDE.md':27 states verbatim 'Each subdir has its own config_path.py (identical files). Keep the import; don't try to package-ify it.' — a documented deliberate project rule forbidding exactly this replacement. Independently, the finding's own plan ('leave existing notebooks untouched') forces all 9 copies to stay on disk since notebook first-cells do 'from config_path import add_to_sys_path'; deleting them breaks every notebook, so the claimed 120 LOC is unrealizable (actual in-repo saving: 0). A conda-env .pth is also machine/env-local, while the in-repo bootstrap survives cloning and env recreation on the user's two machines.

## Cross-repo notes

- Could adopt C2V's cached `wigner_3j/6j/9j` + `delta` block into `matrix_elements.py` (same names/signatures, pick-one-copy, no shared package). STOP for Arian first: Wigner/phase conventions.

## Verification recipe

- `conda run -n structure python -c "import sys; sys.path.insert(0,'Source Code'); import Energy_Levels, matrix_elements, quantum_numbers, hamiltonian_builders"`
- execute `YbOH_Branching.ipynb` via nbconvert, NO --allow-errors, stop at first error
- any change near eigenvector ordering/phase: small ZeemanMap smoke sweep, compare pre/post eigenvalues + overlaps (np.allclose)
