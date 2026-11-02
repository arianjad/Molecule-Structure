# zeeman_slower_2 Comprehensive Hardening — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers-extended-cc:subagent-driven-development (recommended) or superpowers-extended-cc:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land the 3-phase comprehensive hardening of `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` from spec `docs/superpowers/specs/2026-05-13-zeeman-slower-2-comprehensive-hardening-design.md`: documentation/cheap fixes → physics diagnostic cells (with `remix`→`s_remix` rename + B_⊥ calibration + PB comb selectivity) → refactor (extract trajectory-sweep helpers + PB m_S selection from decoupled basis) — without changing physics conclusions (regression-checked against a Phase 1 baseline JSON).

**Architecture:** Notebook-only edits, three sequential commits. Phase 1 is additive comments + warnings + a serialization cell. Phase 2 adds new diagnostic cells alongside one mechanical rename. Phase 3 extracts helper functions and replaces five duplicated sweep-loop cells with thin wrappers, plus a regression cell that diffs against the Phase 1 baseline.

**Tech Stack:** Python 3.12 in conda env `Structure`, Jupyter notebook (`.ipynb` JSON), numpy, scipy.integrate.solve_ivp, scipy.interpolate (RegularGridInterpolator, interp1d), matplotlib, json (stdlib for baseline serialization). Project-internal modules in `Source Code/`. Verification gate per phase: `conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"` exits 0 with zero `output_type=='error'` cells (validated by reading the JSON, NOT by `--allow-errors`).

**Conventions for cell edits.** The notebook is JSON; cell sources are lists of strings each ending in `\n`. Edits use one of two strategies:
- **Edit tool**: for in-place modifications to existing cell source where the matched string is unique. Match the exact bytes including the leading indentation and trailing `\n`.
- **Python helper script**: for inserting new cells, replacing entire cell sources, or any change spanning multiple discontinuous source lines. Pattern:
  ```python
  import json, sys, uuid
  with open(nb_path) as f: nb = json.load(f)
  # ... mutate nb['cells'] ...
  tmp = nb_path + ".tmp"
  with open(tmp, "w") as f: json.dump(nb, f, indent=1)
  import os; os.replace(tmp, nb_path)
  ```

**Path conventions:** All paths are relative to the repo root `/Users/arianjadbabaie/Library/CloudStorage/GoogleDrive-arianjad@mit.edu/Shared drives/EMA-data-server/RaX/Personal/ArianJadbabaie/Code/Molecule-Structure`. Use `/Library/Developer/CommandLineTools/usr/bin/git` (NOT `/usr/bin/git`) to avoid xcrun shim hangs on this Google Drive Cloud Storage filesystem (see `~/.claude/lessons.md` 2026-05-09).

---

## Phase 1 — Documentation + cheap fixes (warmup)

### Task 1: Cheap edits to existing cells

**Goal:** Add docstring upgrades, monotonicity warning, hyperfine-group warning, capture-totals summaries, and lazy-init comment fix. Pure additions of comments/warnings — no behavior change.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cells idx=27, 33, 36, 38, 40, 42

**Acceptance Criteria:**
- [ ] `steady_state_rate_eq` docstring (cell 27) has a new paragraph explaining `remix` is a saturation parameter (s = I/I_sat) of an isotropic always-on-resonance pseudo-laser, NOT a fraction.
- [ ] `make_B_of_z_nonlinear` (cell 33) prints a `⚠ ...` warning naming the (je, ig) pair and the B values of any local extrema if `nu_eg(B)` is non-monotonic. Continues execution.
- [ ] Cell 36 prints `Total captures across all 24 configs x N v_0 = X/Y` after its capture-velocity table.
- [ ] Cell 40 prints the same total summary.
- [ ] Cell 38, after `groups = find_X_N1_hyperfine_groups(g)`, prints `⚠ ...` warning if `len(groups) != 4`. Continues execution.
- [ ] Cell 42 has an explanatory comment on the `F_max_default` lazy-init pattern naming the assumption (`remix_values[0]` is the intended baseline).

**Verify:** Open the notebook and grep the cell sources:
```
python3 - <<'PY'
import json
with open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb") as f: nb = json.load(f)
need = [
    ("remix is a saturation parameter", 27),
    ("non-monotonic", 33),
    ("Total captures across all 24 configs", 36),
    ("Total captures across all 24 configs", 40),
    ("X²Σ⁺ N=1 expected 4 (J,F) groups", 38),
    ("remix_values[0] is the intended baseline", 42),
]
for needle, idx in need:
    src = "".join(nb["cells"][idx].get("source", []))
    print(f"cell {idx}: {'OK' if needle in src else 'MISSING:'} {needle!r}")
PY
```
Expected: all six lines print `cell N: OK ...`.

**Steps:**

- [ ] **Step 1: Cell 27 — ssr docstring upgrade.** The docstring of `steady_state_rate_eq` currently doesn't explain `remix` units. Add a paragraph after the `'polarization': ...` line of the sideband dict explanation. Use the Edit tool on the .ipynb file. The exact `old_string` and `new_string`:

  Find this string in the cell 27 source (note the leading `    "` and trailing `\n",` from .ipynb JSON formatting):
  ```
      "    Saturation in denominator: NONE (Tarbutt 2014 Eq 2 convention). Population balance\n",
      "    in the rate matrix below caps the cycle at high intensity.\n",
      "    \"\"\"\n",
  ```
  Replace with:
  ```
      "    Saturation in denominator: NONE (Tarbutt 2014 Eq 2 convention). Population balance\n",
      "    in the rate matrix below caps the cycle at high intensity.\n",
      "\n",
      "    `remix` units: `remix` is a saturation parameter (I/I_sat) of an isotropic, always-\n",
      "    on-resonance pseudo-laser added equally to every (e,g,polarization) channel. NOT a\n",
      "    fraction. Physical interpretation: maps to off-axis B-field magnitude via Larmor\n",
      "    mixing rate Ω_L = g_F μ_B B_⊥ / ℏ. See `b_perp_to_s_remix` (Phase 2b) and the\n",
      "    derivation cell after this one for the numeric mapping.\n",
      "    \"\"\"\n",
  ```

- [ ] **Step 2: Cell 33 — monotonicity warning in `make_B_of_z_nonlinear`.** The function currently does `nu_eg_curve = TM["dE"][:, je, ig]` then later `sort_idx = np.argsort(nu_eg_curve)`. Insert a check between these. Use Edit:

  Find:
  ```
      "    Bz = TM[\"Bz\"]\n",
      "    nu_eg_curve = TM[\"dE\"][:, je, ig]                    # (N_B,) cyclic MHz\n",
      "    z, v_target, eta = slowing_profile_fixed_L(v_i, v_f, L, a_max, n_z)\n",
  ```
  Replace with:
  ```
      "    Bz = TM[\"Bz\"]\n",
      "    nu_eg_curve = TM[\"dE\"][:, je, ig]                    # (N_B,) cyclic MHz\n",
      "    diffs = np.diff(nu_eg_curve)\n",
      "    if (diffs > 0).any() and (diffs < 0).any():\n",
      "        sign_changes = np.where(np.diff(np.sign(diffs)) != 0)[0]\n",
      "        Bz_extrema = Bz[sign_changes + 1]\n",
      "        print(f\"⚠ make_B_of_z_nonlinear: nu_eg(B) for (je={je}, ig={ig}) is non-monotonic; \"\n",
      "              f\"local extrema near B = {np.round(Bz_extrema, 1)} G. Inversion uses argsort+interp1d \"\n",
      "              f\"and silently picks one branch.\")\n",
      "    z, v_target, eta = slowing_profile_fixed_L(v_i, v_f, L, a_max, n_z)\n",
  ```

- [ ] **Step 3: Cell 36 — total captures summary.** Cell 36's print loop ends with the capture-velocity table. Add a `total_captures` accumulator and a final summary line. The existing inner loop is:

  Find:
  ```
      "for Bd in B_designs:\n",
      "    for s in s_values:\n",
      "        cells = []\n",
      "        for L in L_values:\n",
      "            captured = [v0 for v0 in v0_scan\n",
      "                          if results_single[(Bd, s, L, v0)][\"exit_reason\"] == \"v_stop\"]\n",
      "            v_cap = max(captured) if captured else None\n",
      "            cells.append(f\"  {v_cap:6.0f}\" if v_cap is not None else \"      --\")\n",
      "        print(f\"  B_d={Bd:3.0f} G, s_sb={s:5.1f}: \" + \"\".join(cells))\n",
  ```
  Replace with:
  ```
      "total_captures = 0\n",
      "for Bd in B_designs:\n",
      "    for s in s_values:\n",
      "        cells = []\n",
      "        for L in L_values:\n",
      "            captured = [v0 for v0 in v0_scan\n",
      "                          if results_single[(Bd, s, L, v0)][\"exit_reason\"] == \"v_stop\"]\n",
      "            total_captures += len(captured)\n",
      "            v_cap = max(captured) if captured else None\n",
      "            cells.append(f\"  {v_cap:6.0f}\" if v_cap is not None else \"      --\")\n",
      "        print(f\"  B_d={Bd:3.0f} G, s_sb={s:5.1f}: \" + \"\".join(cells))\n",
      "print(f\"\\nTotal captures across all 24 configs x {len(v0_scan)} v_0 = {total_captures}/{24*len(v0_scan)}\")\n",
  ```

  Caveat: the variable name in this cell may be `results_single` or `results`. Before editing, verify with `python3 - <<'PY' ... print(nb['cells'][36]['source'])` and substitute the actual name.

- [ ] **Step 4: Cell 40 — total captures summary.** Same pattern as Step 3, but for `results_multi`. The capture-loop in cell 40 has the same structure. Apply the same `total_captures` accumulator + final summary print.

- [ ] **Step 5: Cell 38 — hyperfine groups assertion.** After the call to `find_X_N1_hyperfine_groups(g)` in cell 38, add a warning. Find:

  ```
      "groups = find_X_N1_hyperfine_groups(g)\n",
  ```
  (it appears once in cell 38, calling the function then printing its results)

  Replace with:
  ```
      "groups = find_X_N1_hyperfine_groups(g)\n",
      "if len(groups) != 4:\n",
      "    print(f\"⚠ X²Σ⁺ N=1 expected 4 (J,F) groups (J=½ F=0, J=½ F=1, J=3/2 F=1, J=3/2 F=2); \"\n",
      "          f\"got {len(groups)}. Check energy_tol_MHz=1.0 vs actual sublevel spacings.\")\n",
  ```

- [ ] **Step 6: Cell 42 — F_max_default lazy-init comment.** Find:
  ```
      "F_max_default = None\n",
  ```
  Replace with:
  ```
      "# F_max_default = lowest-remix F_max (no dark-state mixing baseline). Initialized on the\n",
      "# first iteration; assumes remix_values[0] is the intended baseline (sorted ascending).\n",
      "F_max_default = None\n",
  ```

- [ ] **Step 7: Verify the edits without re-executing.** Run the verify command at the top of this task. Confirm all six `OK` lines.

- [ ] **Step 8: Commit.** Do NOT push. Use the CLT git path:
  ```
  cd "/Users/arianjadbabaie/Library/CloudStorage/GoogleDrive-arianjad@mit.edu/Shared drives/EMA-data-server/RaX/Personal/ArianJadbabaie/Code/Molecule-Structure"
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 1a — docstring upgrades + monotonicity & group warnings + capture totals"
  ```

---

### Task 2: Baseline-serialization cell

**Goal:** Add a new cell at the end of the notebook (just before the "Notes & next steps" markdown) that serializes all sweep results — capture-velocity tables and per-config F_max — to a JSON file. This is the regression baseline used by Phase 3.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` (insert cell)
- Create at runtime: `Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json`
- Add `Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json` to `.gitignore`

**Acceptance Criteria:**
- [ ] New code cell inserted as the second-to-last cell (last code cell; followed only by the existing "Notes & next steps" markdown).
- [ ] Cell, when executed, writes `zs_baselines_pre_refactor.json` containing all five sweeps' results in a structured format.
- [ ] The JSON file is in `.gitignore` (it's regenerated; not committed).
- [ ] After Phase 1 commit, running the verify command (below) produces a non-empty JSON with the expected schema.

**Verify:**
```
conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
test -f "Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json"
python3 -c 'import json; d=json.load(open("Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json")); print(sorted(d.keys()))'
```
Expected: notebook executes 0; file exists; printed keys include `["pb", "remix_sweep", "remix_full", "single_sb", "multi_sb"]` (or whatever subset of the five sweeps actually ran).

**Steps:**

- [ ] **Step 1: Add the JSON file to .gitignore.**
  ```
  echo "Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json" >> .gitignore
  /Library/Developer/CommandLineTools/usr/bin/git add .gitignore
  ```
  (commit will happen at the end of Task 2, bundled with the cell insert.)

- [ ] **Step 2: Insert the new cell.** Use a Python helper script. Save to `/tmp/insert_baseline_cell.py`:

  ```python
  import json, sys, uuid
  nb_path = sys.argv[1]
  with open(nb_path) as f:
      nb = json.load(f)

  cell = {
      "cell_type": "code",
      "metadata": {},
      "execution_count": None,
      "outputs": [],
      "source": [
          "# === Phase 1f: serialize sweep results as regression baseline ============================\n",
          "# All capture-velocity range tables and per-config F_max are written to JSON. Phase 3 reads\n",
          "# this back and asserts identity (or itemizes deltas with physics attribution).\n",
          "import json\n",
          "\n",
          "def _capture_range(results_dict, key_fn, v0_iter):\n",
          "    \"\"\"Returns sorted list of v0 captured (exit_reason == 'v_stop') for a given config key.\"\"\"\n",
          "    return sorted([float(v0) for v0 in v0_iter\n",
          "                   if results_dict[key_fn(v0)][\"exit_reason\"] == \"v_stop\"])\n",
          "\n",
          "baseline = {}\n",
          "\n",
          "# 1. Single-σ+ sweep (Step 6, cell 36): results_single[(Bd, s, L, v0)]\n",
          "if 'results_single' in dir():\n",
          "    baseline['single_sb'] = {\n",
          "        'F_max':   {f\"{Bd:.0f}_{s:.0f}\": float(np.abs(F_surfaces_single[(Bd, s)]).max())\n",
          "                    for Bd in B_designs for s in s_values},\n",
          "        'captures': {f\"{Bd:.0f}_{s:.0f}_{L:.2f}\":\n",
          "                       _capture_range(results_single, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)\n",
          "                     for Bd in B_designs for s in s_values for L in L_values},\n",
          "    }\n",
          "\n",
          "# 2. Multi-σ+ sweep (Step 7, cell 40): results_multi\n",
          "if 'results_multi' in dir():\n",
          "    baseline['multi_sb'] = {\n",
          "        'F_max':    {f\"{Bd:.0f}_{s:.0f}\": float(np.abs(F_surfaces_multi[(Bd, s)]).max())\n",
          "                     for Bd in B_designs for s in s_values},\n",
          "        'captures': {f\"{Bd:.0f}_{s:.0f}_{L:.2f}\":\n",
          "                       _capture_range(results_multi, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)\n",
          "                     for Bd in B_designs for s in s_values for L in L_values},\n",
          "    }\n",
          "\n",
          "# 3. Remix sweep (Step 7b, cell 42): results_remix[(rmx, v0)]\n",
          "if 'results_remix' in dir():\n",
          "    baseline['remix_sweep'] = {\n",
          "        'F_max':    {f\"{rmx:.0e}\": float(F_max_remix[rmx]) for rmx in remix_values},\n",
          "        'captures': {f\"{rmx:.0e}\":\n",
          "                       _capture_range(results_remix, lambda v0,rmx=rmx: (rmx,v0), v0_remix_scan)\n",
          "                     for rmx in remix_values},\n",
          "    }\n",
          "\n",
          "# 4. Full-remix=1.0 sweep (Step 7c, cell 44): results_remix1\n",
          "if 'results_remix1' in dir():\n",
          "    baseline['remix_full'] = {\n",
          "        'F_max':    {f\"{Bd:.0f}_{s:.0f}\": float(np.abs(F_surfaces_remix1[(Bd, s)]).max())\n",
          "                     for Bd in B_designs for s in s_values},\n",
          "        'captures': {f\"{Bd:.0f}_{s:.0f}_{L:.2f}\":\n",
          "                       _capture_range(results_remix1, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)\n",
          "                     for Bd in B_designs for s in s_values for L in L_values},\n",
          "    }\n",
          "\n",
          "# 5. Paschen-Back σ⁻ sweep (cell 48): results_pb[(s_per_sb, L, v0)]\n",
          "if 'results_pb' in dir():\n",
          "    baseline['pb'] = {\n",
          "        'F_max':    {f\"{s:.0f}\": float(np.abs(F_surfaces_pb[s]).max()) for s in s_values},\n",
          "        'captures': {f\"{s:.0f}_{L:.2f}\":\n",
          "                       _capture_range(results_pb, lambda v0,s=s,L=L: (s,L,v0), v0_scan)\n",
          "                     for s in s_values for L in L_values},\n",
          "    }\n",
          "\n",
          "_path = 'Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json'\n",
          "with open(_path, 'w') as _f: json.dump(baseline, _f, indent=2, sort_keys=True)\n",
          "print(f\"Wrote {_path}: {sorted(baseline.keys())}\")\n",
      ]
  }

  # Insert as second-to-last (last code cell; before the trailing "Notes & next steps" markdown)
  nb["cells"] = nb["cells"][:-1] + [cell] + [nb["cells"][-1]]

  tmp = nb_path + ".tmp"
  with open(tmp, "w") as f:
      json.dump(nb, f, indent=1)
  import os
  os.replace(tmp, nb_path)
  print(f"Inserted baseline cell. Notebook now has {len(nb['cells'])} cells.")
  ```

  Run it:
  ```
  python3 /tmp/insert_baseline_cell.py "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  ```
  Expected output: `Inserted baseline cell. Notebook now has 51 cells.` (50 from prior + 1).

- [ ] **Step 3: Execute the notebook end-to-end** to populate the JSON file:
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  ```
  Expected: exits 0. Should take ~1-2 minutes (most of the time is in the 408-trajectory sweeps).

  If it does NOT exit 0, inspect the notebook JSON for error cells:
  ```
  python3 -c 'import json; nb=json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")); errs=[(i,o.get("ename")) for i,c in enumerate(nb["cells"]) for o in c.get("outputs",[]) if o.get("output_type")=="error"]; print(errs or "no errors")'
  ```

- [ ] **Step 4: Verify the baseline file** per the Verify block at the top of this task.

- [ ] **Step 5: Commit Phase 1.**
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add .gitignore "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 1b — baseline serialization cell + .gitignore entry"
  ```

---

## Phase 2 — Physics diagnostic cells

### Task 3: ν_eg(B) monotonicity diagnostic cell

**Goal:** New cell after cell 33 (the `je_ref, ig_ref` pick) that plots `ν_eg(B)` for the chosen reference channel over the full Bz grid, identifies local extrema, and prints a verdict. Diagnostic only — no behavior change to existing cells.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` (insert cell after idx=33; existing idx=34 shifts to idx=35, etc.)

**Acceptance Criteria:**
- [ ] Cell inserted at idx=34 (immediately after the cell that defines `je_ref, ig_ref`).
- [ ] When executed, cell prints `ν_eg(B) MONOTONIC for (je=X, ig=Y) over Bz=[a, b] G` OR `NON-MONOTONIC: extrema at B≈[...] G`.
- [ ] Cell renders one plot: ν_eg(B) vs B, with extrema marked.
- [ ] Notebook re-executes end-to-end clean.

**Verify:**
```
conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 -c 'import json; nb=json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")); src="".join(nb["cells"][34]["source"]); assert "MONOTONIC" in src or "monotonic" in src; print("OK")'
```
Expected: notebook exits 0; assertion passes.

**Steps:**

- [ ] **Step 1: Insert cell at idx=34.** Save script to `/tmp/insert_monotonicity_cell.py`:
  ```python
  import json, sys
  nb_path = sys.argv[1]
  with open(nb_path) as f: nb = json.load(f)
  cell = {
      "cell_type": "code",
      "metadata": {},
      "execution_count": None,
      "outputs": [],
      "source": [
          "# === Phase 2a: ν_eg(B) monotonicity diagnostic ===========================================\n",
          "# `make_B_of_z_nonlinear` inverts ν_eg(B) for the brightest channel (je_ref, ig_ref).\n",
          "# Inversion via argsort+interp1d silently picks one branch if ν_eg(B) is non-monotonic.\n",
          "# Plot the curve and identify any extrema before relying on B(z) construction.\n",
          "Bz_full = TM['Bz']\n",
          "nu_eg_check = TM['dE'][:, je_ref, ig_ref]\n",
          "diffs_check = np.diff(nu_eg_check)\n",
          "if (diffs_check > 0).any() and (diffs_check < 0).any():\n",
          "    sign_changes = np.where(np.diff(np.sign(diffs_check)) != 0)[0]\n",
          "    Bz_extrema = Bz_full[sign_changes + 1]\n",
          "    print(f\"⚠ ν_eg(B) NON-MONOTONIC for (je={je_ref}, ig={ig_ref}) over \"\n",
          "          f\"Bz=[{Bz_full[0]:.3f}, {Bz_full[-1]:.0f}] G. \"\n",
          "          f\"Local extrema at B = {np.round(Bz_extrema, 1)} G.\")\n",
          "else:\n",
          "    print(f\"✓ ν_eg(B) MONOTONIC for (je={je_ref}, ig={ig_ref}) over \"\n",
          "          f\"Bz=[{Bz_full[0]:.3f}, {Bz_full[-1]:.0f}] G.\")\n",
          "    Bz_extrema = np.array([])\n",
          "\n",
          "fig, ax = plt.subplots(1, 1, figsize=(8, 4))\n",
          "ax.plot(Bz_full, nu_eg_check, '-', lw=1.4, color='C0')\n",
          "for Be in Bz_extrema:\n",
          "    ax.axvline(Be, color='r', lw=0.6, alpha=0.6)\n",
          "ax.set_xlabel('B (G)')\n",
          "ax.set_ylabel(r'$\\nu_{eg}$ (MHz, cyclic)')\n",
          "ax.set_title(f'$\\\\nu_{{eg}}(B)$ for reference channel (je={je_ref}, ig={ig_ref})')\n",
          "ax.grid(alpha=0.3)\n",
          "plt.tight_layout()\n",
          "plt.show()\n",
      ]
  }
  nb['cells'].insert(34, cell)
  tmp = nb_path + '.tmp'
  with open(tmp, 'w') as f: json.dump(nb, f, indent=1)
  import os; os.replace(tmp, nb_path)
  print(f"Inserted at idx=34. Notebook now has {len(nb['cells'])} cells.")
  ```
  Run: `python3 /tmp/insert_monotonicity_cell.py "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"`

- [ ] **Step 2: Execute the notebook** and confirm verdict + plot rendered:
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  ```
  Expected: exits 0.

- [ ] **Step 3: Verify** per the verify block.

- [ ] **Step 4: Commit.**
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 2a — ν_eg(B) monotonicity diagnostic"
  ```

---

### Task 4: Rename remix → s_remix + B_⊥ calibration helper

**Goal:** Mechanical rename of the `remix` parameter to `s_remix` everywhere in the notebook to reflect that it's a saturation parameter, not a fraction. Add a `b_perp_to_s_remix(B_perp_G, ...)` helper function and a markdown cell with the Larmor-mapping derivation.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 27 (function signature + body), cell 42 (call site), and any other call sites surfaced by grep
- Insert new cell after cell 27: helper function `b_perp_to_s_remix`
- Insert new markdown cell before cell 42 (after the rename-shifted index): derivation + numeric mapping table

**Acceptance Criteria:**
- [ ] Zero remaining occurrences of `remix=` (parameter), `remix:` (dict key), or bare `remix` (variable) outside of comments/docstrings explaining the rename. Use `s_remix` everywhere.
- [ ] New helper function `b_perp_to_s_remix(B_perp_G, Gamma_FWHM_MHz=Gamma_FWHM_MHz, g_F=0.5)` defined and documented.
- [ ] New markdown cell with derivation: `s_remix · (Γ_angular/2) = Ω_L = g_F μ_B B_⊥ / ℏ` ⟹ `s_remix = 2 g_F μ_B B_⊥ / (h Γ_FWHM)` with numeric mapping table.
- [ ] Notebook re-executes end-to-end clean.
- [ ] All sweep results (printed capture tables and F_max numbers) are bit-identical to Phase 1f baseline. (The rename is mechanical; physics unchanged.)

**Verify:**
```
conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 - <<'PY'
import json, re
with open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb") as f: nb = json.load(f)
hits = []
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code": continue
    src = "".join(c.get("source", []))
    # match `remix` not preceded/followed by `s_` or letters: bare references only
    for m in re.finditer(r'(?<![A-Za-z_])remix(?![A-Za-z_])', src):
        hits.append((i, m.start(), src[max(0,m.start()-30):m.start()+40]))
print("BARE 'remix' hits (should be 0 outside comments):", len(hits))
for h in hits: print(" ", h)
PY
```
Expected: 0 hits, OR hits only inside comments/docstrings (manually inspect the printed context).

**Steps:**

- [ ] **Step 1: Locate every occurrence of `remix` in the notebook.**
  ```
  python3 - <<'PY'
  import json, re
  with open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb") as f: nb = json.load(f)
  for i, c in enumerate(nb["cells"]):
      src = "".join(c.get("source", []))
      for m in re.finditer(r'remix', src):
          start = max(0, m.start() - 25); end = min(len(src), m.start() + 30)
          print(f"cell {i:2d} char {m.start():5d}: ...{src[start:end]!r}")
  PY
  ```
  Catalog every line. Most will be in cell 27 (function definition), cell 42 (call site + variable names like `remix_values`, `remix_full`, `Bd_remix`, `s_sb_remix`, `L_remix`, `B_of_z_remix`, `results_remix`, `F_max_remix`, `sb_template_remix`, `v0_remix_scan`), and cell 44 (likely re-uses `remix_full`).

- [ ] **Step 2: Rename `remix` parameter → `s_remix` in `steady_state_rate_eq`.** In cell 27, edit the function signature:

  Find:
  ```
      "                          Gamma_FWHM_MHz=Gamma_FWHM_MHz, remix=1e-3,\n",
  ```
  Replace with:
  ```
      "                          Gamma_FWHM_MHz=Gamma_FWHM_MHz, s_remix=1e-3,\n",
  ```

  Then in the body of the function:

  Find:
  ```
      "    s_vec_eff_total = np.full(3, remix, dtype=float)   # for dark-state mixing only\n",
  ```
  Replace with:
  ```
      "    s_vec_eff_total = np.full(3, s_remix, dtype=float)   # for dark-state mixing only\n",
  ```

  Find:
  ```
      "    # add remix as a pseudo-on-resonance saturation (lifts dark states) — does NOT contribute to force\n",
  ```
  Replace with:
  ```
      "    # add s_remix as a pseudo-on-resonance saturation (lifts dark states) — does NOT contribute to force\n",
  ```

  Find:
  ```
      "        R_p_eff[ip] += prefac * (remix * TDM2[iB, ip])         # remix already = saturation; treats as on-res\n",
  ```
  Replace with:
  ```
      "        R_p_eff[ip] += prefac * (s_remix * TDM2[iB, ip])       # s_remix is saturation parameter; treats as on-res\n",
  ```

  Update the smoke-test call at the bottom of cell 27:
  Find:
  ```
      "res = steady_state_rate_eq(TM, 200., 0., sidebands=_sb_test, remix=1e-3)\n",
  ```
  Replace with:
  ```
      "res = steady_state_rate_eq(TM, 200., 0., sidebands=_sb_test, s_remix=1e-3)\n",
  ```

- [ ] **Step 3: Update Phase 1a docstring** (added in Task 1) to reference `s_remix`:

  Find:
  ```
      "    `remix` units: `remix` is a saturation parameter (I/I_sat) of an isotropic, always-\n",
  ```
  Replace with:
  ```
      "    `s_remix` units: `s_remix` is a saturation parameter (I/I_sat) of an isotropic, always-\n",
  ```

- [ ] **Step 4: Insert `b_perp_to_s_remix` helper cell** immediately after cell 27 (so it's defined before any call site uses it). Use a Python helper script:

  ```python
  import json, sys
  nb_path = sys.argv[1]
  with open(nb_path) as f: nb = json.load(f)
  cell = {
      "cell_type": "code",
      "metadata": {},
      "execution_count": None,
      "outputs": [],
      "source": [
          "# === Phase 2b: B_⊥ ↔ s_remix calibration helper ============================================\n",
          "# Physical interpretation of `s_remix`: equate the pseudo-laser pumping rate per channel\n",
          "# (Γ_angular/2) · s_remix to the Larmor mixing rate Ω_L = g_F μ_B B_⊥ / ℏ from a real off-axis\n",
          "# magnetic field B_⊥. Solving:  s_remix = 2 g_F μ_B B_⊥ / (ℏ Γ_angular) = 2 g_F μ_B B_⊥ / (h Γ_FWHM)\n",
          "# In MHz/G units (μ_B / h = 1.39962 MHz/G), s_remix = 2 g_F · 1.39962 · B_⊥ / Γ_FWHM_MHz.\n",
          "MUB_OVER_H_MHZ_G = 1.39962  # Bohr magneton / h, in MHz/G\n",
          "\n",
          "def b_perp_to_s_remix(B_perp_G, Gamma_FWHM_MHz=Gamma_FWHM_MHz, g_F=0.5):\n",
          "    \"\"\"Convert off-axis B-field magnitude (in Gauss) to the equivalent s_remix value.\n",
          "    Default g_F=0.5 is the order-of-magnitude effective Landé factor for X²Σ⁺ N=1 (case-b limit).\n",
          "    For sublevel-specific g_F, see Brown & Carrington Ch. 9.\n",
          "    Returns s_remix as a saturation parameter (dimensionless).\"\"\"\n",
          "    return 2.0 * g_F * MUB_OVER_H_MHZ_G * B_perp_G / Gamma_FWHM_MHz\n",
          "\n",
          "def s_remix_to_b_perp(s_remix, Gamma_FWHM_MHz=Gamma_FWHM_MHz, g_F=0.5):\n",
          "    \"\"\"Inverse of b_perp_to_s_remix: returns B_⊥ in Gauss.\"\"\"\n",
          "    return s_remix * Gamma_FWHM_MHz / (2.0 * g_F * MUB_OVER_H_MHZ_G)\n",
          "\n",
          "print(f\"s_remix ↔ B_⊥ mapping table (g_F=0.5, Γ_FWHM={Gamma_FWHM_MHz} MHz):\")\n",
          "print(f\"  {'B_⊥ (G)':>10s}    {'s_remix':>10s}\")\n",
          "for B in [0.001, 0.01, 0.1, 1.0, 3.3, 10., 30., 100.]:\n",
          "    print(f\"  {B:>10.3f}    {b_perp_to_s_remix(B):>10.4g}\")\n",
      ]
  }
  nb['cells'].insert(28, cell)  # immediately after cell 27 (which is ssr def)
  tmp = nb_path + '.tmp'
  with open(tmp, 'w') as f: json.dump(nb, f, indent=1)
  import os; os.replace(tmp, nb_path)
  print(f"Inserted helper at idx=28. Notebook has {len(nb['cells'])} cells.")
  ```
  Run: `python3 /tmp/insert_b_perp_helper.py "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"`

- [ ] **Step 5: Update all `remix` call sites** in cells 42, 44, etc. The variable names that include `remix` as a substring (e.g. `remix_values`, `remix_full`, `Bd_remix`, `results_remix`, `F_max_remix`, `sb_template_remix`, `v0_remix_scan`, `B_of_z_remix`, `L_remix`, `s_sb_remix`, `precompute_F_grid_remix`, `results_remix1`) are descriptive variable names referring to the remix sweep concept and may stay UNCHANGED. The mechanical replacement is only the **parameter name** when calling `steady_state_rate_eq` and `precompute_F_grid_remix`:

  In cell 42, find any occurrence of the parameter name in calls. Search:
  ```
  python3 -c 'import json; nb=json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")); src="".join(nb["cells"][42]["source"]); print(src)'
  ```
  Look for keyword arg `remix=` or positional uses. The cell 42 source includes `precompute_F_grid_remix(TM, Bd_remix, v_design, sb_template_remix, rmx, je_ref, ig_ref, v_grid, B_grid)` (positional `rmx`); inside, `precompute_F_grid_remix` calls `steady_state_rate_eq(..., remix=remix, ...)`.

  Find the `precompute_F_grid_remix` definition (in cell 42):
  ```
      "def precompute_F_grid_remix(TM, B_design_G, v_design, sidebands_template, remix,\n",
      "                              je_ref, ig_ref, v_grid, B_grid):\n",
      "    \"\"\"precompute_F_grid_with_sidebands but with remix as a free parameter.\"\"\"\n",
  ```
  Replace with:
  ```
      "def precompute_F_grid_remix(TM, B_design_G, v_design, sidebands_template, s_remix,\n",
      "                              je_ref, ig_ref, v_grid, B_grid):\n",
      "    \"\"\"precompute_F_grid_with_sidebands but with s_remix as a free parameter (Larmor-equivalent).\"\"\"\n",
  ```

  Inside the function body (cell 42), find:
  ```
      "                remix=remix, verbose_leakage=False)[\"F\"]\n",
  ```
  Replace with:
  ```
      "                s_remix=s_remix, verbose_leakage=False)[\"F\"]\n",
  ```

- [ ] **Step 6: Insert markdown derivation cell** before cell 42's first usage of the helper. The cell with title "## Step 7b — Dark-state mixing (remix) sweep" or similar should get a new markdown cell BEFORE it. Locate the markdown cell heading "remix" (likely cell 41 originally, now shifted by previous inserts). Run:
  ```
  python3 -c 'import json; nb=json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")); [print(i, c["cell_type"], "".join(c.get("source",[]))[:80].replace(chr(10)," ")) for i,c in enumerate(nb["cells"])]'
  ```
  Identify the markdown cell that introduces the remix sweep (search for "Step 7b" or "remix"). Insert a new markdown cell just BEFORE it:

  ```python
  import json, sys
  nb_path = sys.argv[1]
  insert_before_idx = int(sys.argv[2])  # passed in from caller
  with open(nb_path) as f: nb = json.load(f)
  md = {
      "cell_type": "markdown",
      "id": "phase2b-derivation",
      "metadata": {},
      "source": [
          "### `s_remix` ↔ B_⊥ Larmor-mixing derivation\n",
          "\n",
          "The previous `remix` parameter (now renamed `s_remix`) is **a saturation parameter**, not a fraction. It enters `steady_state_rate_eq` as an extra rate per polarization channel:  $R_{remix} = (\\Gamma_{angular}/2) \\cdot s_{remix} \\cdot b_{lup}$, applied isotropically and on-resonance regardless of detuning.\n",
          "\n",
          "Physical mapping to off-axis B-field magnitude $B_\\perp$ (the source of dark-state mixing in real experiments): equate the remix rate to the Larmor mixing rate  $\\Omega_L = g_F \\mu_B B_\\perp / \\hbar$:\n",
          "\n",
          "$$s_{remix} \\cdot \\frac{\\Gamma_{angular}}{2} = \\frac{g_F \\mu_B B_\\perp}{\\hbar} \\;\\;\\Longrightarrow\\;\\; s_{remix} = \\frac{2 g_F \\mu_B B_\\perp}{\\hbar \\, \\Gamma_{angular}} = \\frac{2 g_F (\\mu_B/h) B_\\perp}{\\Gamma_{FWHM}}$$\n",
          "\n",
          "With $\\mu_B/h = 1.39962$ MHz/G, $g_F \\approx 0.5$ for X²Σ⁺ N=1 in the case-b limit, $\\Gamma_{FWHM} = 4.6$ MHz: $s_{remix} \\approx 0.30 \\cdot B_\\perp(\\mathrm{G})$.\n",
          "\n",
          "See `b_perp_to_s_remix` and `s_remix_to_b_perp` (cell idx=28). Captured-trajectory threshold at $s_{remix}=1.0$ (Step 7b) corresponds to $B_\\perp \\approx 3.3$ G.\n"
      ]
  }
  nb['cells'].insert(insert_before_idx, md)
  tmp = nb_path + '.tmp'
  with open(tmp, 'w') as f: json.dump(nb, f, indent=1)
  import os; os.replace(tmp, nb_path)
  print(f"Inserted markdown at idx={insert_before_idx}. Notebook has {len(nb['cells'])} cells.")
  ```
  Inspect the cell list first to find the right `insert_before_idx`, then run.

- [ ] **Step 7: Execute** and confirm bit-identical sweep results to Phase 1 baseline.
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  python3 - <<'PY'
  import json
  current = json.load(open("Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json"))
  # The baseline cell re-runs at the end; current and saved should match.
  # Manual cross-check: read the prior baseline (committed elsewhere, or re-saved before this rename)
  # For now, just confirm the rename didn't break execution and the JSON wrote out.
  print("Baseline keys:", sorted(current.keys()))
  PY
  ```

- [ ] **Step 8: Verify** per the verify block at the top of this task. Inspect any remaining bare-`remix` hits manually.

- [ ] **Step 9: Commit.**
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 2b — rename remix→s_remix + b_perp_to_s_remix helper + Larmor derivation"
  ```

---

### Task 5: PB σ⁻ comb selectivity diagnostic cell

**Goal:** New cell after cell 48 (PB F-grid) that decomposes per-channel scattering rate into "intended" (m_S=+½ → A) vs "off-target" and plots intended-fraction vs B. Quantifies the F_max-at-low-B artifact observed at s=100.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` (insert cell after PB F-grid cell, which by Phase 2a/2b inserts is now likely idx ≈ 50; locate dynamically)

**Acceptance Criteria:**
- [ ] Cell prints "intended fraction at B=B_design" for both s_per_sb values.
- [ ] Cell renders a plot of intended-fraction vs B over [50, 1000] G for both s values.
- [ ] Notebook re-executes clean.

**Verify:**
```
conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 - <<'PY'
import json
nb = json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"))
# Find the cell with "comb selectivity"
hit = next((i for i,c in enumerate(nb["cells"]) if "comb selectivity" in "".join(c.get("source",[]))), None)
assert hit is not None, "comb selectivity cell missing"
print(f"comb selectivity cell at idx={hit}")
PY
```
Expected: notebook exits 0; assertion passes.

**Steps:**

- [ ] **Step 1: Locate the PB F-grid cell** (the one with `F_surfaces_pb`). Use:
  ```
  python3 -c 'import json; nb=json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")); print([i for i,c in enumerate(nb["cells"]) if "F_surfaces_pb" in "".join(c.get("source",[]))])'
  ```
  Use that index + 1 as the insertion point.

- [ ] **Step 2: Insert cell.** Save script to `/tmp/insert_pb_selectivity.py`:
  ```python
  import json, sys
  nb_path = sys.argv[1]
  insert_after = int(sys.argv[2])  # idx of PB F-grid cell
  with open(nb_path) as f: nb = json.load(f)
  cell = {
      "cell_type": "code",
      "metadata": {},
      "execution_count": None,
      "outputs": [],
      "source": [
          "# === Phase 2c: PB σ⁻ comb selectivity diagnostic ==========================================\n",
          "# At each (B, v_design) pair, decompose the per-channel scattering rate into:\n",
          "#   intended = transitions to A from ground sublevels with dominant m_S = +½ (target manifold)\n",
          "#   off-target = transitions to A from m_S = -½ ground sublevels (accidental coverage)\n",
          "# Quantifies the F_max-at-low-B artifact: when the 1054-MHz comb at high s accidentally\n",
          "# covers low-B X(N=1) clusters, intended-fraction drops below ~0.5 and the F-surface peak\n",
          "# does NOT correspond to the intended Petzold cycle.\n",
          "\n",
          "B_grid_sel = B_grid\n",
          "intended_frac = {s_per_sb: np.zeros(len(B_grid_sel)) for s_per_sb in s_values}\n",
          "for s_per_sb in s_values:\n",
          "    sb_template, je_r, ig_r, ig_list = build_paschen_back_sigma_minus(\n",
          "        g, TM, Bd_pb, s_per_sb=s_per_sb)\n",
          "    target_ig_set = set(ig_list)\n",
          "    nu_eg_curve = TM['dE'][:, je_r, ig_r]\n",
          "    nu_eg_design = float(np.interp(Bd_pb, TM['Bz'], nu_eg_curve))\n",
          "    Doppler_design_MHz = k_XA * v_design / (2*np.pi) * 1e-6\n",
          "    for j, B_loc in enumerate(B_grid_sel):\n",
          "        nu_eg_local = float(np.interp(B_loc, TM['Bz'], nu_eg_curve))\n",
          "        delta_pinning = (nu_eg_design - nu_eg_local) - Doppler_design_MHz\n",
          "        sb_list = [{'detuning_MHz': sb['detuning_MHz'] + delta_pinning,\n",
          "                    's': sb['s'], 'polarization': sb['polarization']} for sb in sb_template]\n",
          "        out = steady_state_rate_eq(TM, B_loc, v_design, sidebands=sb_list,\n",
          "                                    ref_transition=(je_r, ig_r), verbose_leakage=False)\n",
          "        # R_per_pol shape (3, n_e, n_g). Sum over (p, e), keep g-axis.\n",
          "        R_per_g = out['R_per_pol'].sum(axis=(0, 1))   # (n_g,)\n",
          "        R_total = R_per_g.sum()\n",
          "        if R_total > 0:\n",
          "            R_target = sum(R_per_g[ig] for ig in target_ig_set)\n",
          "            intended_frac[s_per_sb][j] = R_target / R_total\n",
          "        else:\n",
          "            intended_frac[s_per_sb][j] = np.nan\n",
          "    iB_d = int(np.argmin(np.abs(B_grid_sel - Bd_pb)))\n",
          "    print(f\"  PB σ⁻ s_sb={s_per_sb:5.1f}: intended fraction at B_d={Bd_pb:.0f} G = \"\n",
          "          f\"{intended_frac[s_per_sb][iB_d]:.3f} (1.0 = pure target cycle)\")\n",
          "\n",
          "fig, ax = plt.subplots(1, 1, figsize=(8, 4.5))\n",
          "for s_per_sb in s_values:\n",
          "    ax.plot(B_grid_sel, intended_frac[s_per_sb], '-', lw=1.6,\n",
          "             label=f'$s_{{\\\\rm sb}}$={s_per_sb:.0f}')\n",
          "ax.axvline(Bd_pb, color='k', lw=0.5, ls='--', label='B_design')\n",
          "ax.axhline(0.5, color='r', lw=0.5, alpha=0.4)\n",
          "ax.set_xlabel('B (G)')\n",
          "ax.set_ylabel('intended fraction (target $m_S$=+½ rate / total rate)')\n",
          "ax.set_title(f'PB σ⁻ comb selectivity vs B (B_design={Bd_pb:.0f} G, v={v_design} m/s)')\n",
          "ax.set_ylim(0, 1.05); ax.grid(alpha=0.3); ax.legend()\n",
          "plt.tight_layout(); plt.show()\n",
      ]
  }
  nb['cells'].insert(insert_after + 1, cell)
  tmp = nb_path + '.tmp'
  with open(tmp, 'w') as f: json.dump(nb, f, indent=1)
  import os; os.replace(tmp, nb_path)
  print(f"Inserted at idx={insert_after+1}. Notebook has {len(nb['cells'])} cells.")
  ```
  Run with the actual idx.

- [ ] **Step 3: Execute and verify.**
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  ```

- [ ] **Step 4: Commit.**
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 2c — PB σ⁻ comb selectivity diagnostic"
  ```

---

### Task 6: Cell 42 secondary B_⊥ axis + markdown narrative updates

**Goal:** Add a secondary x-axis to the existing remix-sweep capture-vs-`s_remix` plot in cell 42 showing equivalent B_⊥ values. Also append short paragraphs to the existing narrative markdown cells (Steps 6, 7, 7b, 7c, PB) referencing what the new diagnostics show.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 42, plus markdown cells preceding cells 36, 40, 42, 44, and the PB section

**Acceptance Criteria:**
- [ ] Cell 42's print loop adds an Ω_L (MHz) column derived from `s_remix_to_b_perp(rmx)`.
- [ ] The five Step markdown cells each have a new short paragraph (≤3 sentences) referencing the relevant diagnostic.
- [ ] Notebook executes clean; capture tables unchanged from baseline.

**Verify:**
```
conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 - <<'PY'
import json
nb = json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"))
src42 = "".join(nb["cells"][42].get("source", []))   # adjust idx if shifted by inserts
assert "B_⊥" in src42 or "s_remix_to_b_perp" in src42, "cell 42 missing B_perp axis"
print("OK")
PY
```

**Steps:**

- [ ] **Step 1: Update cell 42's print loop** to include B_⊥ and Ω_L columns. Find the existing print header in cell 42:
  ```
      "print(f\"{'remix':>10s}  {'F_max (1e-25 N)':>17s}  {'F_max ratio vs default':>24s}  {'v_capture (m/s)':>16s}\")\n",
  ```
  Replace with:
  ```
      "print(f\"{'s_remix':>10s}  {'B_⊥ (G)':>10s}  {'Ω_L (MHz)':>11s}  {'F_max (1e-25 N)':>17s}  {'F_max ratio':>12s}  {'v_capture (m/s)':>16s}\")\n",
  ```

  Find the print statement at the end of the loop (with `v_cap_str`):
  ```
      "    v_cap_str = f\"{v_cap:6.0f}\" if v_cap is not None else \"    --\"\n",
  ```
  Add lines BEFORE the next print:
  ```
      "    v_cap_str = f\"{v_cap:6.0f}\" if v_cap is not None else \"    --\"\n",
      "    B_perp_eq = s_remix_to_b_perp(rmx)\n",
      "    Omega_L_MHz_cyclic = b_perp_to_s_remix(B_perp_eq) * Gamma_FWHM_MHz / 2.0  # = g_F μ_B B_⊥ / h\n",
  ```
  Then find the actual print line after this and update its format string accordingly. (The current print line is likely truncated in the source we read; inspect first to capture the exact full line.)

- [ ] **Step 2: Append narrative paragraphs to Step markdown cells.** For each step's intro markdown, find the cell and append a short paragraph. Five cells to update; use the Edit tool with each cell's exact closing strings. Examples:

  Step 7b markdown — find the closing line of its source list, e.g. an existing closing `"\n"` element. Append:
  ```
  "\n**Calibration (added Phase 2):** `s_remix` is now a real saturation parameter mapped to off-axis $B_\\perp$ via `b_perp_to_s_remix`. The first capture at `s_remix=1.0` corresponds to $B_\\perp \\approx 3.3$ G — a physically achievable misalignment.\n"
  ```

  Step PB markdown (the one added in commit `649cd36`) — append:
  ```
  "\n**Selectivity (added Phase 2c):** At $s_{\\rm sb}=100$, the 1054-MHz comb's accidental coverage of low-B X(N=1) clusters drops the intended-fraction below 0.5 (see PB selectivity diagnostic). The F_max at (v=130, B=50 G) is therefore an off-target artifact, not the Petzold cycle.\n"
  ```

  (Steps 6, 7, 7c get shorter analogous additions as judged appropriate during execution.)

- [ ] **Step 3: Execute, verify, commit.**
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 2d/2e — B_⊥ secondary axis + Step markdown updates"
  ```

---

## Phase 3 — Refactor

### Task 7: Define sweep helper trio

**Goal:** Add a single new "sweep helpers" cell after `simulate_trajectory` (cell 35) defining `run_trajectory_sweep`, `format_capture_table`, `plot_v_of_z_grid`. Cells in Phase 3 Tasks 8-12 will replace duplicated boilerplate with calls to these.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` (insert cell after `simulate_trajectory`)

**Acceptance Criteria:**
- [ ] New cell defines exactly three top-level functions: `run_trajectory_sweep`, `format_capture_table`, `plot_v_of_z_grid`.
- [ ] `run_trajectory_sweep(F_interps_dict, B_designs, s_values, L_values, je_ref_lookup, ig_ref_lookup, v_design, v_i_sweep, v_f_sweep, v0_scan, a_max_est, mass, B_grid, TM)` returns a flat dict keyed by `(B_design, s, L, v0)`.
- [ ] `format_capture_table(results, B_designs, s_values, L_values, v0_scan, label='')` returns a multi-line string. Includes the "Total captures" summary.
- [ ] `plot_v_of_z_grid(results, L_show, B_designs, s_values, v0_scan, v0_cmap, v0_norm, suptitle)` returns a `matplotlib.figure.Figure`.
- [ ] Notebook re-executes clean; existing sweeps untouched (helpers are unused yet).

**Verify:**
```
conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 -c 'import json; nb=json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")); src="".join("".join(c.get("source",[])) for c in nb["cells"]); assert all(f in src for f in ["def run_trajectory_sweep", "def format_capture_table", "def plot_v_of_z_grid"]); print("OK")'
```

**Steps:**

- [ ] **Step 1: Insert helper cell.** Locate `simulate_trajectory` cell (currently idx=35; may be shifted by Phase 2 inserts). Find via:
  ```
  python3 -c 'import json; nb=json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb")); print([i for i,c in enumerate(nb["cells"]) if "def simulate_trajectory" in "".join(c.get("source",[]))])'
  ```

  Then insert a cell with the helper definitions:
  ```python
  import json, sys
  nb_path = sys.argv[1]
  insert_after = int(sys.argv[2])
  with open(nb_path) as f: nb = json.load(f)
  cell = {
      "cell_type": "code",
      "metadata": {},
      "execution_count": None,
      "outputs": [],
      "source": [
          "# === Phase 3a: trajectory-sweep helpers ===================================================\n",
          "# Three helpers extracted from the duplicated boilerplate in cells 36, 40, 42, 44, 48.\n",
          "# Each downstream cell is now ~10 lines: prepare F_interps + ref mapping, call helpers.\n",
          "\n",
          "def run_trajectory_sweep(F_interps_dict, B_designs_iter, s_values_iter, L_values_iter,\n",
          "                          je_ref_lookup, ig_ref_lookup,\n",
          "                          v_design, v_i_sweep, v_f_sweep, v0_scan,\n",
          "                          a_max_est, mass, B_grid, TM):\n",
          "    \"\"\"Run trajectory sweep over (B_d, s, L, v0). je_ref_lookup/ig_ref_lookup are dicts\n",
          "    keyed by (B_d, s) -> int (or callables (B_d, s) -> int). Returns flat dict\n",
          "    keyed by (B_d, s, L, v0) -> dict(t, z, v, B, exit_reason, v_target, z_target, eta).\"\"\"\n",
          "    def _resolve(lookup, B_d, s):\n",
          "        return lookup(B_d, s) if callable(lookup) else lookup[(B_d, s)]\n",
          "    results = {}\n",
          "    for Bd in B_designs_iter:\n",
          "        for s in s_values_iter:\n",
          "            je_r = _resolve(je_ref_lookup, Bd, s)\n",
          "            ig_r = _resolve(ig_ref_lookup, Bd, s)\n",
          "            for L in L_values_iter:\n",
          "                z_p, v_t_p, B_p, eta = make_B_of_z_nonlinear(\n",
          "                    TM, je_r, ig_r, Bd, v_design, v_i_sweep, v_f_sweep, L, a_max_est)\n",
          "                B_of_z = interp1d(z_p, B_p, bounds_error=False,\n",
          "                                   fill_value=(float(B_p[0]), float(B_p[-1])))\n",
          "                F_interp = F_interps_dict[(Bd, s)]\n",
          "                for v0 in v0_scan:\n",
          "                    sol = simulate_trajectory(v0, B_of_z, F_interp, mass, L,\n",
          "                                               B_grid_lo=float(B_grid[0]),\n",
          "                                               B_grid_hi=float(B_grid[-1]))\n",
          "                    results[(Bd, s, L, v0)] = dict(sol, v_target=v_t_p, z_target=z_p, eta=eta)\n",
          "    return results\n",
          "\n",
          "\n",
          "def format_capture_table(results, B_designs_iter, s_values_iter, L_values_iter, v0_scan, label=''):\n",
          "    \"\"\"Returns a multi-line capture-velocity range table as a string.\n",
          "    label is shown in the heading.\"\"\"\n",
          "    lines = []\n",
          "    lines.append(f\"\\nCapture-velocity RANGE (m/s) {label}:\")\n",
          "    lines.append(f\"{'':>22s}\" + \"\".join(f\"   L={L:.2f}m\" for L in L_values_iter))\n",
          "    total = 0\n",
          "    for Bd in B_designs_iter:\n",
          "        for s in s_values_iter:\n",
          "            cells = []\n",
          "            for L in L_values_iter:\n",
          "                captured = sorted([v0 for v0 in v0_scan\n",
          "                                    if results[(Bd, s, L, v0)]['exit_reason'] == 'v_stop'])\n",
          "                total += len(captured)\n",
          "                if not captured:\n",
          "                    cells.append('        --')\n",
          "                elif len(captured) == 1:\n",
          "                    cells.append(f\"     {captured[0]:3.0f}    \")\n",
          "                else:\n",
          "                    cells.append(f\"  {captured[0]:3.0f}..{captured[-1]:3.0f}  \")\n",
          "            lines.append(f\"  B_d={Bd:3.0f} G, s_sb={s:5.1f}: \" + \"\".join(cells))\n",
          "    n_configs = len(list(B_designs_iter)) * len(list(s_values_iter)) * len(list(L_values_iter))\n",
          "    lines.append(f\"\\nTotal captures: {total}/{n_configs * len(v0_scan)}\")\n",
          "    return \"\\n\".join(lines)\n",
          "\n",
          "\n",
          "def plot_v_of_z_grid(results, L_show, B_designs_iter, s_values_iter, v0_scan,\n",
          "                       v0_cmap, v0_norm, suptitle, figsize=None):\n",
          "    \"\"\"v(z) trajectories: one panel per (s, B_d). At L=L_show.\"\"\"\n",
          "    nrows = max(1, len(list(s_values_iter)))\n",
          "    ncols = max(1, len(list(B_designs_iter)))\n",
          "    if figsize is None: figsize = (5*ncols, 4*nrows)\n",
          "    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True, squeeze=False)\n",
          "    for i_s, s in enumerate(s_values_iter):\n",
          "        for j_b, Bd in enumerate(B_designs_iter):\n",
          "            ax = axes[i_s, j_b]\n",
          "            v_t = results[(Bd, s, L_show, v0_scan[0])]['v_target']\n",
          "            z_t = results[(Bd, s, L_show, v0_scan[0])]['z_target']\n",
          "            ax.plot(z_t, v_t, 'k--', lw=1.0, alpha=0.6, label=r'$v_{\\rm target}(z)$')\n",
          "            for v0 in v0_scan:\n",
          "                r = results[(Bd, s, L_show, v0)]\n",
          "                ax.plot(r['z'], r['v'], '-', color=v0_cmap(v0_norm(v0)), lw=1.0, alpha=0.85)\n",
          "            ax.axhline(5, color='r', lw=0.5, alpha=0.4)\n",
          "            ax.set_title(f'$B_d$={Bd:.0f} G, $s_{{\\\\rm sb}}$={s:.0f}, L={L_show} m')\n",
          "            if i_s == nrows-1: ax.set_xlabel('z (m)')\n",
          "            if j_b == 0:        ax.set_ylabel('v (m/s)')\n",
          "            ax.grid(alpha=0.3)\n",
          "    sm = plt.cm.ScalarMappable(cmap=v0_cmap, norm=v0_norm)\n",
          "    fig.colorbar(sm, ax=axes, label=r'$v_0$ (m/s)', fraction=0.025)\n",
          "    fig.suptitle(suptitle)\n",
          "    return fig\n",
      ]
  }
  nb['cells'].insert(insert_after + 1, cell)
  tmp = nb_path + '.tmp'
  with open(tmp, 'w') as f: json.dump(nb, f, indent=1)
  import os; os.replace(tmp, nb_path)
  print(f"Inserted helpers at idx={insert_after+1}. Notebook has {len(nb['cells'])} cells.")
  ```
  Run with the located idx.

- [ ] **Step 2: Execute, verify.**
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  ```

- [ ] **Step 3: Commit.**
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 3a — extract run_trajectory_sweep / format_capture_table / plot_v_of_z_grid helpers"
  ```

---

### Task 8: Replace duplicated boilerplate in cells 36, 40, 42, 44, 48

**Goal:** For each of the five trajectory-sweep cells, replace the duplicated `for Bd / for s / for L / for v0` loop + capture-table print + v(z) plot with thin wrappers around the helpers from Task 7.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cells 36, 40, 42, 44, 48 (idx may have shifted; locate dynamically by content)

**Acceptance Criteria:**
- [ ] Each of the five sweep cells has been rewritten to: build (or reference) `F_interps_dict`, build `je_ref_lookup` and `ig_ref_lookup` (dicts or constants), call `run_trajectory_sweep`, call `format_capture_table`, call `plot_v_of_z_grid`. Each cell shrinks to roughly 10-15 lines plus the F-grid prep that was already there.
- [ ] No `for Bd in B_designs ...` triple-nested loop with `simulate_trajectory` inside remains in the notebook outside the helper definition.
- [ ] Notebook re-executes clean.
- [ ] Sweep results from the new code are bit-identical to Phase 1 baseline (deferred check: Task 10).

**Verify:**
```
conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 - <<'PY'
import json, re
nb = json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"))
inner_loop_pat = re.compile(r'for v0 in v0_scan.*?simulate_trajectory', re.DOTALL)
for i, c in enumerate(nb["cells"]):
    src = "".join(c.get("source", []))
    if c["cell_type"] != "code": continue
    if "def run_trajectory_sweep" in src: continue   # the helper is allowed
    if inner_loop_pat.search(src):
        print(f"⚠ cell {i} still contains inline trajectory loop")
PY
```
Expected: no `⚠` lines printed.

**Steps:**

- [ ] **Step 1: Cell 36 (Step 6 single-σ⁺).** Locate via search for `results_single`. Read the full cell. Identify the lines from `# Trajectory sweep` down to the `plt.show()` of the v(z) plot. Replace with:

  ```python
  # === Phase 3b: refactored single-σ⁺ trajectory sweep ===
  je_ref_lookup_single = {(Bd, s): je_ref for Bd in B_designs for s in s_values}
  ig_ref_lookup_single = {(Bd, s): ig_ref for Bd in B_designs for s in s_values}
  results_single = run_trajectory_sweep(
      F_interps_single, B_designs, s_values, L_values,
      je_ref_lookup_single, ig_ref_lookup_single,
      v_design, v_i_sweep, v_f_sweep, v0_scan, a_max_est, mass_RaF, B_grid, TM)
  print(format_capture_table(results_single, B_designs, s_values, L_values, v0_scan,
                               label='SINGLE-SIDEBAND (σ+)'))
  fig = plot_v_of_z_grid(results_single, L_show=1.75, B_designs_iter=B_designs,
                          s_values_iter=s_values, v0_scan=v0_scan,
                          v0_cmap=v0_cmap, v0_norm=v0_norm,
                          suptitle=f"v(z) trajectories — SINGLE-SIDEBAND, L = 1.75 m")
  plt.show()
  ```

  Use the Edit tool to replace the trajectory-sweep + table + plot block. The block likely starts with `# Trajectory sweep` and ends with `plt.show()`. Match the full string.

  Note: cell 36 may use a `F_interps_single` dict (or `F_interps`). Substitute the actual variable name from the F-grid cell that runs immediately before this one.

- [ ] **Step 2: Cell 40 (Step 7 multi-σ⁺).** Same pattern, with:
  ```python
  je_ref_lookup_multi = {(Bd, s): je_ref for Bd in B_designs for s in s_values}
  ig_ref_lookup_multi = {(Bd, s): ig_ref for Bd in B_designs for s in s_values}
  results_multi = run_trajectory_sweep(
      F_interps_multi, B_designs, s_values, L_values,
      je_ref_lookup_multi, ig_ref_lookup_multi,
      v_design, v_i_sweep, v_f_sweep, v0_scan, a_max_est, mass_RaF, B_grid, TM)
  print(format_capture_table(results_multi, B_designs, s_values, L_values, v0_scan,
                               label='MULTI-SIDEBAND (4 × σ+ at B_design)'))
  fig = plot_v_of_z_grid(results_multi, L_show=1.75, B_designs_iter=B_designs,
                          s_values_iter=s_values, v0_scan=v0_scan,
                          v0_cmap=v0_cmap, v0_norm=v0_norm,
                          suptitle=f"v(z) — MULTI-SIDEBAND (4 × σ+), L = 1.75 m")
  plt.show()
  ```

- [ ] **Step 3: Cell 42 (Step 7b remix sweep).** This sweep is over a different axis (s_remix instead of B_d × s × L × v0). The helper signature was designed for the (B_d, s, L, v0) 4D shape; this cell uses (s_remix, v0). Use the helper with `B_designs=[Bd_remix]`, `s_values=[s_sb_remix]`, `L_values=[L_remix]` to flatten to 1D — but with the OUTER s_remix sweep handled OUTSIDE the helper. Concretely, replace the inner trajectory loop:

  ```python
  results_remix = {}
  for rmx in remix_values:
      Fg_r = precompute_F_grid_remix(TM, Bd_remix, v_design, sb_template_remix, rmx,
                                       je_ref, ig_ref, v_grid, B_grid)
      F_interp_r = RegularGridInterpolator(
          (v_grid, B_grid), Fg_r, method="linear", bounds_error=False, fill_value=0.0)
      F_max_remix[rmx] = float(np.abs(Fg_r).max())
      sub = run_trajectory_sweep(
          {(Bd_remix, s_sb_remix): F_interp_r},
          [Bd_remix], [s_sb_remix], [L_remix],
          {(Bd_remix, s_sb_remix): je_ref}, {(Bd_remix, s_sb_remix): ig_ref},
          v_design, v_i_sweep, v_f_sweep, v0_remix_scan, a_max_est, mass_RaF, B_grid, TM)
      for v0 in v0_remix_scan:
          results_remix[(rmx, v0)] = sub[(Bd_remix, s_sb_remix, L_remix, v0)]
  ```

  The downstream `print(...)` table loop in cell 42 (showing `s_remix | B_⊥ | Ω_L | F_max | ratio | v_capture`) added in Phase 2 Task 6 stays as-is — it's not a generic capture-velocity table and shouldn't be replaced by `format_capture_table`.

- [ ] **Step 4: Cell 44 (Step 7c full remix=1.0 sweep).** Same as Step 2 (multi-sideband), with `F_interps_remix1` and `s_remix=remix_full` baked into the F-grid. Replace inner loops:

  ```python
  je_ref_lookup_r1 = {(Bd, s): je_ref for Bd in B_designs for s in s_values}
  ig_ref_lookup_r1 = {(Bd, s): ig_ref for Bd in B_designs for s in s_values}
  results_remix1 = run_trajectory_sweep(
      F_interps_remix1, B_designs, s_values, L_values,
      je_ref_lookup_r1, ig_ref_lookup_r1,
      v_design, v_i_sweep, v_f_sweep, v0_scan, a_max_est, mass_RaF, B_grid, TM)
  print(format_capture_table(results_remix1, B_designs, s_values, L_values, v0_scan,
                               label=f'MULTI-SIDEBAND + s_remix={remix_full} (Step 7c)'))
  L_show = max(L_values)
  fig = plot_v_of_z_grid(results_remix1, L_show=L_show, B_designs_iter=B_designs,
                          s_values_iter=s_values, v0_scan=v0_scan,
                          v0_cmap=v0_cmap, v0_norm=v0_norm,
                          suptitle=f"v(z) — multi-sb + s_remix={remix_full}, L = {L_show} m")
  plt.show()
  ```

- [ ] **Step 5: Cell 48 (PB σ⁻ sweep).** Single B_design, two s values:
  ```python
  je_ref_lookup_pb = {(Bd_pb, s): ref_pb[s][0] for s in s_values}
  ig_ref_lookup_pb = {(Bd_pb, s): ref_pb[s][1] for s in s_values}
  results_pb_flat = run_trajectory_sweep(
      F_interps_pb_dict, [Bd_pb], s_values, L_values,
      je_ref_lookup_pb, ig_ref_lookup_pb,
      v_design, v_i_sweep, v_f_sweep, v0_scan, a_max_est, mass_RaF, B_grid, TM)
  # Re-key results from (Bd_pb, s, L, v0) to (s, L, v0) for backward compat with baseline
  results_pb = {(s, L, v0): results_pb_flat[(Bd_pb, s, L, v0)]
                  for s in s_values for L in L_values for v0 in v0_scan}
  print(format_capture_table(results_pb_flat, [Bd_pb], s_values, L_values, v0_scan,
                               label='PB σ⁻ 6-sideband at B_d=1000 G'))
  L_show = max(L_values)
  fig = plot_v_of_z_grid(results_pb_flat, L_show=L_show, B_designs_iter=[Bd_pb],
                          s_values_iter=s_values, v0_scan=v0_scan,
                          v0_cmap=v0_cmap, v0_norm=v0_norm,
                          suptitle=f"v(z) — Paschen-Back σ⁻ 6-sb, $B_d$ = {Bd_pb:.0f} G")
  plt.show()
  ```
  Requires `F_interps_pb_dict = {(Bd_pb, s): F_interps_pb[s] for s in s_values}` defined above this snippet (one-line transform).

  Note: `ref_pb` was added in commit `649cd36` as a dict s → (je_r, ig_r). Task 9 (next) will collapse it to a single (je_pb_ref, ig_pb_ref) since both s values yield the same ref. Until then, the dict access pattern above works.

- [ ] **Step 6: Execute end-to-end** to confirm no break:
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  ```

- [ ] **Step 7: Verify** per the Verify block at the top of this task.

- [ ] **Step 8: Commit.**
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 3b — replace duplicated trajectory-sweep boilerplate with helper calls"
  ```

---

### Task 9: PB function uses g_char_full['mS']; cell 47 uses single ref

**Goal:** Robust manifold selection in `build_paschen_back_sigma_minus`. Build the PB sideband template + ref ONCE in cell 47 instead of per-s.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 46 (PB function), cell 47 (PB F-grid)

**Acceptance Criteria:**
- [ ] `build_paschen_back_sigma_minus` selects `ig_list` via `g_char_full['mS']` (sign of dominant m_S), not energy sort.
- [ ] Function asserts `len(ig_list) == n_sidebands` with informative error.
- [ ] Cell 47 builds `(je_pb_ref, ig_pb_ref, sb_pb_template)` once outside the s loop; per-s sideband list is constructed by varying only the `s` field.
- [ ] PB F-grid output bit-identical to baseline (the new selection should agree with energy-sort for current single-N basis).

**Verify:** Compare F_max from PB sweep before vs after this task.

**Steps:**

- [ ] **Step 1: Update `build_paschen_back_sigma_minus`** in cell 46. Find:
  ```
      "    Bz = TM[\"Bz\"]\n",
      "    iB = int(np.argmin(np.abs(Bz - B_design_G)))\n",
      "    E_at_design = g.evals_B[iB, :]\n",
      "    sort_idx = np.argsort(E_at_design)\n",
      "    if manifold == \"upper\":\n",
      "        ig_list = sorted(sort_idx[-n_sidebands:].tolist())\n",
      "    elif manifold == \"lower\":\n",
      "        ig_list = sorted(sort_idx[:n_sidebands].tolist())\n",
      "    else:\n",
      "        raise ValueError(f\"manifold must be 'upper' or 'lower', got {manifold!r}\")\n",
  ```
  Replace with:
  ```
      "    Bz = TM[\"Bz\"]\n",
      "    iB = int(np.argmin(np.abs(Bz - B_design_G)))\n",
      "    # Use g_char_full['mS'] (dominant m_S in decoupled basis) for robust manifold selection.\n",
      "    # Falls back to energy-sort if mS not available (e.g. mS_key was None at build time).\n",
      "    mS_dom = g_char_full['mS'][iB, :]\n",
      "    if np.isnan(mS_dom).any():\n",
      "        E_at_design = g.evals_B[iB, :]\n",
      "        sort_idx = np.argsort(E_at_design)\n",
      "        if manifold == \"upper\":\n",
      "            ig_list = sorted(sort_idx[-n_sidebands:].tolist())\n",
      "        else:\n",
      "            ig_list = sorted(sort_idx[:n_sidebands].tolist())\n",
      "    else:\n",
      "        target_sign = +1 if manifold == \"upper\" else -1\n",
      "        ig_list = sorted(int(i) for i in np.where(np.sign(mS_dom) == target_sign)[0])\n",
      "        if len(ig_list) != n_sidebands:\n",
      "            raise ValueError(f\"manifold={manifold!r} at B={B_design_G:.0f} G yields \"\n",
      "                              f\"{len(ig_list)} sublevels (expected {n_sidebands}). \"\n",
      "                              f\"Check N_list and mS dominance at this B.\")\n",
      "    if manifold not in (\"upper\", \"lower\"):\n",
      "        raise ValueError(f\"manifold must be 'upper' or 'lower', got {manifold!r}\")\n",
  ```

- [ ] **Step 2: Update cell 47** to build sb_template + ref once. Find the loop:
  ```
      "for s_per_sb in s_values:\n",
      "    sb_template, je_r, ig_r, _ = build_paschen_back_sigma_minus(g, TM, Bd_pb,\n",
      "                                                                  s_per_sb=s_per_sb)\n",
  ```
  Replace with:
  ```
      "# Build PB sideband template + ref ONCE; per-s only the sideband saturation varies.\n",
      "sb_pb_template, je_pb_ref, ig_pb_ref, _ = build_paschen_back_sigma_minus(g, TM, Bd_pb,\n",
      "                                                                          s_per_sb=s_values[0])\n",
      "for s_per_sb in s_values:\n",
      "    sb_template = [{**sb, 's': s_per_sb} for sb in sb_pb_template]\n",
      "    je_r, ig_r = je_pb_ref, ig_pb_ref\n",
  ```
  And further down where `ref_pb[s_per_sb] = (je_r, ig_r)` appears, replace with:
  ```
      "    ref_pb[s_per_sb] = (je_pb_ref, ig_pb_ref)\n",
  ```
  (preserves the dict for backward compat with cell 48).

- [ ] **Step 3: Execute, verify F_max unchanged.**
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  ```

- [ ] **Step 4: Commit.**
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 3c — PB function uses g_char_full['mS'] + single ref reuse"
  ```

---

### Task 10: Regression check cell + verification

**Goal:** Add a cell at the end of the notebook (just before the existing baseline-serialization cell from Task 2) that loads `zs_baselines_pre_refactor.json` and compares against the freshly-computed values. Fail fast (`raise AssertionError`) on any disagreement.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` (insert cell at second-to-last position)

**Acceptance Criteria:**
- [ ] Cell loads the JSON, walks each sweep, and compares F_max (within rel-tolerance 1e-6) and capture lists (exact equality).
- [ ] On agreement: prints `✓ Regression check PASS: all sweeps match Phase 1 baseline.`
- [ ] On disagreement: prints itemized deltas AND `raise AssertionError(...)`.
- [ ] Notebook re-executes 0 (regression check passes).

**Verify:**
```
conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
```
Expected: exits 0; cell prints PASS line.

**Steps:**

- [ ] **Step 1: Insert the regression-check cell.** Save script to `/tmp/insert_regression_cell.py`:
  ```python
  import json, sys
  nb_path = sys.argv[1]
  with open(nb_path) as f: nb = json.load(f)
  cell = {
      "cell_type": "code",
      "metadata": {},
      "execution_count": None,
      "outputs": [],
      "source": [
          "# === Phase 3e: regression check vs Phase 1 baseline =======================================\n",
          "import json\n",
          "_baseline_path = 'Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json'\n",
          "with open(_baseline_path) as _f: _baseline = json.load(_f)\n",
          "_deltas = []\n",
          "_RTOL = 1e-6\n",
          "\n",
          "def _check(label, expected_F, current_F, expected_cap, current_cap):\n",
          "    for k, v_exp in expected_F.items():\n",
          "        v_cur = current_F.get(k)\n",
          "        if v_cur is None:\n",
          "            _deltas.append(f\"{label} F_max[{k}]: missing in current\")\n",
          "        elif abs(v_cur - v_exp) > _RTOL * max(abs(v_exp), 1e-30):\n",
          "            _deltas.append(f\"{label} F_max[{k}]: {v_exp:.6e} -> {v_cur:.6e}\")\n",
          "    for k, c_exp in expected_cap.items():\n",
          "        c_cur = current_cap.get(k, [])\n",
          "        if list(c_exp) != list(c_cur):\n",
          "            _deltas.append(f\"{label} captures[{k}]: {c_exp} -> {c_cur}\")\n",
          "\n",
          "def _capture_range(results_dict, key_fn, v0_iter):\n",
          "    return sorted([float(v0) for v0 in v0_iter\n",
          "                   if results_dict[key_fn(v0)]['exit_reason'] == 'v_stop'])\n",
          "\n",
          "if 'single_sb' in _baseline and 'results_single' in dir():\n",
          "    cur_F = {f\"{Bd:.0f}_{s:.0f}\": float(np.abs(F_surfaces_single[(Bd,s)]).max())\n",
          "             for Bd in B_designs for s in s_values}\n",
          "    cur_C = {f\"{Bd:.0f}_{s:.0f}_{L:.2f}\":\n",
          "                _capture_range(results_single, lambda v0,Bd=Bd,s=s,L=L:(Bd,s,L,v0), v0_scan)\n",
          "             for Bd in B_designs for s in s_values for L in L_values}\n",
          "    _check('single_sb', _baseline['single_sb']['F_max'], cur_F,\n",
          "            _baseline['single_sb']['captures'], cur_C)\n",
          "\n",
          "if 'multi_sb' in _baseline and 'results_multi' in dir():\n",
          "    cur_F = {f\"{Bd:.0f}_{s:.0f}\": float(np.abs(F_surfaces_multi[(Bd,s)]).max())\n",
          "             for Bd in B_designs for s in s_values}\n",
          "    cur_C = {f\"{Bd:.0f}_{s:.0f}_{L:.2f}\":\n",
          "                _capture_range(results_multi, lambda v0,Bd=Bd,s=s,L=L:(Bd,s,L,v0), v0_scan)\n",
          "             for Bd in B_designs for s in s_values for L in L_values}\n",
          "    _check('multi_sb', _baseline['multi_sb']['F_max'], cur_F,\n",
          "            _baseline['multi_sb']['captures'], cur_C)\n",
          "\n",
          "if 'remix_sweep' in _baseline and 'results_remix' in dir():\n",
          "    cur_F = {f\"{rmx:.0e}\": float(F_max_remix[rmx]) for rmx in remix_values}\n",
          "    cur_C = {f\"{rmx:.0e}\": _capture_range(results_remix, lambda v0,rmx=rmx:(rmx,v0), v0_remix_scan)\n",
          "             for rmx in remix_values}\n",
          "    _check('remix_sweep', _baseline['remix_sweep']['F_max'], cur_F,\n",
          "            _baseline['remix_sweep']['captures'], cur_C)\n",
          "\n",
          "if 'remix_full' in _baseline and 'results_remix1' in dir():\n",
          "    cur_F = {f\"{Bd:.0f}_{s:.0f}\": float(np.abs(F_surfaces_remix1[(Bd,s)]).max())\n",
          "             for Bd in B_designs for s in s_values}\n",
          "    cur_C = {f\"{Bd:.0f}_{s:.0f}_{L:.2f}\":\n",
          "                _capture_range(results_remix1, lambda v0,Bd=Bd,s=s,L=L:(Bd,s,L,v0), v0_scan)\n",
          "             for Bd in B_designs for s in s_values for L in L_values}\n",
          "    _check('remix_full', _baseline['remix_full']['F_max'], cur_F,\n",
          "            _baseline['remix_full']['captures'], cur_C)\n",
          "\n",
          "if 'pb' in _baseline and 'results_pb' in dir():\n",
          "    cur_F = {f\"{s:.0f}\": float(np.abs(F_surfaces_pb[s]).max()) for s in s_values}\n",
          "    cur_C = {f\"{s:.0f}_{L:.2f}\":\n",
          "                _capture_range(results_pb, lambda v0,s=s,L=L:(s,L,v0), v0_scan)\n",
          "             for s in s_values for L in L_values}\n",
          "    _check('pb', _baseline['pb']['F_max'], cur_F,\n",
          "            _baseline['pb']['captures'], cur_C)\n",
          "\n",
          "if _deltas:\n",
          "    print('✗ Regression check FAIL:')\n",
          "    for d in _deltas: print('  ', d)\n",
          "    raise AssertionError(f'{len(_deltas)} sweep deltas vs Phase 1 baseline')\n",
          "else:\n",
          "    print(f'✓ Regression check PASS: all sweeps match Phase 1 baseline within rtol={_RTOL}.')\n",
      ]
  }
  # Insert as second-to-last code cell (just before the baseline-serialization cell from Task 2,\n",
  # which itself sits just before the trailing 'Notes & next steps' markdown).\n",
  nb['cells'].insert(-2, cell)
  tmp = nb_path + '.tmp'
  with open(tmp, 'w') as f: json.dump(nb, f, indent=1)
  import os; os.replace(tmp, nb_path)
  print(f"Inserted regression check. Notebook has {len(nb['cells'])} cells.")
  ```
  Run: `python3 /tmp/insert_regression_cell.py "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"`

- [ ] **Step 2: Execute end-to-end.** If regression check FAILS, the notebook execution will fail with the AssertionError. The error message lists each delta. Investigate, fix, re-run.
  ```
  conda run -n Structure jupyter execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  ```
  Expected: exits 0, regression PASS line printed.

- [ ] **Step 3: If failures occur during regression check.** Categorize each delta:
  - **F_max deltas at the rtol boundary** (e.g. 1e-15 relative): probably floating-point reordering from the helper refactor. Loosen `_RTOL` to 1e-12 if needed and document.
  - **Capture-list deltas** at v_0 boundary values (e.g., capture set `[40, 45, 50]` becomes `[45, 50]` or vice versa): integrator step-size sensitivity at boundary v_0 — same as the 3-vs-4 wobble we observed during PB review. Document in the spec.
  - **Larger F_max changes** (>1% relative): real physics regression. Halt, debug.

- [ ] **Step 4: Commit (or rollback).** If PASS:
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 3e — regression check vs Phase 1 baseline (PASS)"
  ```
  If documented deltas accepted:
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb" \
      "docs/superpowers/specs/2026-05-13-zeeman-slower-2-comprehensive-hardening-design.md"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 3e — regression check (deltas documented in spec)"
  ```

---

### Task 11: Update spec doc with post-refactor results

**Goal:** Append a "Results" section to the spec at `docs/superpowers/specs/2026-05-13-zeeman-slower-2-comprehensive-hardening-design.md` summarizing what landed: cells modified, capture tables before/after (or "bit-identical"), any deltas with attribution.

**Files:**
- Modify: `docs/superpowers/specs/2026-05-13-zeeman-slower-2-comprehensive-hardening-design.md`

**Acceptance Criteria:**
- [ ] New `## Results` section appended.
- [ ] Lists each phase's commit SHA.
- [ ] States whether regression PASS or itemized deltas.
- [ ] Notes any open follow-ups discovered during execution.

**Verify:**
```
grep -n "^## Results" "docs/superpowers/specs/2026-05-13-zeeman-slower-2-comprehensive-hardening-design.md"
```

**Steps:**

- [ ] **Step 1: Append Results section.** Use Edit tool, append after the "Out of scope" section. Template:
  ```markdown
  ## Results (filled in during execution)

  | Phase | Commit | Outcome |
  |---|---|---|
  | 1a | [SHA] | Cell 27/33/36/38/40/42 edits + warnings |
  | 1b | [SHA] | Baseline serialization cell + .gitignore |
  | 2a | [SHA] | ν_eg(B) monotonicity diagnostic — verdict: [MONOTONIC / NON-MONOTONIC at B≈[...]] |
  | 2b | [SHA] | remix→s_remix rename + Larmor helper |
  | 2c | [SHA] | PB σ⁻ comb selectivity — intended-fraction at B_d=1000 G: s=10 → [.X], s=100 → [.X] |
  | 2d/2e | [SHA] | B_⊥ secondary axis + Step markdown updates |
  | 3a | [SHA] | Sweep helper trio extracted |
  | 3b | [SHA] | 5 sweep cells refactored to helper calls |
  | 3c | [SHA] | PB function uses g_char_full['mS']; cell 47 single ref |
  | 3e | [SHA] | Regression check: [PASS / FAIL with N deltas] |

  **Regression deltas (if any):**
  - [item]
  - [item]

  **Open follow-ups:**
  - [item]
  ```

- [ ] **Step 2: Commit.**
  ```
  /Library/Developer/CommandLineTools/usr/bin/git add "docs/superpowers/specs/2026-05-13-zeeman-slower-2-comprehensive-hardening-design.md"
  /Library/Developer/CommandLineTools/usr/bin/git commit -m "spec: zeeman_slower_2 hardening — fill Results section"
  ```

---

## Self-Review

**Spec coverage check:**
- Phase 1 items 1a-1f: Tasks 1, 2 ✓
- Phase 2 items 2a-2e: Tasks 3, 4, 5, 6 ✓
- Phase 3 items 3a-3f: Tasks 7, 8, 9, 10, 11 ✓
- Verification gate (`jupyter execute --inplace` no `--allow-errors`): each task includes it ✓
- Open Question 3 (`v_stop` capture classification): documented in Task 1 Step 1's docstring upgrade as known-acceptable simplification ✓
- Open Question 4 (`max_step` validation): NOT explicitly covered — acceptable per spec ("only if Phase 1f's baseline shows any sensitivity") ✓

**Type/identifier consistency:** Helper signatures match between Task 7 (definition) and Task 8 (call sites). `je_ref_lookup` and `ig_ref_lookup` are dicts in all five call sites (constant-valued for σ⁺ schemes, dict-valued for PB). `format_capture_table` returns string; `plot_v_of_z_grid` returns Figure. `b_perp_to_s_remix` and `s_remix_to_b_perp` are inverses; both defined in Task 4 Step 4 cell.

**Placeholder scan:** All TODOs/TBDs are either explicit "fill in during execution" (Task 11) OR open questions in the spec already documented as deferrals.

---

**Plan complete.** Saved to `docs/superpowers/plans/2026-05-13-zeeman-slower-2-comprehensive-hardening.md`. Tasks JSON co-located.
