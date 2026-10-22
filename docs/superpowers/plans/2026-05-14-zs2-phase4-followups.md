# zeeman_slower_2 Phase 4 follow-ups — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers-extended-cc:subagent-driven-development` (recommended) or `superpowers-extended-cc:executing-plans` to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply the three Phase 4 fixes from [the spec](../specs/2026-05-14-zs2-phase4-followups-design.md): (Fix 1) extract `_captured_at` helper; (Fix 3) two-phase split + lazy-init removal in cell 46; (Fix 2) env-var bless + git-tracked baseline.

**Architecture:** Sequential commits in execution order Fix 1 → Fix 3 → Fix 2 (Fix 3 consumes Fix 1's helper; Fix 2 is independent). Notebook-only (no `Source Code/` changes). Each commit must produce a clean `jupyter nbconvert --to notebook --execute --inplace` run with zero `output_type=='error'` cells, validated by reading the .ipynb JSON. **Never use `--allow-errors`** (cf. `feedback_no_allow_errors.md`).

**Tech Stack:** Python 3 (conda env `Structure`), numpy, scipy, matplotlib, jupyter. macOS git path: `/Library/Developer/CommandLineTools/usr/bin/git` (xcrun shim hangs the standard path; cf. `~/.claude/lessons.md` 2026-05-09).

---

## File map

| File | Touched by |
|---|---|
| `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 38 | Task 1 |
| `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 46 | Task 2 |
| `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cells 54, 55 | Task 3 |
| `Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json` (renamed) | Task 3 |
| `.gitignore` (line 16) | Task 3 |

**Cell-number drift guard:** indices above are valid at HEAD `3ea7493`. If they have shifted, locate by content:
- Helper-trio cell: `# === Phase 3a: trajectory-sweep helpers`
- s_remix sweep cell: `def precompute_F_grid_remix`
- Regression-check cell: `# === Phase 3e: regression check vs prior-run baseline`
- Baseline-serializer cell: `# === Phase 1f: serialize sweep results as regression baseline`

---

### Task 1: Extract `_captured_at` helper; refactor `format_capture_table` to use it

**Goal:** Add a fourth helper to cell 38 (`_captured_at`) and refactor `format_capture_table` to call it. Establishes the shared primitive that Fix 3 (cell 46) will also use.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 38

**Acceptance Criteria:**
- [ ] Cell 38 defines `_captured_at(results, key_fn, v0_iter)` returning `sorted([float(v0) for v0 in v0_iter if results[key_fn(v0)]["exit_reason"] == "v_stop"])`.
- [ ] `format_capture_table`'s inline `sorted([v0 for v0 in v0_scan if ...])` is replaced by a call to `_captured_at`.
- [ ] No other cells touched in this task.
- [ ] Notebook re-executes clean (zero error cells; cell 54 regression PASS).

**Verify:**
```bash
conda run -n Structure jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=900 \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 - <<'PY'
import json
nb = json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"))
errs = sum(1 for c in nb["cells"] if c["cell_type"]=="code"
            for o in c.get("outputs", []) if o.get("output_type")=="error")
src = "".join("".join(c.get("source", [])) for c in nb["cells"])
assert errs == 0, f"{errs} error cells"
assert "def _captured_at" in src, "_captured_at missing"
# Confirm format_capture_table calls it
assert "_captured_at(results, lambda v0" in src, "format_capture_table not refactored"
# Confirm cell 54 regression PASS line is in outputs
saw_pass = False
for c in nb["cells"]:
    for o in c.get("outputs", []):
        if o.get("output_type") == "stream":
            if "Regression check PASS" in "".join(o.get("text", [])):
                saw_pass = True
assert saw_pass, "regression PASS line not found in outputs"
print("OK")
PY
```
Expected: `OK`.

**Steps:**

- [ ] **Step 1: Edit cell 38 source via Python script.** Save to `/tmp/zs_p4_task1.py`:

```python
"""Phase 4 Fix 1 — add _captured_at helper, refactor format_capture_table."""
import json, ast, os

NB = "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"

CELL_38 = """\
# === Phase 3a: trajectory-sweep helpers ===================================================
# Three helpers extracted from the duplicated boilerplate in cells 39, 43, 48, 53.
# (Cell 46 keeps its inline loop because its sweep axis is s_remix and its print/plot
# layout is bespoke - wrapping the helper around it would add scaffolding, not remove it.
# Phase 4 Fix 1 added _captured_at as a shared primitive used by both
# format_capture_table and cell 46's bespoke print table.)
# Each refactored downstream cell is now ~10 lines: prepare F_interps + ref mapping,
# call helpers.

def _captured_at(results, key_fn, v0_iter):
    \"\"\"Sorted list of v0 (float) captured (exit_reason=='v_stop') for a given key_fn.
    Used by format_capture_table and by cell 46's bespoke print table.\"\"\"
    return sorted([float(v0) for v0 in v0_iter
                   if results[key_fn(v0)]["exit_reason"] == "v_stop"])


def run_trajectory_sweep(F_interps_dict, B_designs_iter, s_values_iter, L_values_iter,
                          je_ref_lookup, ig_ref_lookup,
                          v_design, v_i_sweep, v_f_sweep, v0_scan,
                          a_max_est, mass, B_grid, TM):
    \"\"\"Run trajectory sweep over (B_d, s, L, v0).
    je_ref_lookup / ig_ref_lookup are dicts keyed by (B_d, s) -> int.
    Returns flat dict keyed by (B_d, s, L, v0) -> dict(t, z, v, B, exit_reason,
    v_target, z_target, eta).\"\"\"
    results = {}
    for Bd in B_designs_iter:
        for s in s_values_iter:
            je_r = je_ref_lookup[(Bd, s)]
            ig_r = ig_ref_lookup[(Bd, s)]
            for L in L_values_iter:
                z_p, v_t_p, B_p, eta = make_B_of_z_nonlinear(
                    TM, je_r, ig_r, Bd, v_design, v_i_sweep, v_f_sweep, L, a_max_est)
                B_of_z = interp1d(z_p, B_p, bounds_error=False,
                                   fill_value=(float(B_p[0]), float(B_p[-1])))
                F_interp = F_interps_dict[(Bd, s)]
                for v0 in v0_scan:
                    sol = simulate_trajectory(v0, B_of_z, F_interp, mass, L,
                                               B_grid_lo=float(B_grid[0]),
                                               B_grid_hi=float(B_grid[-1]))
                    results[(Bd, s, L, v0)] = dict(sol, v_target=v_t_p, z_target=z_p, eta=eta)
    return results


def format_capture_table(results, B_designs_iter, s_values_iter, L_values_iter, v0_scan, label=""):
    \"\"\"Multi-line capture-velocity RANGE table as a string. Captured = exit_reason 'v_stop'.\"\"\"
    Bds = list(B_designs_iter); ss = list(s_values_iter); Ls = list(L_values_iter)
    lines = []
    lines.append(f"\\nCapture-velocity RANGE (m/s) {label}:")
    lines.append(f"{'':>22s}" + "".join(f"   L={L:.2f}m" for L in Ls))
    total = 0
    for Bd in Bds:
        for s in ss:
            cells = []
            for L in Ls:
                captured = _captured_at(results,
                                          lambda v0, Bd=Bd, s=s, L=L: (Bd, s, L, v0),
                                          v0_scan)
                total += len(captured)
                if not captured:
                    cells.append("        --")
                elif len(captured) == 1:
                    cells.append(f"     {captured[0]:3.0f}    ")
                else:
                    cells.append(f"  {captured[0]:3.0f}..{captured[-1]:3.0f}  ")
            lines.append(f"  B_d={Bd:3.0f} G, s_sb={s:5.1f}: " + "".join(cells))
    n_configs = len(Bds) * len(ss) * len(Ls)
    lines.append(f"\\nTotal captures: {total}/{n_configs * len(v0_scan)}")
    return "\\n".join(lines)


def plot_v_of_z_grid(results, L_show, B_designs_iter, s_values_iter, v0_scan,
                       v0_cmap, v0_norm, suptitle, figsize=None):
    \"\"\"v(z) trajectories: one panel per (s, B_d) at L=L_show. Returns Figure.\"\"\"
    Bds = list(B_designs_iter); ss = list(s_values_iter)
    nrows, ncols = max(1, len(ss)), max(1, len(Bds))
    if figsize is None: figsize = (5*ncols, 4*nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                              sharex=True, sharey=True, squeeze=False)
    for i_s, s in enumerate(ss):
        for j_b, Bd in enumerate(Bds):
            ax = axes[i_s, j_b]
            v_t = results[(Bd, s, L_show, v0_scan[0])]["v_target"]
            z_t = results[(Bd, s, L_show, v0_scan[0])]["z_target"]
            ax.plot(z_t, v_t, "k--", lw=1.0, alpha=0.6, label=r"$v_{\\rm target}(z)$")
            for v0 in v0_scan:
                r = results[(Bd, s, L_show, v0)]
                ax.plot(r["z"], r["v"], "-", color=v0_cmap(v0_norm(v0)), lw=1.0, alpha=0.85)
            ax.axhline(5, color="r", lw=0.5, alpha=0.4)
            ax.set_title(f"$B_d$={Bd:.0f} G, $s_{{\\\\rm sb}}$={s:.0f}, L={L_show} m")
            if i_s == nrows-1: ax.set_xlabel("z (m)")
            if j_b == 0:        ax.set_ylabel("v (m/s)")
            ax.grid(alpha=0.3)
    sm = plt.cm.ScalarMappable(cmap=v0_cmap, norm=v0_norm)
    fig.colorbar(sm, ax=axes, label=r"$v_0$ (m/s)", fraction=0.025)
    fig.suptitle(suptitle)
    return fig
"""

ast.parse(CELL_38)
print(f"OK syntax cell 38: {len(CELL_38.splitlines())} lines")

# Locate cell containing "Phase 3a: trajectory-sweep helpers" header
nb = json.load(open(NB))
target = None
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] == "code" and "Phase 3a: trajectory-sweep helpers" in "".join(c.get("source", [])):
        target = i; break
assert target is not None, "helper-trio cell not found"
print(f"Targeting cell {target}")

nb["cells"][target]["source"] = CELL_38.splitlines(keepends=True)
nb["cells"][target]["outputs"] = []
nb["cells"][target]["execution_count"] = None

tmp = NB + ".tmp"
with open(tmp, "w") as f: json.dump(nb, f, indent=1)
os.replace(tmp, NB)
print(f"Cell {target} replaced. Total cells: {len(nb['cells'])}")
```

Run: `python3 /tmp/zs_p4_task1.py`
Expected: prints `OK syntax cell 38: ~95 lines`, `Targeting cell 38`, `Cell 38 replaced. Total cells: 57`.

- [ ] **Step 2: Re-execute notebook end-to-end.**

```bash
conda run -n Structure jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=900 \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
```
Expected: prints `[NbConvertApp] Writing ... bytes to ...`. No tracebacks.

- [ ] **Step 3: Verify per the Verify block at the top of this task.** Run the inline Python check. Expected: `OK`.

- [ ] **Step 4: Commit.**

If you get `fatal: Unable to create '.git/index.lock'` (Google Drive sync race), wait 3 seconds and retry — the lock disappears on its own.

```bash
/Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
/Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 4 Fix 1 — extract _captured_at helper, refactor format_capture_table"
```

---

### Task 2: Cell 46 two-phase split using `_captured_at`

**Goal:** Restructure cell 46's main loop into two phases (precompute all F-grids, then trajectory sweep + print). Removes the `F_max_default = None` lazy-init. Composes with Task 1's `_captured_at`.

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cell 46

**Acceptance Criteria:**
- [ ] Cell 46 contains exactly one `for rmx in remix_values:` loop in each of two phases (precompute + trajectory). No `F_max_default = None` line.
- [ ] `F_max_default = F_max_remix[remix_values[0]]` appears explicitly between the two phases.
- [ ] The captured-v0 list is built via `_captured_at(results_remix, lambda v0, rmx=rmx: (rmx, v0), v0_remix_scan)`.
- [ ] Cell 46's print-table output is bit-identical to pre-refactor (same s_remix / B_⊥ / Ω_L / F_max / ratio / v_capture rows).
- [ ] The v(z)/v(t) overlay plot below the table is unchanged.
- [ ] Notebook re-executes clean; cell 54 regression PASS.

**Verify:**
```bash
conda run -n Structure jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=900 \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 - <<'PY'
import json, re
nb = json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"))
# Find cell 46 by content
target = None
for i, c in enumerate(nb["cells"]):
    if c["cell_type"]=="code" and "def precompute_F_grid_remix" in "".join(c.get("source", [])):
        target = i; break
src = "".join(nb["cells"][target]["source"])
assert "F_max_default = None" not in src, "lazy-init still present"
assert "F_max_default = F_max_remix[remix_values[0]]" in src, "explicit anchor missing"
assert "F_grids_remix = {}" in src or "F_grids_remix =" in src, "two-phase split missing"
assert "_captured_at(results_remix" in src, "_captured_at not used in cell 46"
errs = sum(1 for c in nb["cells"] if c["cell_type"]=="code"
            for o in c.get("outputs", []) if o.get("output_type")=="error")
assert errs == 0, f"{errs} error cells"
saw_pass = any("Regression check PASS" in "".join(o.get("text", []))
                for c in nb["cells"] for o in c.get("outputs", [])
                if o.get("output_type") == "stream")
assert saw_pass, "regression PASS line not found"
print("OK")
PY
```
Expected: `OK`.

**Steps:**

- [ ] **Step 1: Edit cell 46 via Python script.** Save to `/tmp/zs_p4_task2.py`:

```python
"""Phase 4 Fix 3 — cell 46 two-phase split using _captured_at."""
import json, ast, os

NB = "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"

CELL_46 = """\
def precompute_F_grid_remix(TM, B_design_G, v_design, sidebands_template, s_remix,
                              je_ref, ig_ref, v_grid, B_grid):
    \"\"\"precompute_F_grid_with_sidebands but with s_remix as a free parameter (Larmor-equivalent).\"\"\"
    Bz = TM["Bz"]
    nu_eg_curve = TM["dE"][:, je_ref, ig_ref]
    nu_eg_design = float(np.interp(B_design_G, Bz, nu_eg_curve))
    Doppler_design_MHz = k_XA * v_design / (2*np.pi) * 1e-6
    F_grid = np.zeros((len(v_grid), len(B_grid)))
    for j, B_loc in enumerate(B_grid):
        nu_eg_local = float(np.interp(B_loc, Bz, nu_eg_curve))
        delta_pinning = (nu_eg_design - nu_eg_local) - Doppler_design_MHz
        sb_list = [{"detuning_MHz": sb["detuning_MHz"] + delta_pinning,
                    "s": sb["s"], "polarization": sb["polarization"]}
                    for sb in sidebands_template]
        for i, v in enumerate(v_grid):
            F_grid[i, j] = steady_state_rate_eq(
                TM, B_loc, v, sidebands=sb_list,
                ref_transition=(je_ref, ig_ref),
                s_remix=s_remix, verbose_leakage=False)["F"]
    return F_grid


# Most-aggressive config from Step 7
Bd_remix    = 150.
s_sb_remix  = 100.
L_remix     = 2.75
remix_values = [1e-3, 1e-2, 1e-1, 1.0, 10.0]
v0_remix_scan = np.arange(40., 121., 5.)

# Build sidebands once (don't depend on remix)
sb_template_remix, _ = build_X_N1_sidebands(g, TM, je_ref, ig_ref, Bd_remix,
                                              "sigma+", s_sb_remix)

# Build B(z) profile once
z_p_r, v_t_r, B_p_r, eta_r = make_B_of_z_nonlinear(
    TM, je_ref, ig_ref, Bd_remix, v_design, v_i_sweep, v_f_sweep, L_remix, a_max_est)
B_of_z_remix = interp1d(z_p_r, B_p_r, bounds_error=False,
                          fill_value=(float(B_p_r[0]), float(B_p_r[-1])))

# === Phase 4 Fix 3: two-phase split ============================================
# Phase 1 — precompute all F-grids and F_max values
t0 = time.time()
F_grids_remix = {}
for rmx in remix_values:
    F_grids_remix[rmx] = precompute_F_grid_remix(TM, Bd_remix, v_design,
                                                    sb_template_remix, rmx,
                                                    je_ref, ig_ref, v_grid, B_grid)
F_max_remix   = {rmx: float(np.abs(Fg).max()) for rmx, Fg in F_grids_remix.items()}
F_max_default = F_max_remix[remix_values[0]]   # explicit "ascending-sorted" anchor
print(f"Remix F-grid precompute ({len(remix_values)} grids): {time.time()-t0:.1f} s")

# Phase 2 — trajectory sweep + print table
print(f"\\nRemix sweep at (B_d={Bd_remix:.0f} G, s_sb={s_sb_remix:.0f}, L={L_remix} m, η={eta_r:.3f}):")
print(f"{'s_remix':>10s}  {'B_⊥ (G)':>10s}  {'Ω_L (MHz)':>11s}  {'F_max (1e-25 N)':>17s}  {'ratio':>10s}  {'v_capture (m/s)':>16s}")
results_remix = {}
t0 = time.time()
for rmx in remix_values:
    F_interp_r = RegularGridInterpolator(
        (v_grid, B_grid), F_grids_remix[rmx], method="linear",
        bounds_error=False, fill_value=0.0)
    for v0 in v0_remix_scan:
        sol = simulate_trajectory(v0, B_of_z_remix, F_interp_r, mass_RaF, L_remix,
                                    B_grid_lo=float(B_grid[0]),
                                    B_grid_hi=float(B_grid[-1]))
        results_remix[(rmx, v0)] = sol
    captured = _captured_at(results_remix, lambda v0, rmx=rmx: (rmx, v0), v0_remix_scan)
    v_cap = max(captured) if captured else None
    v_cap_str = f"{v_cap:6.0f}" if v_cap is not None else "    --"
    ratio = F_max_remix[rmx] / F_max_default
    B_perp_eq = s_remix_to_b_perp(rmx)
    Omega_L_MHz = b_perp_to_s_remix(B_perp_eq) * Gamma_FWHM_MHz / 2.0
    print(f"  {rmx:8.1e}  {B_perp_eq:10.3f}  {Omega_L_MHz:11.3f}  {F_max_remix[rmx]*1e25:>17.2f}  "
          f"{ratio:>10.2f}  {v_cap_str:>16s}")
print(f"Total remix trajectory sweep: {time.time()-t0:.1f} s")

# v(z) plot, overlay all remix values (unchanged)
remix_cmap = plt.get_cmap("plasma")
fig, axes = plt.subplots(1, 2, figsize=(15, 5))
ax_vz, ax_vt = axes
for ir, rmx in enumerate(remix_values):
    color = remix_cmap(ir / max(1, len(remix_values)-1))
    for v0 in v0_remix_scan:
        r = results_remix[(rmx, v0)]
        ax_vz.plot(r["z"], r["v"], "-", color=color, lw=0.8, alpha=0.6)
        ax_vt.plot(r["t"]*1e3, r["v"], "-", color=color, lw=0.8, alpha=0.6)
    # Add a labeled trace at the highest v_0 for legend
    r0 = results_remix[(rmx, v0_remix_scan[-1])]
    ax_vz.plot(r0["z"], r0["v"], "-", color=color, lw=2.0, label=f"s_remix={rmx:.0e}")
ax_vz.plot(z_p_r, v_t_r, "k--", lw=1.0, alpha=0.6, label=r"$v_{\\rm target}(z)$")
ax_vz.axhline(5, color="r", lw=0.5, alpha=0.4)
ax_vt.axhline(5, color="r", lw=0.5, alpha=0.4)
ax_vz.set_xlabel("z (m)"); ax_vz.set_ylabel("v (m/s)")
ax_vt.set_xlabel("t (ms)"); ax_vt.set_ylabel("v (m/s)")
ax_vz.set_title(f"v(z) — remix sweep at $B_d$={Bd_remix:.0f}, $s_{{\\\\rm sb}}$={s_sb_remix:.0f}, L={L_remix} m")
ax_vt.set_title("v(t) — same configs")
ax_vz.legend(fontsize=8, loc="best")
ax_vz.grid(alpha=0.3); ax_vt.grid(alpha=0.3)
fig.suptitle("Dark-state-mixing sweep — does lifting m_F dark states unlock capture?")
fig.tight_layout()
plt.show()
"""

ast.parse(CELL_46)
print(f"OK syntax cell 46: {len(CELL_46.splitlines())} lines")

nb = json.load(open(NB))
target = None
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] == "code" and "def precompute_F_grid_remix" in "".join(c.get("source", [])):
        target = i; break
assert target is not None, "s_remix sweep cell not found"
print(f"Targeting cell {target}")

nb["cells"][target]["source"] = CELL_46.splitlines(keepends=True)
nb["cells"][target]["outputs"] = []
nb["cells"][target]["execution_count"] = None

tmp = NB + ".tmp"
with open(tmp, "w") as f: json.dump(nb, f, indent=1)
os.replace(tmp, NB)
print(f"Cell {target} replaced. Total cells: {len(nb['cells'])}")
```

Run: `python3 /tmp/zs_p4_task2.py`
Expected: prints `OK syntax cell 46: ~85 lines`, `Targeting cell 46`.

- [ ] **Step 2: Re-execute notebook.**

```bash
conda run -n Structure jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=900 \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
```

- [ ] **Step 3: Verify per the Verify block at the top.** Run the inline Python check. Expected: `OK`.

- [ ] **Step 4: Spot-check the print table is unchanged.** Read cell 46's stdout from the notebook JSON and confirm the row format (`s_remix | B_⊥ | Ω_L | F_max | ratio | v_capture`) is intact. Compare to git HEAD~1 outputs:

```bash
python3 - <<'PY'
import json
nb = json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"))
for i, c in enumerate(nb["cells"]):
    if c["cell_type"]=="code" and "def precompute_F_grid_remix" in "".join(c.get("source", [])):
        for o in c.get("outputs", []):
            if o.get("output_type") == "stream" and o.get("name") == "stdout":
                print("".join(o.get("text", [])))
        break
PY
```

Expected: a multi-line table with the same column headers and 5 data rows for the 5 `remix_values`. The `ratio` column for `s_remix=1.0e-03` should still be `1.00` (baseline anchor unchanged).

- [ ] **Step 5: Commit.**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
/Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 4 Fix 3 — cell 46 two-phase split (precompute F-grids first, explicit F_max_default anchor); use _captured_at"
```

---

### Task 3: Env-var bless + git-tracked baseline (rename + cell 54 logic)

**Goal:** Rename `zs_baselines_pre_refactor.json` → `zs_baselines.json`; remove the `.gitignore` exemption and commit the file; rewrite cell 54 to read `ZS_BASELINE_BLESS` env var; update cell 55 path. Eliminates `_capture_range_local` from cell 54 (uses Task 1's `_captured_at` instead).

**Files:**
- Modify: `Jupyter Notebooks/RaX/zeeman_slower_2.ipynb` cells 54, 55
- Modify: `.gitignore` (remove line 16)
- Rename + commit: `Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json` → `Jupyter Notebooks/RaX/zs_baselines.json`

**Acceptance Criteria:**
- [ ] `Jupyter Notebooks/RaX/zs_baselines.json` exists, is tracked in git, contains the current sweep results.
- [ ] `Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json` does not exist on disk.
- [ ] `.gitignore` does not mention `zs_baselines_pre_refactor.json`.
- [ ] Cell 54 reads `os.environ.get('ZS_BASELINE_BLESS')`, prints the bless message and skips the check when set, falls through to current diff-and-raise logic when unset.
- [ ] Cell 54 calls `_captured_at` (from Task 1) instead of defining `_capture_range_local`.
- [ ] Cell 55's path string is `'zs_baselines.json'` (no `_pre_refactor`).
- [ ] Notebook re-executes clean (no env var → cell 54 PASS).
- [ ] **Bless smoke test:** running with `ZS_BASELINE_BLESS=1` produces `git diff` of zero on the JSON.

**Verify (main path):**
```bash
conda run -n Structure jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=900 \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
python3 - <<'PY'
import json, os
assert os.path.exists("Jupyter Notebooks/RaX/zs_baselines.json"), "renamed file missing"
assert not os.path.exists("Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json"), \
    "old filename still present"
gi = open(".gitignore").read()
assert "zs_baselines_pre_refactor.json" not in gi, ".gitignore still has old entry"
nb = json.load(open("Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"))
src_54 = src_55 = None
for i, c in enumerate(nb["cells"]):
    if c["cell_type"]!="code": continue
    s = "".join(c.get("source", []))
    if "regression check vs prior-run baseline" in s or "ZS_BASELINE_BLESS" in s:
        src_54 = s
    if "serialize sweep results as regression baseline" in s:
        src_55 = s
assert src_54 is not None and "ZS_BASELINE_BLESS" in src_54, "cell 54 not updated"
assert "_capture_range_local" not in src_54, "cell 54 still defines _capture_range_local"
assert "_captured_at(" in src_54, "cell 54 not using _captured_at"
assert src_55 is not None and "zs_baselines.json" in src_55, "cell 55 path not updated"
assert "zs_baselines_pre_refactor" not in src_55, "cell 55 still has old path"
errs = sum(1 for c in nb["cells"] if c["cell_type"]=="code"
            for o in c.get("outputs", []) if o.get("output_type")=="error")
assert errs == 0, f"{errs} error cells"
saw_pass = any("Regression check PASS" in "".join(o.get("text", []))
                for c in nb["cells"] for o in c.get("outputs", [])
                if o.get("output_type") == "stream")
assert saw_pass, "regression PASS line not found"
print("OK")
PY
```
Expected: `OK`.

**Steps:**

- [ ] **Step 1: Rename file on disk.**

```bash
mv "Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json" \
   "Jupyter Notebooks/RaX/zs_baselines.json"
```
If the source file doesn't exist (e.g., never generated locally), skip — the next notebook execute will create it.

- [ ] **Step 2: Update `.gitignore` line 16.**

Use the Edit tool:
- old_string: `Jupyter Notebooks/RaX/zs_baselines_pre_refactor.json` (with surrounding context: the line above and below it to ensure uniqueness — read `.gitignore` first to capture exact context).
- Replace it with nothing (delete the line).

If preferred, use `sed` (BSD sed on macOS — note the `-i ''` requirement):
```bash
sed -i '' '/^Jupyter Notebooks\/RaX\/zs_baselines_pre_refactor\.json$/d' .gitignore
grep "zs_baselines" .gitignore && echo "ERROR: still present" || echo "OK removed"
```
Expected: `OK removed`.

- [ ] **Step 3: Edit cells 54 and 55 via Python script.** Save to `/tmp/zs_p4_task3.py`:

```python
"""Phase 4 Fix 2 — env-var bless + path rename in cells 54, 55."""
import json, ast, os

NB = "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"

CELL_54 = """\
# === Phase 4 Fix 2: regression check vs git-tracked baseline =============================
# Loads zs_baselines.json (committed in git, updated only via the bless workflow below)
# and diffs current sweeps. Drift in F_max (rtol=1e-6) or capture-velocity sets fails fast.
#
# To bless an intentional physics change:
#   1. Run normally; cell 54 FAILS with deltas printed. Confirm intentional.
#   2. ZS_BASELINE_BLESS=1 conda run -n Structure jupyter nbconvert --to notebook \\\\
#          --execute --inplace "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
#   3. git diff "Jupyter Notebooks/RaX/zs_baselines.json"
#   4. git add ... && git commit -m "physics: <reason>"
#
# Mirrors the pytest-regressions / Jest --updateSnapshot / cargo-insta-accept pattern.
import json as _json, os as _os

_baseline_path = "zs_baselines.json"
_bless = bool(_os.environ.get("ZS_BASELINE_BLESS"))

if _bless:
    print("⚡ ZS_BASELINE_BLESS set — skipping regression check; cell 55 will write the new baseline.")
elif not _os.path.exists(_baseline_path):
    print(f"⚠ {_baseline_path} not found — skipping regression check (cell 55 will create it).")
else:
    with open(_baseline_path) as _f: _baseline = _json.load(_f)
    _deltas = []
    _RTOL = 1e-6

    def _check(label, expected_F, current_F, expected_cap, current_cap):
        for k, v_exp in expected_F.items():
            v_cur = current_F.get(k)
            if v_cur is None:
                _deltas.append(f"{label} F_max[{k}]: missing in current")
            elif abs(v_cur - v_exp) > _RTOL * max(abs(v_exp), 1e-30):
                _deltas.append(f"{label} F_max[{k}]: {v_exp:.6e} -> {v_cur:.6e}")
        for k, c_exp in expected_cap.items():
            c_cur = current_cap.get(k, [])
            if list(c_exp) != list(c_cur):
                _deltas.append(f"{label} captures[{k}]: {c_exp} -> {c_cur}")

    if "single_sb" in _baseline and "results_single" in dir():
        cur_F = {f"{Bd:.0f}_{s:.0f}": float(np.abs(F_surfaces[(Bd, s)]).max())
                  for Bd in B_designs for s in s_values}
        cur_C = {f"{Bd:.0f}_{s:.0f}_{L:.2f}":
                    _captured_at(results_single, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)
                  for Bd in B_designs for s in s_values for L in L_values}
        _check("single_sb", _baseline["single_sb"]["F_max"], cur_F,
                _baseline["single_sb"]["captures"], cur_C)

    if "multi_sb" in _baseline and "results_multi" in dir():
        cur_F = {f"{Bd:.0f}_{s:.0f}": float(np.abs(F_surfaces_multi[(Bd, s)]).max())
                  for Bd in B_designs for s in s_values}
        cur_C = {f"{Bd:.0f}_{s:.0f}_{L:.2f}":
                    _captured_at(results_multi, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)
                  for Bd in B_designs for s in s_values for L in L_values}
        _check("multi_sb", _baseline["multi_sb"]["F_max"], cur_F,
                _baseline["multi_sb"]["captures"], cur_C)

    if "remix_sweep" in _baseline and "results_remix" in dir():
        cur_F = {f"{rmx:.0e}": float(F_max_remix[rmx]) for rmx in remix_values}
        cur_C = {f"{rmx:.0e}":
                    _captured_at(results_remix, lambda v0,rmx=rmx: (rmx,v0), v0_remix_scan)
                  for rmx in remix_values}
        _check("remix_sweep", _baseline["remix_sweep"]["F_max"], cur_F,
                _baseline["remix_sweep"]["captures"], cur_C)

    if "remix_full" in _baseline and "results_remix1" in dir():
        cur_F = {f"{Bd:.0f}_{s:.0f}": float(np.abs(F_surfaces_remix1[(Bd, s)]).max())
                  for Bd in B_designs for s in s_values}
        cur_C = {f"{Bd:.0f}_{s:.0f}_{L:.2f}":
                    _captured_at(results_remix1, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)
                  for Bd in B_designs for s in s_values for L in L_values}
        _check("remix_full", _baseline["remix_full"]["F_max"], cur_F,
                _baseline["remix_full"]["captures"], cur_C)

    if "pb" in _baseline and "results_pb" in dir():
        cur_F = {f"{s:.0f}": float(np.abs(F_surfaces_pb[s]).max()) for s in s_values}
        cur_C = {f"{s:.0f}_{L:.2f}":
                    _captured_at(results_pb, lambda v0,s=s,L=L: (s,L,v0), v0_scan)
                  for s in s_values for L in L_values}
        _check("pb", _baseline["pb"]["F_max"], cur_F,
                _baseline["pb"]["captures"], cur_C)

    if _deltas:
        print("✗ Regression check FAIL:")
        for _d in _deltas: print("  ", _d)
        raise AssertionError(f"{len(_deltas)} sweep deltas vs git-tracked baseline")
    else:
        print(f"✓ Regression check PASS: all sweeps match baseline within rtol={_RTOL}.")
"""

CELL_55 = """\
# === Phase 1f: serialize sweep results as regression baseline ============================
# All capture-velocity range tables and per-config F_max are written to JSON.
# Phase 4 Fix 2 made the file git-tracked — see cell 54 for the bless workflow.
import json as _json

_baseline = {}

# 1. Single-σ+ sweep: results_single[(Bd, s, L, v0)]
if 'results_single' in dir():
    _baseline["single_sb"] = {
        "F_max":   {f"{Bd:.0f}_{s:.0f}": float(np.abs(F_surfaces[(Bd, s)]).max())
                    for Bd in B_designs for s in s_values},
        "captures": {f"{Bd:.0f}_{s:.0f}_{L:.2f}":
                       _captured_at(results_single, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)
                     for Bd in B_designs for s in s_values for L in L_values},
    }

# 2. Multi-σ+ sweep: results_multi
if 'results_multi' in dir():
    _baseline["multi_sb"] = {
        "F_max":    {f"{Bd:.0f}_{s:.0f}": float(np.abs(F_surfaces_multi[(Bd, s)]).max())
                     for Bd in B_designs for s in s_values},
        "captures": {f"{Bd:.0f}_{s:.0f}_{L:.2f}":
                       _captured_at(results_multi, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)
                     for Bd in B_designs for s in s_values for L in L_values},
    }

# 3. Remix sweep: results_remix[(rmx, v0)]
if 'results_remix' in dir():
    _baseline["remix_sweep"] = {
        "F_max":    {f"{rmx:.0e}": float(F_max_remix[rmx]) for rmx in remix_values},
        "captures": {f"{rmx:.0e}":
                       _captured_at(results_remix, lambda v0,rmx=rmx: (rmx,v0), v0_remix_scan)
                     for rmx in remix_values},
    }

# 4. Full-remix=1.0 sweep: results_remix1
if 'results_remix1' in dir():
    _baseline["remix_full"] = {
        "F_max":    {f"{Bd:.0f}_{s:.0f}": float(np.abs(F_surfaces_remix1[(Bd, s)]).max())
                     for Bd in B_designs for s in s_values},
        "captures": {f"{Bd:.0f}_{s:.0f}_{L:.2f}":
                       _captured_at(results_remix1, lambda v0,Bd=Bd,s=s,L=L: (Bd,s,L,v0), v0_scan)
                     for Bd in B_designs for s in s_values for L in L_values},
    }

# 5. Paschen-Back σ⁻ sweep: results_pb[(s_per_sb, L, v0)]
if 'results_pb' in dir():
    _baseline["pb"] = {
        "F_max":    {f"{s:.0f}": float(np.abs(F_surfaces_pb[s]).max()) for s in s_values},
        "captures": {f"{s:.0f}_{L:.2f}":
                       _captured_at(results_pb, lambda v0,s=s,L=L: (s,L,v0), v0_scan)
                     for s in s_values for L in L_values},
    }

_path = 'zs_baselines.json'
with open(_path, 'w') as _f: _json.dump(_baseline, _f, indent=2, sort_keys=True)
print(f"Wrote {_path}: {sorted(_baseline.keys())}")
"""

ast.parse(CELL_54)
ast.parse(CELL_55)
print(f"OK syntax cell 54 ({len(CELL_54.splitlines())} lines), cell 55 ({len(CELL_55.splitlines())} lines)")

nb = json.load(open(NB))
target_54 = target_55 = None
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code": continue
    s = "".join(c.get("source", []))
    if "regression check vs prior-run baseline" in s:
        target_54 = i
    if "serialize sweep results as regression baseline" in s:
        target_55 = i
assert target_54 is not None and target_55 is not None, \
    f"cells not found: 54={target_54}, 55={target_55}"
print(f"Targeting cells {target_54} (regression) and {target_55} (serializer)")

nb["cells"][target_54]["source"] = CELL_54.splitlines(keepends=True)
nb["cells"][target_55]["source"] = CELL_55.splitlines(keepends=True)
for i in (target_54, target_55):
    nb["cells"][i]["outputs"] = []
    nb["cells"][i]["execution_count"] = None

tmp = NB + ".tmp"
with open(tmp, "w") as f: json.dump(nb, f, indent=1)
os.replace(tmp, NB)
print(f"Cells {target_54} + {target_55} replaced. Total cells: {len(nb['cells'])}")
```

Run: `python3 /tmp/zs_p4_task3.py`
Expected: prints `OK syntax ...`, `Targeting cells 54 (regression) and 55 (serializer)`.

- [ ] **Step 4: Re-execute notebook (main path, no env var).**

```bash
conda run -n Structure jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.timeout=900 \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
```
Expected: prints `[NbConvertApp] Writing ...`. No tracebacks.

- [ ] **Step 5: Verify per the Verify (main path) block.** Run the inline Python check. Expected: `OK`.

- [ ] **Step 6: Stage everything.**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb" \
    "Jupyter Notebooks/RaX/zs_baselines.json" \
    .gitignore
/Library/Developer/CommandLineTools/usr/bin/git status -s
```
Expected status (paraphrased): `M  .gitignore`, `M  Jupyter Notebooks/RaX/zeeman_slower_2.ipynb`, `A  Jupyter Notebooks/RaX/zs_baselines.json`. The old filename should NOT appear (the `mv` removes it from the working tree, and since it was previously gitignored, git is unaware of the deletion — there's no `D` line for it).

- [ ] **Step 7: Commit.**

```bash
/Library/Developer/CommandLineTools/usr/bin/git commit -m "zeeman_slower_2: Phase 4 Fix 2 — env-var bless + git-tracked baseline (rename to zs_baselines.json; cell 54 reads ZS_BASELINE_BLESS; reuse _captured_at)"
```

- [ ] **Step 8: Bless smoke test.** Re-execute with env var set; confirm zero diff on the JSON.

```bash
ZS_BASELINE_BLESS=1 conda run -n Structure jupyter nbconvert --to notebook \
    --execute --inplace --ExecutePreprocessor.timeout=900 \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
```
Expected: cell 54 stdout includes `⚡ ZS_BASELINE_BLESS set — skipping regression check`.

```bash
/Library/Developer/CommandLineTools/usr/bin/git diff --stat \
    "Jupyter Notebooks/RaX/zs_baselines.json"
```
Expected: empty output (no diff).

If diff is non-empty: investigate. The bless smoke test should be reproducible (same code → same baseline values). A non-empty diff means cell 55 produced a different value than what's committed — possibly FP non-determinism or an earlier change in the notebook that affected the sweeps.

- [ ] **Step 9: Discard the re-execution outputs from the smoke test.** The notebook itself was modified by the smoke-test execute (new exec counts, new image hashes); revert those.

```bash
/Library/Developer/CommandLineTools/usr/bin/git checkout \
    "Jupyter Notebooks/RaX/zeeman_slower_2.ipynb"
```
Expected: working tree clean. Verify with `git status` — should show no modifications.

---

## Self-Review

**Spec coverage (3 sections of the spec, each maps to a task):**
- Spec §"Three-fix structure / Fix 1" → Task 1 ✓
- Spec §"Three-fix structure / Fix 2" → Task 3 ✓
- Spec §"Three-fix structure / Fix 3" → Task 2 ✓
- Spec §"Sequencing" (1 → 3 → 2 commit order) → Tasks ordered 1, 2, 3 in this plan with commit messages reflecting Fix 1 / Fix 3 / Fix 2 ✓
- Spec §"Verification" (no `--allow-errors`, error-cell count == 0, regression PASS) → each task's Verify block ✓
- Spec §"Out of scope" → no tasks; nothing to map.
- Spec §"Pattern provenance" → embedded as comment in Task 3's cell 54 source ✓

**Placeholder scan:** No "TBD"/"TODO"/"fill in"/"similar to" patterns. Code blocks complete in every step. The spec's `<reason>` placeholder in the bless workflow is a user-input slot, not a missing plan piece.

**Type/identifier consistency:**
- `_captured_at(results, key_fn, v0_iter)` defined in Task 1, used in Task 2 and Task 3 with identical signature ✓
- `F_max_default` used inside cell 46 in Task 2; computed before Phase 2 loop ✓
- `_baseline_path = "zs_baselines.json"` consistent between cells 54 and 55 in Task 3 ✓
- Cell 54's removal of `_capture_range_local` matches the ban in the Verify check ✓

**Bonus improvement found during plan-write:** Cell 54 had its own `_capture_range_local` that's identical to Task 1's `_captured_at`. Folded into Task 3 (eliminate duplicate, use the global helper). Mentioned in spec § "Pattern provenance" as the natural consequence of having a shared primitive.
