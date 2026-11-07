"""Mutation testing of the revised analyze.py against the revised test_analyze.py (author's harness, not part of the
study). Usage: python3 -I mutate2.py CODE_DIR OUT_DIR"""
import os
import shutil
import subprocess
import sys

CODE, OUT = sys.argv[1], sys.argv[2]
SRC = open(os.path.join(CODE, "analyze.py")).read()

M = [
    ("NEED7", "NEED = 8 ", "NEED = 7 "),
    ("NEED9", "NEED = 8 ", "NEED = 9 "),
    ("hits>need", "verdict(hits >= need))", "verdict(hits > need))"),
    ("gate>need", 'gate = {key: acc[key]["correct"] >= need', 'gate = {key: acc[key]["correct"] > need'),
    ("H3_or", 'v3 = NOT_MET if NOT_MET in clause_verdicts else', 'v3 = NOT_MET if clause_verdicts.count(NOT_MET) == 2 else'),
    ("H3_testable_first", 'v3 = NOT_MET if NOT_MET in clause_verdicts else NOT_TESTABLE if NOT_TESTABLE in clause_verdicts else MET',
     'v3 = NOT_TESTABLE if NOT_TESTABLE in clause_verdicts else NOT_MET if NOT_MET in clause_verdicts else MET'),
    ("H3_only66", 'clause_verdicts = (h3a["verdict"], h3b["verdict"])', 'clause_verdicts = (h3a["verdict"],)'),
    ("H3_gate_std_ring", '("ring_yaw_inc66",)), shape_hyp("yaw_inc0", "ring", ("ring_yaw_inc0",))',
     '("ring",)), shape_hyp("yaw_inc0", "ring", ("ring",))'),
    ("H3_gates_swapped", '("ring_yaw_inc66",)), shape_hyp("yaw_inc0", "ring", ("ring_yaw_inc0",))',
     '("ring_yaw_inc0",)), shape_hyp("yaw_inc0", "ring", ("ring_yaw_inc66",))'),
    ("H3_why_empty", 'why3 = h3a["not_testable_because"] + h3b["not_testable_because"] if v3 == NOT_TESTABLE else []', "why3 = []"),
    ("H3_why_always", 'why3 = h3a["not_testable_because"] + h3b["not_testable_because"] if v3 == NOT_TESTABLE else []',
     'why3 = h3a["not_testable_because"] + h3b["not_testable_because"]'),
    ("gate_target_torus", 'GATE_RINGS = tuple(f"ring_{c}" for c in YAW)', 'GATE_RINGS = tuple(f"ring_{c}" for c in YAW[::-1])'),
    ("atlas_target_torus", '**{g: "ring" for g in GATE_RINGS}}', '**{g: "torus" for g in GATE_RINGS}}'),
    ("Q1_with_gates", "verdict=verdict(all(gate[k] for k in SHAPES)))", "verdict=verdict(all(gate[k] for k in ATLAS_KEYS)))"),
    ("blas_none_ok", 'if meta.get("blas_threads") != [1]:', 'if meta.get("blas_threads") not in ([1], None):'),
    ("blas_any_one", 'if meta.get("blas_threads") != [1]:', 'if 1 not in (meta.get("blas_threads") or [1]):'),
    ("null_only", 'if (params.get("rel_factor"), params.get("null_factor")) != FROZEN_FACTORS:',
     'if params.get("rel_factor") != FROZEN_FACTORS[0]:'),
    ("rel_only", 'if (params.get("rel_factor"), params.get("null_factor")) != FROZEN_FACTORS:',
     'if params.get("null_factor") != FROZEN_FACTORS[1]:'),
    ("load_unreadable_skip", 'problems.append(f"{name}: unreadable")', 'pass'),
    ("main_no_provisional", '"" if res["run_valid"] else "(provisional)")', '"")'),
    ("H4_eq_upper", "t_ci(c0)[0] >= CONF_EQUATOR_MIN", "t_ci(c0)[1] >= CONF_EQUATOR_MIN"),
    ("H4_66_lower", "t_ci(c66)[1] <= CONF_66_MAX", "t_ci(c66)[0] <= CONF_66_MAX"),
    ("H4_eq_mean", "bool(t_ci(c0)[0] >= CONF_EQUATOR_MIN)", "bool(np.mean(c0) >= CONF_EQUATOR_MIN)"),
    ("H4_66_mean", "bool(t_ci(c66)[1] <= CONF_66_MAX)", "bool(np.mean(c66) <= CONF_66_MAX)"),
    ("H4_or", "verdict=verdict(all(clauses.values())))", "verdict=verdict(any(clauses.values())))"),
    ("EQ_MIN_0.45", "CONF_EQUATOR_MIN, CONF_66_MAX = 0.40, 0.05", "CONF_EQUATOR_MIN, CONF_66_MAX = 0.45, 0.05"),
    ("EQ_MIN_0.32", "CONF_EQUATOR_MIN, CONF_66_MAX = 0.40, 0.05", "CONF_EQUATOR_MIN, CONF_66_MAX = 0.32, 0.05"),
    ("MAX66_0.03", "CONF_EQUATOR_MIN, CONF_66_MAX = 0.40, 0.05", "CONF_EQUATOR_MIN, CONF_66_MAX = 0.40, 0.03"),
    ("MAX66_0.065", "CONF_EQUATOR_MIN, CONF_66_MAX = 0.40, 0.05", "CONF_EQUATOR_MIN, CONF_66_MAX = 0.40, 0.065"),
    ("CI_level_0.90", "def t_ci(v, level=0.95):", "def t_ci(v, level=0.90):"),
    ("CI_level_0.99", "def t_ci(v, level=0.95):", "def t_ci(v, level=0.99):"),
    ("CI_ddof0", "v.std(ddof=1)", "v.std(ddof=0)"),
    ("CI_normal", "stats.t.ppf(0.5 + level / 2.0, len(v) - 1)", "stats.norm.ppf(0.5 + level / 2.0)"),
    ("needed_const", "return math.ceil(NEED / 10 * n_seeds)", "return NEED"),
    ("needed_floor", "return math.ceil(NEED / 10 * n_seeds)", "return math.floor(NEED / 10 * n_seeds)"),
    ("I1_phys_1e-7", "antipodal_physics=1e-12", "antipodal_physics=1e-7"),
    ("I1_pop_1e-7", "antipodal_population=1e-9", "antipodal_population=1e-7"),
    ("I1_or", 'within(physics.get("antipodal_max_diff"), TOL["antipodal_physics"]) and',
     'within(physics.get("antipodal_max_diff"), TOL["antipodal_physics"]) or'),
    ("I1_any", 'all(within(a, TOL["antipodal_population"]) for a in antip)', 'any(within(a, TOL["antipodal_population"]) for a in antip)'),
    ("I2_1e-7", "odd_power=1e-10", "odd_power=1e-7"),
    ("I3_radicalpy_5e-3", "radicalpy=1e-6", "radicalpy=5e-3"),
    ("I3_resolvent_1e-6", "resolvent=1e-10", "resolvent=1e-6"),
    ("I3_zeeman_1e-6", "zeeman=1e-12,", "zeeman=1e-6,"),
    ("I3_any", 'verdict="ok" if all(i3.values()) else "bug")', 'verdict="ok" if any(i3.values()) else "bug")'),
    ("within_missing_ok", "return is_num(x) and x <= tol", "return (not is_num(x)) or x <= tol"),
    ("is_num_nan_ok", "return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)",
     "return isinstance(x, (int, float)) and not isinstance(x, bool)"),
    ("all_finite_true", "        return math.isfinite(x)\n    return True", "        return True\n    return True"),
    ("I5_finite_default", 'if rec.get("finite") is not True:', 'if rec.get("finite", True) is not True:'),
    ("I5_no_topo_keys", "missing = [k for k in TOPO_KEYS if k not in rec]", "missing = []"),
    ("I5_no_label", 'if rec["shape"] not in LABELS:', "if False:"),
    ("I5_no_finite_bars", 'if not all_finite([rec["bars"], rec["null"], rec["thresholds"]]):', "if False:"),
    ("I5_no_yaw", "if yaw and not (", "if False and not ("),
    ("I5_skip_gate_rings", "for key in ATLAS_KEYS:\n            if key not in atlas:", "for key in SHAPES:\n            if key not in atlas:"),
    ("I5_skip_design", "for key in DESIGN_KEYS:\n            if as_dict", "for key in DESIGN_KEYS[:0]:\n            if as_dict"),
    ("I5_skip_last_cond", "for name in CONDITIONS:", "for name in CONDITIONS[:-1]:"),
    ("I5_no_antip", 'if not is_num(s.get("noiseless_antipodal_max_diff")):', "if False:"),
    ("I5_no_early_return", '    if problems:\n        for k in ("Q1"', '    if False:\n        for k in ("Q1"'),
    ("I5_extra_ignored", "problems = list(extra_problems) + structure_problems(seeds)", "problems = structure_problems(seeds)"),
    ("I5_no_seeds_ok", 'return ["no seed records"]', "return []"),
    ("Q1_hyp_only", "verdict=verdict(all(gate[k] for k in SHAPES)))",
     'verdict=verdict(all(gate[k] for k in ("point", "ring", "sphere", "rp2"))))'),
    ("Q1_any", "verdict=verdict(all(gate[k] for k in SHAPES)))", "verdict=verdict(any(gate[k] for k in SHAPES)))"),
    ("Q2_0.006", "anisotropy_min=0.01", "anisotropy_min=0.006"),
    ("Q3_0.015", "fit_oos_max=0.01)", "fit_oos_max=0.015)"),
    ("Q3_in_sample", 'physics.get("fit_oos_rms_over_sd")', 'physics.get("fit_rms_over_sd")'),
    ("H1_cond_d025", 'shape_hyp("rp_d0", "rp2", ("rp2",))', 'shape_hyp("rp_d0.25", "rp2", ("rp2",))'),
    ("H1_gate_also_ring", 'shape_hyp("rp_d0", "rp2", ("rp2",))', 'shape_hyp("rp_d0", "rp2", ("rp2", "ring"))'),
    ("H2_gate_also_ring", 'shape_hyp("polar", "sphere", ("sphere",))', 'shape_hyp("polar", "sphere", ("sphere", "ring"))'),
    ("run_valid_true", 'res["run_valid"] = bool(all(res[i]["verdict"] == "ok" for i in INVARIANTS) and res["provenance"]["ok"])',
     'res["run_valid"] = True'),
    ("run_valid_no_prov", 'res["run_valid"] = bool(all(res[i]["verdict"] == "ok" for i in INVARIANTS) and res["provenance"]["ok"])',
     'res["run_valid"] = bool(all(res[i]["verdict"] == "ok" for i in INVARIANTS))'),
    ("INVARIANTS_no_I5", 'INVARIANTS = ("I1", "I2", "I3", "I4", "I5")', 'INVARIANTS = ("I1", "I2", "I3", "I4")'),
    ("confirmatory_any_seeds",
     'res["confirmatory"] = bool(res["run_valid"] and got == CONFIRMATORY_SEEDS and len(seeds) == len(CONFIRMATORY_SEEDS))',
     'res["confirmatory"] = bool(res["run_valid"])'),
    ("prov_no_code", '(("code_sha256", "mag_s2.py"), ("analysis_sha256"', '(("analysis_sha256"'),
    ("prov_no_analysis", '("analysis_sha256", "analyze.py"), ', ""),
    ("prov_no_physics", ', ("physics_sha256", "physics.json")):', "):"),
    ("prov_no_factors", 'if (params.get("rel_factor"), params.get("null_factor")) != FROZEN_FACTORS:', "if False:"),
    ("prov_no_blas", 'if meta.get("blas_threads") != [1]:', "if False:"),
    ("prov_no_seeds", 'if sorted(meta.get("seeds") or []) != got:', "if False:"),
    ("prov_no_mode", 'if meta.get("mode") != "run":', "if False:"),
    ("load_no_name_check", 'if m is None or s.get("seed") != int(m.group(1)):', "if False:"),
    ("E1_wrong_key", 'counts([s["conditions"][f"rp_d{d:g}"]["shape"] for s in seeds]) for d in D_LEVELS',
     'counts([s["conditions"]["rp_d0"]["shape"] for s in seeds]) for d in D_LEVELS'),
    ("E2_drop", 'res["E2_design"] = {k: counts([s["design"][k]["shape"] for s in seeds]) for k in DESIGN_KEYS}',
     'res["E2_design"] = {k: counts([s["design"][k]["shape"] for s in seeds]) for k in DESIGN_KEYS[:9]}'),
    ("E3_inverse", "float(bars[0] / tau)", "float(tau / bars[0])"),
    ("E3_no_gate", 'margins(lambda s, k=k: s["atlas"][k]) for k in ATLAS_KEYS}', 'margins(lambda s, k=k: s["atlas"][k]) for k in SHAPES}'),
    ("main_no_notconf", '    elif not res["confirmatory"]:', "    elif False:"),
    ("main_no_invalid", '    if not res["run_valid"]:\n        bad', "    if False:\n        bad"),
]


def run_one(mid, old, new):
    n = SRC.count(old)
    if n != 1:
        return mid, f"PATTERN COUNT {n}"
    d = os.path.join(OUT, mid)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "analyze.py"), "w") as fh:
        fh.write(SRC.replace(old, new))
    for f in ("test_analyze.py", "mag_s2.py", "physics.json"):
        shutil.copy(os.path.join(CODE, f), d)
    try:
        p = subprocess.run([sys.executable, "-I", "-B", os.path.join(d, "test_analyze.py")], capture_output=True, text=True,
                           cwd=d, timeout=300)
    except subprocess.TimeoutExpired:
        return mid, "killed (timeout)"
    fails = [ln[5:] for ln in p.stdout.splitlines() if ln.startswith("FAIL ")]
    status = "SURVIVED" if p.returncode == 0 else ("killed" if fails else f"killed (crash rc={p.returncode})")
    return mid, status + (f"  by: {fails[0]}" + (f" (+{len(fails) - 1})" if len(fails) > 1 else "") if fails else "")


if __name__ == "__main__":
    shutil.rmtree(OUT, ignore_errors=True)
    rows = [run_one(*m) for m in M]
    for mid, st in rows:
        print("%-24s %s" % (mid, st), flush=True)
    surv = [m for m, st in rows if st.startswith("SURVIVED") or st.startswith("PATTERN")]
    print(f"\n{len(rows)} mutations, {len(rows) - len(surv)} killed, survivors or bad patterns: {surv}")
