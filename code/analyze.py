#!/usr/bin/env python3
"""
MAG-S2 confirmatory analysis, protocol version 0.1.

  I1-I5   code and physics invariants: "bug" invalidates the run
  Q1-Q3   expectations: the atlas check of the pipeline in the run itself and the physics of the compass
  H1-H4   hypotheses: "met", "not met" or "not testable"

A seed is one simulated bird (its own retina, neurons, field directions and noise) together with one atlas of known
shapes drawn at the same settings. A hypothesis about a shape is testable only when the pipeline recognises that shape
in the atlas of the same run in at least 8 of 10 seeds. For H3 the gate is a pair of atlas rings whose neurons carry
the noise share of the two yaw conditions, one ring per clause; H3 is "not met" as soon as one testable clause fails,
and "not testable" when a gate is closed and no testable clause fails.

A run is valid when every invariant holds and its metadata match the files next to this script: the hashes of
mag_s2.py, analyze.py and physics.json, the calibrated threshold factors, one BLAS thread and the seeds the run was
started with. A valid run is confirmatory only when its seed files are exactly seeds 1-10.

Usage
  python3 -I analyze.py RUN_DIR      (reads physics.json and mag_s2.py next to this file and RUN_DIR/meta_run.json;
                                      writes RUN_DIR/results.json)
"""
from __future__ import annotations

import glob
import hashlib
import json
import math
import os
import re
import sys

import numpy as np
from scipy import stats

MET, NOT_MET, NOT_TESTABLE = "met", "not met", "not testable"
CONFIRMATORY_SEEDS = list(range(1, 11))
SHAPES = ("point", "ring", "sphere", "torus", "rp2", "klein")
LABELS = SHAPES + ("other",)
YAW = ("yaw_inc66", "yaw_inc0")
GATE_RINGS = tuple(f"ring_{c}" for c in YAW)          # atlas rings carrying the noise share of each yaw condition
ATLAS_KEYS = SHAPES + GATE_RINGS
ATLAS_TARGET = {**{s: s for s in SHAPES}, **{g: "ring" for g in GATE_RINGS}}
D_LEVELS = (0.0, 0.25, 0.5, 0.75, 1.0)
CONDITIONS = tuple(f"rp_d{d:g}" for d in D_LEVELS) + ("rp_above", "polar") + YAW
DESIGN_KEYS = tuple(f"N{n}_noise{z:g}_{k}" for n in (25, 50, 100) for z in (0.25, 0.5, 1.0) for k in ("rp", "polar"))
TOPO_KEYS = ("n_points", "bars", "null", "thresholds", "betti", "signature", "shape")
NEED = 8                       # seeds out of 10 for shape claims and for the atlas gates
CONF_EQUATOR_MIN, CONF_66_MAX = 0.40, 0.05
FROZEN_FACTORS = (0.4, 4.0)    # rel_factor and null_factor chosen by the calibration on the atlas
TOL = dict(antipodal_physics=1e-12, antipodal_population=1e-9, odd_power=1e-10, radicalpy=1e-6, resolvent=1e-10,
           zeeman=1e-12, anisotropy_min=0.01, fit_oos_max=0.01)
ORDER = ["I1", "I2", "I3", "I4", "I5", "Q1", "Q2", "Q3", "H1", "H2", "H3", "H4"]
INVARIANTS = ("I1", "I2", "I3", "I4", "I5")


def verdict(ok):
    return MET if ok else NOT_MET


def is_num(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def within(x, tol):
    """x is a finite number no larger than tol; a missing or non-finite value fails."""
    return is_num(x) and x <= tol


def all_finite(x):
    """Every number inside nested dicts and lists is finite."""
    if isinstance(x, dict):
        return all(all_finite(v) for v in x.values())
    if isinstance(x, (list, tuple)):
        return all(all_finite(v) for v in x)
    if isinstance(x, (int, float)) and not isinstance(x, bool):
        return math.isfinite(x)
    return True


def as_dict(x):
    return x if isinstance(x, dict) else {}


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def t_ci(v, level=0.95):
    v = np.asarray(v, float)
    if len(v) < 2 or not np.all(np.isfinite(v)):
        return [float("nan"), float("nan")]
    h = stats.t.ppf(0.5 + level / 2.0, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
    return [float(v.mean() - h), float(v.mean() + h)]


def needed(n_seeds):
    return math.ceil(NEED / 10 * n_seeds)


def counts(labels):
    return {lab: labels.count(lab) for lab in sorted(set(labels), key=str)}


# --------------------------------------------------------------------------------------------------------------
# Completeness of the records (I5)
# --------------------------------------------------------------------------------------------------------------
def record_problems(where, rec, topo=True, yaw=False):
    if not isinstance(rec, dict):
        return [f"{where}: missing"]
    out = []
    if rec.get("finite") is not True:
        out.append(f"{where}: population not finite")
    if topo:
        missing = [k for k in TOPO_KEYS if k not in rec]
        if missing:
            out.append(f"{where}: missing {', '.join(missing)}")
        else:
            if rec["shape"] not in LABELS:
                out.append(f"{where}: unknown label {rec['shape']!r}")
            if not all_finite([rec["bars"], rec["null"], rec["thresholds"]]):
                out.append(f"{where}: non-finite bars or thresholds")
    if yaw and not (is_num(rec.get("confusion_180")) and 0.0 <= rec["confusion_180"] <= 1.0):
        out.append(f"{where}: heading confusion missing or outside [0, 1]")
    return out


def structure_problems(seeds):
    if not seeds:
        return ["no seed records"]
    out = []
    for s in seeds:
        sid = s.get("seed")
        cond, atlas, design = as_dict(s.get("conditions")), as_dict(s.get("atlas")), as_dict(s.get("design"))
        if not is_num(s.get("noiseless_antipodal_max_diff")):
            out.append(f"seed {sid}: noiseless antipodal difference missing or not finite")
        for name in CONDITIONS:
            if name not in cond:
                out.append(f"seed {sid}: condition {name} missing")
            elif name == "rp_above":
                out += record_problems(f"seed {sid} rp_above", cond[name], topo=False)
            else:
                out += record_problems(f"seed {sid} {name}", cond[name], yaw=name in YAW)
        for key in ATLAS_KEYS:
            if key not in atlas:
                out.append(f"seed {sid}: atlas {key} missing")
            else:
                out += record_problems(f"seed {sid} atlas {key}", atlas[key])
        for key in DESIGN_KEYS:
            if as_dict(design.get(key)).get("shape") not in LABELS:
                out.append(f"seed {sid}: design point {key} missing or unlabelled")
    return out


# --------------------------------------------------------------------------------------------------------------
# Verdicts
# --------------------------------------------------------------------------------------------------------------
def analyze(seeds, physics, extra_problems=()):
    res = {}
    n = len(seeds)
    need = needed(n)

    # ---------------- invariants
    antip = [s.get("noiseless_antipodal_max_diff") for s in seeds]
    ok1 = within(physics.get("antipodal_max_diff"), TOL["antipodal_physics"]) and \
        all(within(a, TOL["antipodal_population"]) for a in antip)
    res["I1"] = dict(what="the radical-pair yield and the noiseless population are identical for n and -n",
                     physics_max_diff=physics.get("antipodal_max_diff"), population_per_seed=antip,
                     verdict="ok" if ok1 else "bug")
    ok2 = within(physics.get("odd_power_fraction"), TOL["odd_power"])
    res["I2"] = dict(what="the spherical-harmonic fit adds no odd degrees (share of odd degrees in the anisotropic power)",
                     odd_power_fraction=physics.get("odd_power_fraction"), verdict="ok" if ok2 else "bug")
    i3 = dict(radicalpy=within(physics.get("radicalpy_max_diff"), TOL["radicalpy"]),
              resolvent=within(physics.get("resolvent_max_diff"), TOL["resolvent"]),
              zeeman=within(physics.get("zeeman_cartesian_max_rel_diff"), TOL["zeeman"]))
    res["I3"] = dict(what="the exact yield agrees with a time-domain propagation on RadicalPy's propagator and with a "
                          "resolvent solution, and the Zeeman term agrees with its Cartesian form",
                     radicalpy_max_diff=physics.get("radicalpy_max_diff"), resolvent_max_diff=physics.get("resolvent_max_diff"),
                     zeeman_cartesian_max_rel_diff=physics.get("zeeman_cartesian_max_rel_diff"), clauses=i3,
                     verdict="ok" if all(i3.values()) else "bug")
    mism = 0
    for s in seeds:
        c = as_dict(s.get("conditions"))
        above, intact = as_dict(c.get("rp_above")).get("fingerprint"), as_dict(c.get("rp_d0")).get("fingerprint")
        mism += int(above is None or above != intact)
    res["I4"] = dict(what="noise above the threshold leaves the population identical to the intact compass",
                     mismatches=mism, verdict="ok" if mism == 0 else "bug")
    problems = list(extra_problems) + structure_problems(seeds)
    res["I5"] = dict(what="every population is finite and every record is present and complete", problems=len(problems),
                     examples=problems[:20], verdict="ok" if not problems else "bug")

    # ---------------- expectations from physics
    aniso, fit = physics.get("relative_anisotropy"), physics.get("fit_oos_rms_over_sd")
    res["Q2"] = dict(what="the compass works: relative anisotropy of the singlet yield >= 1%", relative_anisotropy=aniso,
                     verdict=verdict(is_num(aniso) and aniso >= TOL["anisotropy_min"]))
    res["Q3"] = dict(what="the spherical-harmonic fit reproduces the yield at 400 directions outside the fit grid "
                          "within 1% of its SD", fit_oos_rms_over_sd=fit, verdict=verdict(within(fit, TOL["fit_oos_max"])))
    if problems:
        for k in ("Q1", "H1", "H2", "H3", "H4"):
            res[k] = dict(verdict=NOT_TESTABLE, not_testable_because=["records incomplete or not finite (I5)"])
        return res

    # ---------------- atlas
    acc = {}
    for key in ATLAS_KEYS:
        labels = [s["atlas"][key]["shape"] for s in seeds]
        acc[key] = dict(correct=int(sum(lab == ATLAS_TARGET[key] for lab in labels)), labels=counts(labels))
    gate = {key: acc[key]["correct"] >= need for key in ATLAS_KEYS}
    res["Q1"] = dict(what="the pipeline recognises every atlas shape in at least 8 of 10 seeds",
                     atlas={k: acc[k] for k in SHAPES}, seeds_needed=need, verdict=verdict(all(gate[k] for k in SHAPES)))

    # ---------------- hypotheses
    def shape_hyp(cond, target, gates):
        labels = [s["conditions"][cond]["shape"] for s in seeds]
        hits = int(sum(lab == target for lab in labels))
        why = [f"atlas {g}: recognised in {acc[g]['correct']} of {n} seeds" for g in gates if not gate[g]]
        return dict(condition=cond, target=target, hits=hits, seeds_needed=need, labels=counts(labels), gates=list(gates),
                    not_testable_because=why, verdict=NOT_TESTABLE if why else verdict(hits >= need))

    res["H1"] = shape_hyp("rp_d0", "rp2", ("rp2",))
    res["H2"] = shape_hyp("polar", "sphere", ("sphere",))
    h3a, h3b = shape_hyp("yaw_inc66", "ring", ("ring_yaw_inc66",)), shape_hyp("yaw_inc0", "ring", ("ring_yaw_inc0",))
    rings = {g: dict(**acc[g], median_snr=[s["atlas"][g].get("median_snr") for s in seeds]) for g in GATE_RINGS}
    # H3 is the conjunction of its two clauses: one testable clause that fails decides it; otherwise a closed gate
    # leaves it not testable
    clause_verdicts = (h3a["verdict"], h3b["verdict"])
    v3 = NOT_MET if NOT_MET in clause_verdicts else NOT_TESTABLE if NOT_TESTABLE in clause_verdicts else MET
    why3 = h3a["not_testable_because"] + h3b["not_testable_because"] if v3 == NOT_TESTABLE else []
    res["H3"] = dict(inclination_66=h3a, equator=h3b, gate_rings=rings, not_testable_because=why3, verdict=v3)
    c0 = np.array([s["conditions"]["yaw_inc0"]["confusion_180"] for s in seeds], float)
    c66 = np.array([s["conditions"]["yaw_inc66"]["confusion_180"] for s in seeds], float)
    clauses = dict(equator_ci95_low_at_least_0_40=bool(t_ci(c0)[0] >= CONF_EQUATOR_MIN),
                   inc66_ci95_high_at_most_0_05=bool(t_ci(c66)[1] <= CONF_66_MAX))
    res["H4"] = dict(confusion_equator=dict(mean=float(np.mean(c0)), ci=t_ci(c0), per_seed=c0.tolist()),
                     confusion_66=dict(mean=float(np.mean(c66)), ci=t_ci(c66), per_seed=c66.tolist()),
                     clauses=clauses, verdict=verdict(all(clauses.values())))

    # ---------------- exploratory
    res["E1_rf_dose"] = {f"d{d:g}": counts([s["conditions"][f"rp_d{d:g}"]["shape"] for s in seeds]) for d in D_LEVELS}
    res["E2_design"] = {k: counts([s["design"][k]["shape"] for s in seeds]) for k in DESIGN_KEYS}

    def margins(rec_of):
        out = {}
        for z in ("Z2", "Z3"):
            for h in ("H1", "H2"):
                vals = []
                for s in seeds:
                    rec = rec_of(s)
                    bars, tau = rec["bars"][z][h], rec["thresholds"][z][h]
                    vals.append(float(bars[0] / tau) if bars and tau > 0 else float("nan"))
                out[f"{z}_{h}"] = float(np.nanmedian(vals)) if np.any(np.isfinite(vals)) else float("nan")
        return out

    res["E3_bar_margins"] = {c: margins(lambda s, c=c: s["conditions"][c]) for c in ("rp_d0", "polar") + YAW}
    res["E3_bar_margins"].update({f"atlas_{k}": margins(lambda s, k=k: s["atlas"][k]) for k in ATLAS_KEYS})
    return res


# --------------------------------------------------------------------------------------------------------------
# Run directory: seed files, provenance, validity
# --------------------------------------------------------------------------------------------------------------
def load(run_dir):
    seeds, problems = [], []
    for f in sorted(glob.glob(os.path.join(run_dir, "seed_*.json"))):
        name = os.path.basename(f)
        m = re.fullmatch(r"seed_(\d+)\.json", name)
        try:
            with open(f) as fh:
                s = json.load(fh)
        except (OSError, ValueError):
            problems.append(f"{name}: unreadable")
            continue
        if not isinstance(s, dict):
            problems.append(f"{name}: not a seed record")
            continue
        if m is None or s.get("seed") != int(m.group(1)):
            problems.append(f"{name}: seed field {s.get('seed')!r} differs from the file name")
        seeds.append(s)
    seeds.sort(key=lambda s: s.get("seed") if isinstance(s.get("seed"), int) else -1)
    return seeds, problems


def provenance(run_dir, code_dir, got):
    problems = []
    try:
        with open(os.path.join(run_dir, "meta_run.json")) as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        return dict(ok=False, problems=["meta_run.json missing or unreadable"])
    meta = as_dict(meta)
    if meta.get("mode") != "run":
        problems.append("meta_run.json does not come from run mode")
    for key, name in (("code_sha256", "mag_s2.py"), ("analysis_sha256", "analyze.py"), ("physics_sha256", "physics.json")):
        path = os.path.join(code_dir, name)
        if not os.path.exists(path) or meta.get(key) != sha256(path):
            problems.append(f"{key} of the run differs from {name} next to analyze.py")
    params = as_dict(meta.get("params"))
    if (params.get("rel_factor"), params.get("null_factor")) != FROZEN_FACTORS:
        problems.append("threshold factors differ from the calibrated 0.4 and 4")
    if meta.get("blas_threads") != [1]:
        problems.append("BLAS did not run on one thread")
    if sorted(meta.get("seeds") or []) != got:
        problems.append("the seed files differ from the seeds the run was started with")
    return dict(ok=not problems, problems=problems)


def evaluate(run_dir, code_dir):
    seeds, load_problems = load(run_dir)
    with open(os.path.join(code_dir, "physics.json")) as fh:
        physics = json.load(fh)
    got = sorted(s["seed"] for s in seeds if isinstance(s.get("seed"), int))
    res = analyze(seeds, physics, load_problems)
    res["seeds"] = got
    res["n_seeds"] = len(seeds)
    res["physics_sha256"] = sha256(os.path.join(code_dir, "physics.json"))
    res["provenance"] = provenance(run_dir, code_dir, got)
    res["run_valid"] = bool(all(res[i]["verdict"] == "ok" for i in INVARIANTS) and res["provenance"]["ok"])
    res["confirmatory"] = bool(res["run_valid"] and got == CONFIRMATORY_SEEDS and len(seeds) == len(CONFIRMATORY_SEEDS))
    return res


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: python3 -I analyze.py RUN_DIR")
    run_dir = sys.argv[1]
    here = os.path.dirname(os.path.abspath(__file__))
    res = evaluate(run_dir, here)
    with open(os.path.join(run_dir, "results.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    if not res["run_valid"]:
        bad = [i for i in INVARIANTS if res[i]["verdict"] != "ok"] + res["provenance"]["problems"]
        print("RUN INVALID:", "; ".join(bad))
    elif not res["confirmatory"]:
        print(f"NOT CONFIRMATORY: seeds {res['seeds']} differ from 1-10; the verdicts below check the pipeline only")
    for h in ORDER:
        why = res[h].get("not_testable_because")
        print(h, res[h]["verdict"], ("| " + "; ".join(why)) if why else "", "" if res["run_valid"] else "(provisional)")


if __name__ == "__main__":
    main()
