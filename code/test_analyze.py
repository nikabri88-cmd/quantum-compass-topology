"""MAG-S2: synthetic tests of the analysis logic on fabricated records with known answers.
No model output is read. Each world changes one thing against a baseline in which every criterion is met, and most
thresholds are probed just inside and just outside. The last block builds run folders on disk and tests evaluate()
and main(): provenance, seed sets and the validity flags.
Exits with status 1 when any check fails. Usage: python3 -I test_analyze.py (from the folder that holds analyze.py)"""
import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyze as A  # noqa: E402

SIG = {"point": (0, 0, 0, 0), "ring": (1, 0, 1, 0), "sphere": (0, 1, 0, 1), "torus": (2, 1, 2, 1),
       "rp2": (1, 1, 0, 0), "klein": (2, 1, 1, 0), "other": (3, 0, 0, 0)}
# the record layout, written out here instead of being taken from analyze.py
SHAPES6 = ("point", "ring", "sphere", "torus", "rp2", "klein")
EXPECTED_CONDITIONS = ("rp_d0", "rp_d0.25", "rp_d0.5", "rp_d0.75", "rp_d1", "rp_above", "polar", "yaw_inc66", "yaw_inc0")
GATES = ("ring_yaw_inc66", "ring_yaw_inc0")
EXPECTED_ATLAS = {**{s: s for s in SHAPES6}, "ring_yaw_inc66": "ring", "ring_yaw_inc0": "ring"}
EXPECTED_DESIGN = tuple(f"N{n}_noise{z}_{k}" for n in (25, 50, 100) for z in ("0.25", "0.5", "1") for k in ("rp", "polar"))
fails, done = [], []


def check(name, cond):
    done.append(name)
    if not cond:
        fails.append(name)
    print(("ok   " if cond else "FAIL ") + name)


def topo(shape, fp="a"):
    b = SIG[shape]
    bars = {"Z2": {"H1": [2.0] * b[0] + [0.1], "H2": [2.0] * b[1] + [0.1]},
            "Z3": {"H1": [2.0] * b[2] + [0.1], "H2": [2.0] * b[3] + [0.1]}}
    thr = {z: {"H1": 0.5, "H2": 0.5} for z in ("Z2", "Z3")}
    return dict(n_points=540, bars=bars, null=copy.deepcopy(thr), thresholds=thr,
                betti={"Z2": {"H1": b[0], "H2": b[1]}, "Z3": {"H1": b[2], "H2": b[3]}},
                signature=list(b), shape=shape, fingerprint=fp, finite=True)


BASE = {"rp_d0": "rp2", "rp_d0.25": "rp2", "rp_d0.5": "rp2", "rp_d0.75": "point", "rp_d1": "point",
        "polar": "sphere", "yaw_inc66": "ring", "yaw_inc0": "ring"}


def world(seed, labels=None, atlas=None, conf0=0.5, conf66=0.0):
    labels = {**BASE, **(labels or {})}
    cond = {name: topo(lab, fp=f"fp{seed}{name}") for name, lab in labels.items()}
    cond["rp_above"] = dict(fingerprint=cond["rp_d0"]["fingerprint"], finite=True)
    cond["yaw_inc0"]["confusion_180"] = conf0
    cond["yaw_inc66"]["confusion_180"] = conf66
    at = {key: topo((atlas or {}).get(key, target)) for key, target in EXPECTED_ATLAS.items()}
    for g in GATES:
        at[g]["median_snr"] = 1.4
    design = {k: {"shape": "rp2" if k.endswith("_rp") else "sphere"} for k in EXPECTED_DESIGN}
    return dict(seed=seed, noiseless_antipodal_max_diff=1e-15, conditions=cond, atlas=at, design=design)


PHYS = dict(antipodal_max_diff=1e-16, odd_power_fraction=1e-29, radicalpy_max_diff=2e-8, resolvent_max_diff=1e-14,
            zeeman_cartesian_max_rel_diff=5e-16, relative_anisotropy=0.2, fit_rms_over_sd=0.003, fit_oos_rms_over_sd=0.003)


def seeds_of(per_seed=None, n=10, edit=None):
    seeds = [world(i, **(per_seed(i) if per_seed else {})) for i in range(1, n + 1)]
    if edit:
        edit(seeds)
    return seeds


def run(per_seed=None, physics=None, edit=None, n=10):
    return A.analyze(seeds_of(per_seed, n, edit), {**PHYS, **(physics or {})})


def first(k):
    """per-seed edit for the first k seeds"""
    return lambda kw: (lambda i: kw if i <= k else {})


def only(r, changed, expect):
    others = all(r[h]["verdict"] == A.MET for h in ("H1", "H2", "H3", "H4") if h not in changed)
    return all(r[c]["verdict"] == expect for c in changed) and others


def inv(r, name):
    return r[name]["verdict"]


# ---------------------------------------------------------------- record layout: analysis, model and this file agree
check("analysis expects the conditions, atlas keys, gate targets and design keys written out here",
      tuple(A.CONDITIONS) == EXPECTED_CONDITIONS and dict(A.ATLAS_TARGET) == EXPECTED_ATLAS
      and tuple(A.DESIGN_KEYS) == EXPECTED_DESIGN and tuple(A.GATE_RINGS) == GATES)
import mag_s2 as M  # noqa: E402

with open(os.path.join(HERE, "physics.json")) as _fh:
    _phys = json.load(_fh)
_P = M.replace(M.Params(), n_dirs=12, n_neurons=6, n_receptors=20)
_pops, _psi, _antip, _frac = M.compass_populations(_P, M.Pattern(_phys, _P), 9308)
check("the model writes exactly these conditions and yaw noise shares",
      tuple(_pops) == EXPECTED_CONDITIONS and tuple(_frac) == ("yaw_inc66", "yaw_inc0"))
check("model and analysis agree on shapes, yaw names, dose levels, design grid and calibrated factors",
      M.SHAPES == SHAPES6 == A.SHAPES and M.YAW_CONDITIONS == A.YAW and tuple(M.D_LEVELS) == tuple(A.D_LEVELS)
      and tuple(f"N{n}_noise{z:g}_{k}" for n in M.DESIGN_NEURONS for z in M.DESIGN_NOISE for k in ("rp", "polar"))
      == EXPECTED_DESIGN and (M.Params().rel_factor, M.Params().null_factor) == A.FROZEN_FACTORS == (0.4, 4.0))

# ---------------------------------------------------------------- baseline
r = run()
for h in A.ORDER:
    print(" ", h, r[h]["verdict"])
check("baseline: invariants ok", all(r[i]["verdict"] == "ok" for i in A.INVARIANTS))
check("baseline: expectations met", all(r[q]["verdict"] == A.MET for q in ("Q1", "Q2", "Q3")))
check("baseline: hypotheses met", all(r[h]["verdict"] == A.MET for h in ("H1", "H2", "H3", "H4")))
check("needed(): 8 of 10, 8 of 9, 4 of 5, 2 of 2", [A.needed(n) for n in (10, 9, 5, 2)] == [8, 8, 4, 2])

# ---------------------------------------------------------------- H1, H2: counts and gates at the boundary
check("H1: projective plane in 8 of 10 -> met", only(run(first(2)(dict(labels={"rp_d0": "sphere"}))), [], A.MET))
check("H1: projective plane in 7 of 10 -> not met",
      only(run(first(3)(dict(labels={"rp_d0": "sphere"}))), ["H1"], A.NOT_MET))
check("H1: 7 of 9 seeds -> not met (8 needed)",
      run(first(2)(dict(labels={"rp_d0": "sphere"})), n=9)["H1"]["verdict"] == A.NOT_MET)
check("H1: 4 of 5 seeds -> met (4 needed)", run(first(1)(dict(labels={"rp_d0": "sphere"})), n=5)["H1"]["verdict"] == A.MET)
check("H1: atlas RP2 recognised in 8 of 10 -> gate open, met",
      only(run(first(2)(dict(atlas={"rp2": "other"}))), [], A.MET))
r = run(first(3)(dict(atlas={"rp2": "other"})))
check("H1: atlas RP2 recognised in 7 of 10 -> not testable, with a reason",
      only(r, ["H1"], A.NOT_TESTABLE) and len(r["H1"]["not_testable_because"]) == 1)
check("H1 and H2 are gated only by their own shape: atlas ring and torus in 7 of 10 -> H1, H2 met",
      all(run(first(3)(dict(atlas={"ring": "other", "torus": "other"})))[h]["verdict"] == A.MET for h in ("H1", "H2")))
check("H1 reads the intact compass: RF dose d = 0.25 in 0 of 10 -> H1 met",
      run(lambda i: dict(labels={"rp_d0.25": "point"}))["H1"]["verdict"] == A.MET)
check("H2: sphere in 8 of 10 -> met", only(run(first(2)(dict(labels={"polar": "rp2"}))), [], A.MET))
check("H2: sphere in 7 of 10 -> not met", only(run(first(3)(dict(labels={"polar": "rp2"}))), ["H2"], A.NOT_MET))
check("H2: atlas sphere recognised in 8 of 10 -> met", only(run(first(2)(dict(atlas={"sphere": "other"}))), [], A.MET))
check("H2: atlas sphere recognised in 7 of 10 -> not testable",
      only(run(first(3)(dict(atlas={"sphere": "other"}))), ["H2"], A.NOT_TESTABLE))

# ---------------------------------------------------------------- H3: counts and the noise-matched gate rings
check("H3: ring in 8 of 10 at both inclinations -> met",
      only(run(first(2)(dict(labels={"yaw_inc66": "point", "yaw_inc0": "other"}))), [], A.MET))
check("H3: equator ring in 7 of 10 -> not met", only(run(first(3)(dict(labels={"yaw_inc0": "other"}))), ["H3"], A.NOT_MET))
check("H3: inclination-66 ring in 7 of 10 -> not met",
      only(run(first(3)(dict(labels={"yaw_inc66": "point"}))), ["H3"], A.NOT_MET))
check("H3: gate ring of inclination 66 in 8 of 10 -> met",
      only(run(first(2)(dict(atlas={"ring_yaw_inc66": "point"}))), [], A.MET))
r = run(first(3)(dict(atlas={"ring_yaw_inc66": "point"})))
check("H3: gate ring of inclination 66 in 7 of 10 -> not testable, with a reason",
      only(r, ["H3"], A.NOT_TESTABLE) and len(r["H3"]["not_testable_because"]) == 1)
check("H3: gate ring of the equator in 7 of 10 -> not testable",
      only(run(first(3)(dict(atlas={"ring_yaw_inc0": "other"}))), ["H3"], A.NOT_TESTABLE))
r = run(first(3)(dict(atlas={"ring": "other", "point": "ring"})))
check("H3 is gated by the noise-matched rings: atlas ring and point in 7 of 10 -> H3 met, Q1 not met",
      r["H3"]["verdict"] == A.MET and r["Q1"]["verdict"] == A.NOT_MET)
check("H3 reports the gate rings and their median SNR",
      set(run()["H3"]["gate_rings"]) == set(GATES) and run()["H3"]["gate_rings"]["ring_yaw_inc0"]["median_snr"][0] == 1.4)
r = run(first(3)(dict(atlas={"ring_yaw_inc66": "point"}, labels={"yaw_inc0": "other"})))
check("H3: inclination-66 gate closed and the equator clause testable but failed -> not met",
      only(r, ["H3"], A.NOT_MET) and r["H3"]["not_testable_because"] == [])
r = run(first(3)(dict(atlas={"ring_yaw_inc66": "point"}, labels={"yaw_inc66": "point"})))
check("H3: each clause has its own gate (inclination-66 gate closed, that clause short, equator fine) -> not testable",
      only(r, ["H3"], A.NOT_TESTABLE))
r = run(first(3)(dict(atlas={"ring_yaw_inc0": "other"})))
check("Q1 counts only the six atlas shapes: a gate ring in 7 of 10 -> Q1 met, H3 not testable",
      r["Q1"]["verdict"] == A.MET and r["H3"]["verdict"] == A.NOT_TESTABLE)

# ---------------------------------------------------------------- H4: cut-offs and the interval
check("H4: equator 0.401 in every seed -> met", run(lambda i: dict(conf0=0.401))["H4"]["verdict"] == A.MET)
check("H4: equator 0.399 in every seed -> not met", only(run(lambda i: dict(conf0=0.399)), ["H4"], A.NOT_MET))
check("H4: inclination 66 at 0.049 in every seed -> met", run(lambda i: dict(conf66=0.049))["H4"]["verdict"] == A.MET)
check("H4: inclination 66 at 0.051 in every seed -> not met", only(run(lambda i: dict(conf66=0.051)), ["H4"], A.NOT_MET))
alt = lambda a, d: (lambda i: a + d if i % 2 else a - d)  # noqa: E731
r = run(lambda i: dict(conf0=alt(0.51, 0.15)(i)))
check("H4: equator t-interval low 0.397 (normal, 90% or ddof 0 would pass) -> not met", r["H4"]["verdict"] == A.NOT_MET)
r = run(lambda i: dict(conf0=alt(0.53, 0.15)(i)))
check("H4: equator t-interval low 0.417 (a 99% interval would fail) -> met", r["H4"]["verdict"] == A.MET)
r = run(lambda i: dict(conf66=alt(0.0288, 0.0288)(i)))
check("H4: inclination-66 t-interval high 0.0505 (normal, 90% or ddof 0 would pass) -> not met",
      r["H4"]["verdict"] == A.NOT_MET)
r = run(lambda i: dict(conf66=alt(0.026, 0.026)(i)))
check("H4: inclination-66 t-interval high 0.0456 (a 99% interval would fail) -> met", r["H4"]["verdict"] == A.MET)
check("H4: equator mean 0.40 with a wide interval -> not met",
      run(lambda i: dict(conf0=alt(0.40, 0.15)(i)))["H4"]["verdict"] == A.NOT_MET)
check("H4: inclination-66 mean 0.04 with a wide interval -> not met",
      run(lambda i: dict(conf66=alt(0.04, 0.04)(i)))["H4"]["verdict"] == A.NOT_MET)
check("H4: equator 0.475, inclination 66 at 0.02 -> met",
      run(lambda i: dict(conf0=0.5 - 0.05 * (i % 2), conf66=0.02))["H4"]["verdict"] == A.MET)

# ---------------------------------------------------------------- expectations
for shape in SHAPES6:
    r = run(first(3)(dict(atlas={shape: "other"})))
    check(f"Q1: atlas {shape} in 7 of 10 -> Q1 not met", r["Q1"]["verdict"] == A.NOT_MET)
check("Q1: atlas Klein bottle in 8 of 10 -> met", run(first(2)(dict(atlas={"klein": "other"})))["Q1"]["verdict"] == A.MET)
check("Q1: atlas torus in 7 of 10 -> hypotheses still met", only(run(first(3)(dict(atlas={"torus": "other"}))), [], A.MET))
check("Q2: anisotropy 1.01% -> met", run(physics=dict(relative_anisotropy=0.0101))["Q2"]["verdict"] == A.MET)
check("Q2: anisotropy 0.99% -> not met", run(physics=dict(relative_anisotropy=0.0099))["Q2"]["verdict"] == A.NOT_MET)
check("Q2: anisotropy missing -> not met", run(physics=dict(relative_anisotropy=None))["Q2"]["verdict"] == A.NOT_MET)
check("Q3: out-of-sample residual 0.99% of SD -> met", run(physics=dict(fit_oos_rms_over_sd=0.0099))["Q3"]["verdict"] == A.MET)
check("Q3: out-of-sample residual 1.01% of SD -> not met",
      run(physics=dict(fit_oos_rms_over_sd=0.0101))["Q3"]["verdict"] == A.NOT_MET)
check("Q3 reads the out-of-sample residual: in-sample 5%, out-of-sample 0.3% -> met",
      run(physics=dict(fit_rms_over_sd=0.05))["Q3"]["verdict"] == A.MET)
check("Q3: out-of-sample residual missing -> not met", run(physics=dict(fit_oos_rms_over_sd=None))["Q3"]["verdict"] == A.NOT_MET)

# ---------------------------------------------------------------- invariants I1-I4 at the boundary
check("I1: physics antipodal 0.9e-12 -> ok", inv(run(physics=dict(antipodal_max_diff=0.9e-12)), "I1") == "ok")
check("I1: physics antipodal 2e-12 -> bug", inv(run(physics=dict(antipodal_max_diff=2e-12)), "I1") == "bug")
check("I1: physics antipodal NaN -> bug", inv(run(physics=dict(antipodal_max_diff=float("nan"))), "I1") == "bug")


def set_antip(value, idx=4):
    def edit(seeds):
        seeds[idx]["noiseless_antipodal_max_diff"] = value
    return edit


check("I1: population antipodal 0.9e-9 in seed 5 -> ok", inv(run(edit=set_antip(0.9e-9)), "I1") == "ok")
check("I1: population antipodal 2e-9 in seed 5 -> bug", inv(run(edit=set_antip(2e-9)), "I1") == "bug")
r = run(edit=set_antip(float("nan")))
check("I1 and I5: population antipodal NaN in seed 5 (not the first) -> bug", inv(r, "I1") == "bug" and inv(r, "I5") == "bug")
check("I2: odd power 0.9e-10 -> ok", inv(run(physics=dict(odd_power_fraction=0.9e-10)), "I2") == "ok")
check("I2: odd power 2e-10 -> bug", inv(run(physics=dict(odd_power_fraction=2e-10)), "I2") == "bug")
check("I2: odd power missing -> bug", inv(run(physics=dict(odd_power_fraction=None)), "I2") == "bug")
check("I3: time-domain difference 0.9e-6 -> ok", inv(run(physics=dict(radicalpy_max_diff=0.9e-6)), "I3") == "ok")
check("I3: time-domain difference 2e-6 -> bug", inv(run(physics=dict(radicalpy_max_diff=2e-6)), "I3") == "bug")
check("I3: time-domain difference missing -> bug", inv(run(physics=dict(radicalpy_max_diff=None)), "I3") == "bug")
check("I3: resolvent difference 2e-10 -> bug", inv(run(physics=dict(resolvent_max_diff=2e-10)), "I3") == "bug")
check("I3: Cartesian Zeeman difference 2e-12 -> bug", inv(run(physics=dict(zeeman_cartesian_max_rel_diff=2e-12)), "I3") == "bug")


def edit_cond(seed_idx, cond, key, value, where="conditions"):
    def edit(seeds):
        if value is KeyError:
            del seeds[seed_idx][where][cond][key]
        else:
            seeds[seed_idx][where][cond][key] = value
    return edit


check("I4: above-threshold fingerprint differs -> bug",
      inv(run(edit=edit_cond(2, "rp_above", "fingerprint", "different")), "I4") == "bug")
check("I4: above-threshold fingerprint missing -> bug", inv(run(edit=edit_cond(2, "rp_above", "fingerprint", KeyError)), "I4") == "bug")

# ---------------------------------------------------------------- I5: every completeness sub-check


def i5_bug(edit):
    r = run(edit=edit)
    return inv(r, "I5") == "bug" and all(r[h]["verdict"] == A.NOT_TESTABLE for h in ("Q1", "H1", "H2", "H3", "H4"))


def drop(where, key, sub=None):
    def edit(seeds):
        if sub is None:
            del seeds[3][where][key]
        else:
            del seeds[3][where][key][sub]
    return edit


check("I5: non-finite condition population -> bug, hypotheses not testable",
      i5_bug(edit_cond(6, "polar", "finite", False)))
check("I5: 'finite' flag missing -> bug", i5_bug(edit_cond(6, "polar", "finite", KeyError)))
check("I5: non-finite above-threshold population -> bug", i5_bug(edit_cond(6, "rp_above", "finite", False)))
check("I5: non-finite atlas point population -> bug", i5_bug(edit_cond(1, "point", "finite", False, where="atlas")))
check("I5: non-finite gate ring population -> bug", i5_bug(edit_cond(1, "ring_yaw_inc0", "finite", False, where="atlas")))
check("I5: condition record without its shape -> bug", i5_bug(edit_cond(1, "yaw_inc66", "shape", KeyError)))
check("I5: atlas record without its Betti numbers -> bug", i5_bug(edit_cond(1, "klein", "betti", KeyError, where="atlas")))
check("I5: unknown label -> bug", i5_bug(edit_cond(1, "rp_d0", "shape", "donut")))
check("I5: a missing condition -> bug, no crash", i5_bug(drop("conditions", "yaw_inc0")))
check("I5: a missing atlas gate ring -> bug, no crash", i5_bug(drop("atlas", "ring_yaw_inc66")))
check("I5: a missing design point -> bug", i5_bug(drop("design", "N25_noise1_polar")))
check("I5: heading confusion NaN -> bug", i5_bug(edit_cond(5, "yaw_inc0", "confusion_180", float("nan"))))
check("I5: heading confusion 1.5 -> bug", i5_bug(edit_cond(5, "yaw_inc66", "confusion_180", 1.5)))
check("I5: heading confusion missing -> bug", i5_bug(edit_cond(5, "yaw_inc66", "confusion_180", KeyError)))


def nan_bar(seeds):
    seeds[8]["conditions"]["rp_d0"]["bars"]["Z3"]["H2"][0] = float("nan")


def nan_threshold(seeds):
    seeds[8]["atlas"]["sphere"]["thresholds"]["Z2"]["H2"] = float("inf")


check("I5: NaN bar -> bug", i5_bug(nan_bar))
check("I5: infinite threshold -> bug", i5_bug(nan_threshold))
check("I5: noiseless antipodal difference missing -> bug", i5_bug(lambda seeds: seeds[7].pop("noiseless_antipodal_max_diff")))
r = A.analyze([], PHYS)
check("I5: no seeds -> bug, hypotheses not testable",
      inv(r, "I5") == "bug" and all(r[h]["verdict"] == A.NOT_TESTABLE for h in ("H1", "H2", "H3", "H4")))
r = A.analyze(seeds_of(), PHYS, extra_problems=["seed_3.json: seed field 4 differs from the file name"])
check("I5: a load problem passed in -> bug", inv(r, "I5") == "bug")

# ---------------------------------------------------------------- exploratory outputs
r = run(first(4)(dict(labels={"rp_d0.75": "rp2"})))
check("E1: dose 0.75 counted per label", r["E1_rf_dose"]["d0.75"] == {"point": 6, "rp2": 4} and r["E1_rf_dose"]["d0"] == {"rp2": 10})
check("E2: all 18 design points counted", len(r["E2_design"]) == 18 and r["E2_design"]["N25_noise1_rp"] == {"rp2": 10})
check("E3: bar margin is the longest bar over its threshold (2.0 / 0.5)",
      r["E3_bar_margins"]["rp_d0"]["Z2_H1"] == 4.0 and r["E3_bar_margins"]["atlas_ring_yaw_inc0"]["Z2_H1"] == 4.0)

# ---------------------------------------------------------------- run folders: evaluate(), provenance and main()
TMP = tempfile.mkdtemp(prefix="mag_s2_test_")


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def make_run(tag, seeds=range(1, 11), meta=None, drop_meta=False, rename=None, edit=None, copy_analysis=False):
    """Code folder (stand-in files, or a copy of analyze.py) and run folder whose metadata match it."""
    code, rd = os.path.join(TMP, tag, "code"), os.path.join(TMP, tag, "run")
    os.makedirs(code)
    os.makedirs(rd)
    with open(os.path.join(code, "mag_s2.py"), "w") as fh:
        fh.write("# stand-in model code\n")
    if copy_analysis:
        shutil.copy(os.path.join(HERE, "analyze.py"), code)
    else:
        with open(os.path.join(code, "analyze.py"), "w") as fh:
            fh.write("# stand-in analysis\n")
    with open(os.path.join(code, "physics.json"), "w") as fh:
        json.dump(PHYS, fh)
    recs = [world(i) for i in seeds]
    if edit:
        edit(recs)
    for rec in recs:
        name = (rename or {}).get(rec["seed"], f"seed_{rec['seed']}.json")
        with open(os.path.join(rd, name), "w") as fh:
            json.dump(rec, fh)
    m = dict(mode="run", code_sha256=sha(os.path.join(code, "mag_s2.py")), analysis_sha256=sha(os.path.join(code, "analyze.py")),
             physics_sha256=sha(os.path.join(code, "physics.json")), params=dict(rel_factor=0.4, null_factor=4.0),
             blas_threads=[1], seeds=list(seeds))
    m.update(meta or {})
    if not drop_meta:
        with open(os.path.join(rd, "meta_run.json"), "w") as fh:
            json.dump(m, fh)
    return rd, code


def ev(tag, **kw):
    rd, code = make_run(tag, **kw)
    return A.evaluate(rd, code)


r = ev("valid")
check("run folder: seeds 1-10 with matching metadata -> valid and confirmatory", r["run_valid"] and r["confirmatory"])
r = ev("nine", seeds=range(1, 10))
check("run folder: seeds 1-9 -> valid, not confirmatory", r["run_valid"] and not r["confirmatory"])
r = ev("stray", seeds=list(range(1, 11)) + [9201])
check("run folder: seeds 1-10 plus 9201 -> not confirmatory", r["run_valid"] and not r["confirmatory"])
r = ev("pipe", seeds=[9201, 9202])
check("run folder: pipeline-test seeds -> valid, not confirmatory", r["run_valid"] and not r["confirmatory"])
for tag, kw, what in (("nometa", dict(drop_meta=True), "metadata missing"),
                      ("code", dict(meta=dict(code_sha256="0" * 64)), "code hash differs"),
                      ("analysis", dict(meta=dict(analysis_sha256="0" * 64)), "analysis hash differs"),
                      ("physics", dict(meta=dict(physics_sha256="0" * 64)), "physics hash differs"),
                      ("factors", dict(meta=dict(params=dict(rel_factor=0.35, null_factor=3.0))), "threshold factors differ"),
                      ("relonly", dict(meta=dict(params=dict(rel_factor=0.35, null_factor=4.0))), "the share factor alone differs"),
                      ("nullonly", dict(meta=dict(params=dict(rel_factor=0.4, null_factor=3.0))), "the null factor alone differs"),
                      ("blas", dict(meta=dict(blas_threads=[2])), "two BLAS threads"),
                      ("blasmix", dict(meta=dict(blas_threads=[1, 2])), "BLAS libraries on one and two threads"),
                      ("blasnone", dict(meta=dict(blas_threads=None)), "BLAS thread count unknown"),
                      ("mode", dict(meta=dict(mode="calibrate")), "metadata from another mode"),
                      ("started", dict(meta=dict(seeds=list(range(1, 12)))), "seed 11 started but its file is missing"),
                      ("renamed", dict(rename={3: "seed_33.json"}), "a seed file named after another seed"),
                      ("invariant", dict(edit=lambda recs: recs[0]["conditions"]["rp_above"].update(fingerprint="x")),
                       "an invariant fails"),
                      ("incomplete", dict(edit=lambda recs: recs[0]["conditions"].pop("polar")), "a record is incomplete")):
    r = ev(tag, **kw)
    check(f"run folder: {what} -> invalid, not confirmatory", not r["run_valid"] and not r["confirmatory"])
rd, code = make_run("unreadable")
with open(os.path.join(rd, "seed_5.json"), "w") as fh:
    fh.write("{ not json")
r = A.evaluate(rd, code)
check("run folder: an unreadable seed file -> invalid, I5 bug", not r["run_valid"] and r["I5"]["verdict"] == "bug")

rd, code = make_run("main", copy_analysis=True)
p = subprocess.run([sys.executable, "-I", os.path.join(code, "analyze.py"), rd], capture_output=True, text=True, timeout=120)
with open(os.path.join(rd, "results.json")) as fh:
    written = json.load(fh)
check("main(): writes results.json, valid and confirmatory, no warning",
      p.returncode == 0 and written["confirmatory"] and "INVALID" not in p.stdout and "NOT CONFIRMATORY" not in p.stdout)
check("main(): prints every verdict", all(any(line.startswith(h + " ") for line in p.stdout.splitlines()) for h in A.ORDER))
rd, code = make_run("main_pipe", seeds=[9201, 9202], copy_analysis=True)
p = subprocess.run([sys.executable, "-I", os.path.join(code, "analyze.py"), rd], capture_output=True, text=True, timeout=120)
check("main(): pipeline-test seeds print NOT CONFIRMATORY", p.returncode == 0 and "NOT CONFIRMATORY" in p.stdout)
rd, code = make_run("main_bad", meta=dict(code_sha256="0" * 64), copy_analysis=True)
p = subprocess.run([sys.executable, "-I", os.path.join(code, "analyze.py"), rd], capture_output=True, text=True, timeout=120)
check("main(): a hash mismatch prints RUN INVALID and marks the verdicts provisional",
      p.returncode == 0 and "RUN INVALID" in p.stdout and "(provisional)" in p.stdout)
shutil.rmtree(TMP, ignore_errors=True)

print(f"\nFAILED {len(fails)} of {len(done)}:" if fails else f"\nall {len(done)} synthetic checks passed", fails if fails else "")
sys.exit(1 if fails else 0)
