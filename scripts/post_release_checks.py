"""Checks added in v0.1.1 after an external review of v0.1.0. They are exploratory and leave the frozen protocol, code,
run and verdicts unchanged.

A. Readout rank. Symmetry under field reversal makes every neuron an even function of the field direction, so the
   population factors through RP². The population is RP² itself when, in addition, the readout tells apart all
   directions except opposite ones. The yield is a quadratic form of the direction up to 3e-5 of its anisotropic
   power, and quadratic forms without their trace span five dimensions (the Veronese map, an embedding of RP²). So for
   the dominant part the condition is that the map from these five dimensions to the neuron inputs has rank 5. For
   the polar control the receptor signals are linear in the direction, and the condition is rank 3. The script
   reports the singular values of the noiseless neuron inputs over 2000 directions, and the smallest distance between
   population states of two directions whose axes differ by at least 10 degrees (for the polar control, whose
   directions differ by at least 10 degrees), in units of the median distance. A single axial signal such as n_z^2
   is shown for contrast: its states form a segment, and the smallest distance falls to zero.
B. Radio-frequency dose without neural noise. The noiseless population at d = 0.75 and d = 0.9 goes through the frozen
   classifier. For d < 1 the compression keeps the noiseless map invertible, so the E1 labels at d = 0.75 measure
   what the pipeline can detect at fixed neural noise.
C. Exact bounds for the heading errors at inclination 66 degrees. The frozen t-interval over ten zero shares has
   zero width. The script gives Clopper-Pearson upper bounds for two different quantities. For the error probability
   of a single test heading it gives the bound from the 300 headings of one bird and the bound from all 3000 headings
   pooled; the second assumes that headings are independent across birds. For the share of birds that would show at
   least one error among their 300 test headings it gives the bound from 0 of 10 birds. A bound on the error
   probability per heading that allows for differences between birds needs a model of those differences.

Usage: python3 -I scripts/post_release_checks.py   (about 6 minutes; writes results/checks/post_release_checks.json)"""
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "code"))
import mag_s2 as m  # noqa: E402
import numpy as np  # noqa: E402

SEEDS = range(1, 11)
P = m.Params()
with open(os.path.join(ROOT, "code", "physics.json")) as fh:
    PHYS = json.load(fh)
PATTERN = m.Pattern(PHYS, P)


def zscore(X):
    return (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)


def separation(states, dirs, axial, min_deg=10.0):
    """Smallest distance between two states whose directions (axes, when axial) differ by at least min_deg,
    in units of the median distance over all such pairs."""
    Z = zscore(states)
    d2 = np.maximum((Z ** 2).sum(1)[:, None] + (Z ** 2).sum(1)[None, :] - 2 * Z @ Z.T, 0.0)
    c = dirs @ dirs.T
    ang = np.degrees(np.arccos(np.clip(np.abs(c) if axial else c, -1.0, 1.0)))
    mask = np.triu(ang >= min_deg, 1)
    dist = np.sqrt(d2[mask])
    return float(dist.min() / np.median(dist))


def readout_rank():
    dirs = m.fibonacci_sphere(2000)
    out = {}
    for s in SEEDS:
        v, R, eyes = m.retina(P, m.rng_for(s, 1))
        W, sign, bias = m.neurons(P, m.rng_for(s, 2), v, eyes)
        row = {}
        for kind, rank in (("rp", 5), ("polar", 3)):
            z = m.receptor_signals(kind, dirs, R, PATTERN)
            U = sign[:, None] * (W @ z)
            sv = np.linalg.svd(U - U.mean(1, keepdims=True), compute_uv=False)
            rates = m.rates(z, W, sign, bias, P)
            row[kind] = dict(expected_rank=rank, singular_values_over_largest=[float(x) for x in sv[:rank + 2] / sv[0]],
                             smallest_kept_over_largest=float(sv[rank - 1] / sv[0]),
                             next_over_smallest_kept=float(sv[rank] / sv[rank - 1]),
                             min_separation_over_median=separation(rates, dirs, axial=(kind == "rp")))
        out[f"seed_{s}"] = row
    nz2 = (dirs[:, 2] ** 2)[:, None] * np.ones((1, 3))       # contrast: one axial signal copied into three neurons
    out["contrast_single_axial_signal_nz2"] = dict(min_separation_over_median=separation(nz2, dirs, axial=True))
    return out


def noiseless_dose():
    out = {}
    for s in SEEDS:
        v, R, eyes = m.retina(P, m.rng_for(s, 1))
        W, sign, bias = m.neurons(P, m.rng_for(s, 2), v, eyes)
        dirs = m.uniform_dirs(m.rng_for(s, 3), P.n_dirs)
        z = m.receptor_signals("rp", dirs, R, PATTERN)
        row = {}
        for d in (0.75, 0.9):
            X = m.rates((1.0 - d) * z, W, sign, bias, P)
            rec = m.classify(X, P, m.rng_for(s, 30, int(round(100 * d))))
            row[f"d{d:g}"] = dict(shape=rec["shape"], signature=rec["signature"],
                                  orientability_ratio=rec["orientability_ratio"])
        out[f"seed_{s}"] = row
    return out


def clopper_pearson_upper_for_zero(n, alpha=0.05):
    """Two-sided exact interval for 0 events in n trials is [0, 1 - (alpha/2)^(1/n)]."""
    return float(1.0 - (alpha / 2.0) ** (1.0 / n))


def heading_bounds():
    with open(os.path.join(ROOT, "results", "confirmatory_run_v01", "results.json")) as fh:
        res = json.load(fh)
    shares = res["H4"]["confusion_66"]["per_seed"]
    n_test = P.n_dirs - P.n_dirs // 2                      # the decoder tests on the second half of the headings
    assert all(x == 0.0 for x in shares)
    return dict(per_seed_share=shares, test_headings_per_bird=n_test, birds=len(shares),
                frozen_t_interval=res["H4"]["confusion_66"]["ci"],
                per_heading_error_upper_95_one_bird=clopper_pearson_upper_for_zero(n_test),
                per_heading_error_upper_95_pooled_independent=clopper_pearson_upper_for_zero(n_test * len(shares)),
                share_of_birds_with_any_error_upper_95=clopper_pearson_upper_for_zero(len(shares)))


def main():
    t0 = time.time()
    out = dict(readout_rank=readout_rank(), noiseless_dose=noiseless_dose(), heading_errors_66=heading_bounds())
    out["seconds"] = round(time.time() - t0, 1)
    path = os.path.join(ROOT, "results", "checks", "post_release_checks.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    rr = out["readout_rank"]
    for kind in ("rp", "polar"):
        lo = min(rr[f"seed_{s}"][kind]["smallest_kept_over_largest"] for s in SEEDS)
        nx = max(rr[f"seed_{s}"][kind]["next_over_smallest_kept"] for s in SEEDS)
        sep = min(rr[f"seed_{s}"][kind]["min_separation_over_median"] for s in SEEDS)
        print(f"A {kind}: rank {rr['seed_1'][kind]['expected_rank']} in every bird; smallest kept singular value "
              f">= {lo:.3f} of the largest, next one <= {nx:.2e} of it; smallest separation >= {sep:.3f} of the median")
    print(f"A contrast n_z^2: smallest separation {rr['contrast_single_axial_signal_nz2']['min_separation_over_median']:.2e}"
          " of the median")
    for d in ("d0.75", "d0.9"):
        labels = [out["noiseless_dose"][f"seed_{s}"][d]["shape"] for s in SEEDS]
        print(f"B noiseless {d}: " + ", ".join(f"{k} {labels.count(k)}" for k in sorted(set(labels))))
    hb = out["heading_errors_66"]
    print(f"C 66 degrees: 0 errors in {hb['birds']} x {hb['test_headings_per_bird']} test headings; frozen t-interval "
          f"{hb['frozen_t_interval']}; exact 95% upper bounds for the error probability per heading: one bird "
          f"{hb['per_heading_error_upper_95_one_bird']:.4f}, pooled assuming independence "
          f"{hb['per_heading_error_upper_95_pooled_independent']:.5f}; for the share of birds with any error "
          f"{hb['share_of_birds_with_any_error_upper_95']:.3f}")
    print(f"written {os.path.relpath(path, ROOT)} in {out['seconds']} s")


if __name__ == "__main__":
    main()
