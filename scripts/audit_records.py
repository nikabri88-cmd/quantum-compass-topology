"""Audit of the seed records, added in v0.1.1 after an external review. The frozen analysis (code/analyze.py) is
unchanged; this script adds the two checks the review found missing from it.

1. Every number anywhere in every seed record, nested fields included, is finite.
2. Every label stored together with its persistence bars is derived again from the stored bars and shuffle maxima
   with the frozen decision rule (mag_s2.decide with the frozen factors 0.4 and 4). Thresholds, Betti numbers,
   orientability ratio, signature and shape must all match the stored values.
3. Every label of the requirements map, stored without bars, agrees with its stored signature and orientability
   ratio under the frozen signature table.

Usage: python3 -I scripts/audit_records.py RUN_DIR
Prints every problem found and exits with status 1 when there is any."""
import glob
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "code"))
import mag_s2  # noqa: E402

REL, NULL = 0.4, 4.0   # frozen factors, as checked by analyze.py against the run metadata


def walk_numbers(obj, path=""):
    """Yield (path, value) for every number in a nested JSON structure."""
    if isinstance(obj, bool) or obj is None or isinstance(obj, str):
        return
    if isinstance(obj, (int, float)):
        yield path, obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk_numbers(v, f"{path}.{k}" if path else str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from walk_numbers(v, f"{path}[{i}]")


def same(a, b):
    if isinstance(a, float) or isinstance(b, float):
        return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-15)
    return a == b


def check_classified(rec, where, problems):
    """Derive the label of a record that carries bars and shuffle maxima and compare it with the stored one."""
    redo = mag_s2.decide(rec, REL, NULL)
    for key in ("orientability_ratio", "shape", "longest"):
        if not same(redo[key], rec.get(key)):
            problems.append(f"{where}: stored {key}={rec.get(key)!r}, derived from the bars {redo[key]!r}")
    if [int(x) for x in rec.get("signature", [])] != redo["signature"]:
        problems.append(f"{where}: stored signature={rec.get('signature')!r}, derived {redo['signature']!r}")
    for z in ("Z2", "Z3"):
        for h in ("H1", "H2"):
            if not same(redo["thresholds"][z][h], rec["thresholds"][z][h]):
                problems.append(f"{where}: threshold {z} {h} stored {rec['thresholds'][z][h]!r}, "
                                f"derived {redo['thresholds'][z][h]!r}")
            if redo["betti"][z][h] != rec["betti"][z][h]:
                problems.append(f"{where}: Betti {z} {h} stored {rec['betti'][z][h]!r}, derived {redo['betti'][z][h]!r}")


def check_design(rec, where, problems):
    """A design-map record stores its signature and ratio only; its label must follow from them."""
    b1, b2 = int(rec["signature"][0]), int(rec["signature"][1])
    ratio = rec["orientability_ratio"]
    orientable = None if b2 == 0 else bool(ratio >= mag_s2.ORIENTABLE_MIN)
    expected = mag_s2.SIGNATURES.get((b1, b2, orientable), "other")
    if rec["shape"] != expected:
        problems.append(f"{where}: stored shape={rec['shape']!r}, signature and ratio give {expected!r}")


def audit(run_dir):
    problems, counts = [], dict(records=0, classified=0, design=0, numbers=0)
    files = sorted(glob.glob(os.path.join(run_dir, "seed_*.json")), key=lambda p: int(p.split("_")[-1][:-5]))
    for path in files:
        name = os.path.basename(path)
        with open(path) as fh:
            seed = json.load(fh)
        counts["records"] += 1
        for where, value in walk_numbers(seed):
            counts["numbers"] += 1
            if isinstance(value, float) and not math.isfinite(value):
                problems.append(f"{name}: {where} is not finite ({value!r})")
        for group in ("conditions", "atlas"):
            for key, rec in seed.get(group, {}).items():
                if "bars" in rec:
                    counts["classified"] += 1
                    try:
                        check_classified(rec, f"{name}: {group}.{key}", problems)
                    except (KeyError, TypeError, ValueError, IndexError) as exc:
                        problems.append(f"{name}: {group}.{key} cannot be derived again ({exc!r})")
        for key, rec in seed.get("design", {}).items():
            counts["design"] += 1
            try:
                check_design(rec, f"{name}: design.{key}", problems)
            except (KeyError, TypeError, ValueError, IndexError) as exc:
                problems.append(f"{name}: design.{key} cannot be checked ({exc!r})")
    return files, counts, problems


def main():
    run_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "confirmatory_run_v01")
    files, counts, problems = audit(run_dir)
    print(f"seed files: {len(files)}; numbers checked: {counts['numbers']}; labels derived again from bars: "
          f"{counts['classified']}; design labels checked against their signature: {counts['design']}")
    for p in problems:
        print("PROBLEM", p)
    print(f"problems: {len(problems)}")
    sys.exit(1 if problems or not files else 0)


if __name__ == "__main__":
    main()
