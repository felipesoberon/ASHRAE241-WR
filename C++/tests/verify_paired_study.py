#!/usr/bin/env python3
"""Run representative paired/baseline comparisons in a new external directory."""

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import struct
import subprocess
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("read_paired", ROOT / "analysis" / "read_paired.py")
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)


def execute(command, expect_success=True):
    result = subprocess.run([str(part) for part in command], text=True, capture_output=True)
    if (result.returncode == 0) != expect_success:
        raise RuntimeError(f"unexpected exit {result.returncode}: {command}\n{result.stdout}\n{result.stderr}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--updated", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True,
                        help="new directory under the separate data folder")
    parser.add_argument("--simulations", type=int, default=10000)
    parser.add_argument("--replicates", type=int, default=50)
    args = parser.parse_args()
    args.destination.mkdir(parents=False, exist_ok=False)
    targets = "0.001,0.01,0.05,0.1,0.25,0.5"
    evidence = []
    for category, slug in (("Office", "office"), ("Patient room", "patient"),
                           ("Convention", "convention")):
        for policy in ("community", "risk-averse"):
            stem = f"{slug}-{policy}"
            common = ["--category", category, "--simulations", args.simulations,
                      "--replicates", args.replicates, "--infector-policy", policy]
            baseline_dir = args.destination / (stem + "-baseline")
            paired_dir = args.destination / (stem + "-paired")
            execute([args.baseline, *common, "--target-probability", "0.001",
                     "--output-directory", baseline_dir, "--write-replicates"])
            execute([args.updated, *common, "--target-probabilities", targets,
                     "--raw-output-directory", paired_dir])
            with (baseline_dir / "replicates.csv").open(newline="") as stream:
                baseline = list(csv.DictReader(stream))
            output = io.StringIO()
            reader.read_percentiles(paired_dir, [0, 0.5, 0.96, 1], output)
            rows = list(csv.DictReader(io.StringIO(output.getvalue())))
            assert len(rows) == args.replicates * 6 * 4
            comparisons = 0
            for row in rows:
                if float(row["target_probability"]) == 0.001 and float(row["percentile"]) == 0.96:
                    idx = int(row["replicate_index"])
                    assert float(row["ecai_l_s_person"]) == float(baseline[idx]["raw_percentile_l_s_person"]), (stem, idx)
                    comparisons += 1
            assert comparisons == args.replicates
            # Every higher risk target must require no more ECAi for the same row.
            with (paired_dir / "ecai.bin").open("rb") as stream:
                assert stream.read(8) == b"ECAIP01\n"
                count, n, r = struct.unpack("<QQQ", stream.read(24))
                assert (count, n, r) == (6, args.simulations, args.replicates)
                stream.read(count * 8)
                for _ in range(r):
                    previous = None
                    for _ in range(count):
                        current = struct.unpack(f"<{n}d", stream.read(n * 8))
                        assert all(math.isfinite(value) and value >= 0 for value in current)
                        if previous is not None:
                            assert all(a >= b for a, b in zip(previous, current)), stem
                        previous = current
            evidence.append({"category": category, "policy": policy,
                             "baseline_percentiles_exact": comparisons,
                             "monotonic_rows": args.simulations * args.replicates})
            print(f"PASS {stem}: {comparisons} exact replicate percentiles; "
                  f"{args.simulations * args.replicates} paired rows", flush=True)

    # The list tail order must not alter the primary seed or matching columns.
    common = ["--category", "Patient room", "--simulations", 200, "--replicates", 2,
              "--infector-policy", "risk-averse"]
    for name, order in (("order-a", "0.001,0.01,0.5"), ("order-b", "0.001,0.5,0.01"),
                        ("order-c", "0.01,0.001,0.5"), ("one-target", "0.001")):
        execute([args.updated, *common, "--target-probabilities", order,
                 "--raw-output-directory", args.destination / name])
    def matrix(name):
        content = io.StringIO()
        reader.read_percentiles(args.destination / name, [0, 0.5, 0.96, 1], content)
        return {(r["replicate_index"], r["target_probability"], r["percentile"]):
                r["ecai_l_s_person"] for r in csv.DictReader(io.StringIO(content.getvalue()))}
    assert matrix("order-a") == matrix("order-b")
    print("PASS tail target ordering preserves paired results", flush=True)
    assert matrix("order-a") != matrix("order-c")
    assert all(matrix("order-a")[key] == value for key, value in matrix("one-target").items())
    print("PASS first target controls sampling and one-target list matches", flush=True)

    # Existing output must not be overwritten, including on a repeat invocation.
    output = args.destination / "order-a"
    before = hashlib.sha256((output / "ecai.bin").read_bytes()).digest()
    execute([args.updated, *common, "--target-probabilities", "0.001,0.01,0.5",
             "--raw-output-directory", output], expect_success=False)
    assert hashlib.sha256((output / "ecai.bin").read_bytes()).digest() == before
    print("PASS existing paired output is not overwritten", flush=True)
    broken = args.destination / "truncated-dataset"
    broken.mkdir()
    shutil.copyfile(output / "manifest.txt", broken / "manifest.txt")
    with (output / "ecai.bin").open("rb") as source, (broken / "ecai.bin").open("wb") as dest:
        dest.write(source.read(40))
    try:
        reader.read_percentiles(broken, [0.96], io.StringIO())
    except ValueError:
        pass
    else:
        raise AssertionError("reader accepted truncated dataset")
    incomplete = args.destination / "incomplete-dataset"
    incomplete.mkdir()
    try:
        reader.read_percentiles(incomplete, [0.96], io.StringIO())
    except ValueError:
        pass
    else:
        raise AssertionError("reader accepted incomplete dataset")
    print("PASS incomplete and truncated datasets are rejected", flush=True)
    (args.destination / "verification.json").write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, OSError, RuntimeError, ValueError) as exc:
        print(f"verification failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
