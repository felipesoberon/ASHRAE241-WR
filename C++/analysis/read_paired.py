#!/usr/bin/env python3
"""Read a completed ECAIP01 dataset and write per-replicate percentiles as CSV."""

import argparse
import csv
import math
from pathlib import Path
import struct
import sys


def percentiles_argument(text):
    try:
        values = [float(value) for value in text.split(",")]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated percentiles") from exc
    if not values or any(not math.isfinite(p) or not 0 <= p <= 1 for p in values):
        raise argparse.ArgumentTypeError("percentiles must be finite fractions in [0, 1]")
    return values


def exact(stream, count):
    data = stream.read(count)
    if len(data) != count:
        raise ValueError("truncated paired dataset")
    return data


def read_percentiles(directory, percentiles, output):
    directory = Path(directory)
    manifest_path = directory / "manifest.txt"
    if not manifest_path.is_file():
        raise ValueError("missing completed manifest.txt; dataset may be incomplete")
    manifest = {}
    for line in manifest_path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition("=")
        if not separator or not key or key in manifest:
            raise ValueError("malformed or duplicate manifest field")
        manifest[key] = value
    try:
        manifest_targets = [float(x) for x in manifest["targets"].split(",")]
        manifest_simulations = int(manifest["simulations"])
        manifest_replicates = int(manifest["replicates"])
        if manifest["format"] != "ECAIP01" or manifest["layout"] != "replicate,target,simulation":
            raise ValueError("unsupported manifest format or layout")
    except (KeyError, ValueError) as exc:
        raise ValueError("incomplete or invalid paired manifest") from exc
    path = directory / "ecai.bin"
    with path.open("rb") as stream:
        if exact(stream, 8) != b"ECAIP01\n":
            raise ValueError("unsupported paired dataset format")
        target_count, simulations, replicates = struct.unpack("<QQQ", exact(stream, 24))
        if not 1 <= target_count <= 64 or not 1 <= simulations <= (1 << 29) or replicates < 1:
            raise ValueError("invalid paired dataset dimensions")
        if simulations > (1 << 30) // (8 * target_count):
            raise ValueError("paired dataset exceeds per-replicate limit")
        expected_size = 32 + 8 * target_count + 8 * target_count * simulations * replicates
        if path.stat().st_size != expected_size:
            raise ValueError("paired dataset has incorrect size or is incomplete")
        targets = struct.unpack(f"<{target_count}d", exact(stream, target_count * 8))
        if len(set(targets)) != len(targets) or any(not 0 < t < 1 for t in targets):
            raise ValueError("invalid target probabilities in dataset")
        if (list(targets) != manifest_targets or simulations != manifest_simulations or
                replicates != manifest_replicates):
            raise ValueError("manifest and paired dataset dimensions or targets disagree")
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(("replicate_index", "target_probability", "percentile", "ecai_l_s_person"))
        for replicate in range(replicates):
            for target in targets:
                values = list(struct.unpack(f"<{simulations}d", exact(stream, simulations * 8)))
                if any(not math.isfinite(value) or value < 0 for value in values):
                    raise ValueError("invalid ECAi value in paired dataset")
                values.sort()
                for p in percentiles:
                    rank = p * (simulations - 1)
                    lo = math.floor(rank)
                    hi = math.ceil(rank)
                    fraction = rank - lo
                    value = values[lo] * (1.0 - fraction) + values[hi] * fraction
                    writer.writerow((replicate, repr(target), repr(p), repr(value)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", help="completed paired-output directory")
    parser.add_argument("--percentiles", type=percentiles_argument, required=True,
                        help="comma-separated fractions, for example 0.75,0.9,0.96")
    args = parser.parse_args()
    try:
        read_percentiles(args.directory, args.percentiles, sys.stdout)
    except (OSError, ValueError, OverflowError) as exc:
        parser.exit(1, f"read_paired: {exc}\n")


if __name__ == "__main__":
    main()
