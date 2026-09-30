#!/usr/bin/env python3
"""Aggregate per-CpG calls into per-window methylation summary stats.

Reads, on stdin, the output of `bedtools intersect -wa -wb -a windows.bed -b
cpg.bed`: each line is a 3-column window (chrom, start, end) followed by the
full CpG bed record it overlaps. Writes one aggregated row per window (in
first-seen order) to stdout: chrom, start, end, num_cpgs, total_coverage,
weighted_mean_methylation, unweighted_mean_methylation.
"""

import argparse
import sys

WINDOW_COLS = 3


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--meth-col", type=int, default=4,
        help="1-based column of the methylation percentage (0-100) in the CpG bed "
             "(default: 4, pb-CpG-tools' modification_probability column)",
    )
    parser.add_argument(
        "--cov-col", type=int, default=6,
        help="1-based column of the coverage in the CpG bed "
             "(default: 6, pb-CpG-tools' coverage column)",
    )
    args = parser.parse_args()

    meth_idx = WINDOW_COLS + args.meth_col - 1
    cov_idx = WINDOW_COLS + args.cov_col - 1

    windows = {}
    order = []
    for line in sys.stdin:
        fields = line.rstrip("\n").split("\t")
        key = (fields[0], fields[1], fields[2])
        meth = float(fields[meth_idx])
        cov = int(float(fields[cov_idx]))

        if key not in windows:
            windows[key] = {"n": 0, "cov_sum": 0, "weighted_sum": 0.0, "meth_sum": 0.0}
            order.append(key)

        stats = windows[key]
        stats["n"] += 1
        stats["cov_sum"] += cov
        stats["weighted_sum"] += meth * cov
        stats["meth_sum"] += meth

    for key in order:
        stats = windows[key]
        weighted_mean = stats["weighted_sum"] / stats["cov_sum"] if stats["cov_sum"] > 0 else "NA"
        unweighted_mean = stats["meth_sum"] / stats["n"]
        print("\t".join([
            key[0], key[1], key[2],
            str(stats["n"]),
            str(stats["cov_sum"]),
            f"{weighted_mean:.4f}" if weighted_mean != "NA" else "NA",
            f"{unweighted_mean:.4f}",
        ]))


if __name__ == "__main__":
    main()
