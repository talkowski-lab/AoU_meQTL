#!/usr/bin/env python3
"""Combine reduced bin batch statistics and write the final BED summary."""

from __future__ import annotations

import argparse
import gc

import polars as pl

from summary_common import coord_name, empty_bin_stats, merge_bin_stats
from table_io import NA, coordinate_sort, log, read_manifest


def read_bin_reduced(path: str) -> pl.DataFrame:
    return (
        pl.read_csv(path, separator="\t", has_header=False, infer_schema_length=0)
        .rename(
            {
                "column_1": "chrom",
                "column_2": "start",
                "column_3": "end",
                "column_4": "n",
                "column_5": "sum_cpgs",
                "column_6": "sum_val",
                "column_7": "sumsq_val",
                "column_8": "min_val",
                "column_9": "max_val",
            }
        )
        .with_columns(
            pl.col("n").cast(pl.Int64),
            pl.col("sum_cpgs").cast(pl.Float64),
            pl.col("sum_val").cast(pl.Float64),
            pl.col("sumsq_val").cast(pl.Float64),
            pl.col("min_val").cast(pl.Float64),
            pl.col("max_val").cast(pl.Float64),
        )
    )


def write_bin_summary(stats: pl.DataFrame, output: str, total_samples: int) -> None:
    stats = coordinate_sort(stats)
    with open(output, "wt", encoding="utf-8") as out:
        out.write(
            "#chrom\tstart\tend\tname\tn\tmean_value\tvariance_value\tmin_value\tmax_value\t"
            "delta\tmean_num_cpgs\tpresence\n"
        )
        for chrom, start, end, n, sum_cpgs, sum_val, sumsq_val, min_val, max_val in stats.iter_rows():
            mean_value = sum_val / n
            variance_value = (sumsq_val - (sum_val * sum_val) / n) / (n - 1) if n > 1 else NA
            delta = max_val - min_val
            mean_num_cpgs = sum_cpgs / n
            presence = n / total_samples
            out.write(
                f"{chrom}\t{start}\t{end}\t{coord_name(chrom, start, end)}\t{n}\t{mean_value}\t"
                f"{variance_value}\t{min_val}\t{max_val}\t{delta}\t{mean_num_cpgs}\t{presence}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--total-samples", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    paths = [row[0] for row in read_manifest(args.manifest)]
    if not paths:
        raise ValueError("manifest is empty")

    log(f"Combining bin statistics from {len(paths)} batch files")

    stats = empty_bin_stats()
    for idx, path in enumerate(paths, start=1):
        log(f"Reading bin batch stats {idx}/{len(paths)}: {path}")
        batch = read_bin_reduced(path)
        log(f"Merging bin batch stats {idx}/{len(paths)}: {batch.height} rows")
        stats = merge_bin_stats(stats, batch)
        log(f"Running bin combined stats after batch {idx}/{len(paths)}: {stats.height} coordinate rows")
        del batch
        gc.collect()

    write_bin_summary(stats, args.output, args.total_samples)
    log(f"Wrote bin summary BED with {stats.height} rows to {args.output}")


if __name__ == "__main__":
    main()
