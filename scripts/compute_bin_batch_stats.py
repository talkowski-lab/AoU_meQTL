#!/usr/bin/env python3
"""Compute reduced bin statistics for one batch of BED files."""

from __future__ import annotations

import argparse
import gc

import polars as pl

from summary_common import merge_bin_stats, write_reduced_stats
from table_io import log, numeric, read_bed, read_manifest


def bin_contribution(path: str, value_col: int) -> pl.DataFrame:
    frame = read_bed(path).select(
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").alias("start"),
        pl.col("column_3").alias("end"),
        numeric(pl.col("column_5")).alias("num_cpgs"),
        numeric(pl.col(f"column_{value_col}")).alias("value"),
    )
    return frame.filter(pl.col("value").is_not_null() & pl.col("num_cpgs").is_not_null()).select(
        "chrom",
        "start",
        "end",
        pl.lit(1).alias("n"),
        pl.col("num_cpgs").alias("sum_cpgs"),
        pl.col("value").alias("sum_val"),
        (pl.col("value") * pl.col("value")).alias("sumsq_val"),
        pl.col("value").alias("min_val"),
        pl.col("value").alias("max_val"),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--value-col", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    paths = [row[0] for row in read_manifest(args.manifest)]
    if not paths:
        raise ValueError("manifest is empty")

    log(f"Computing bin batch statistics from {len(paths)} BED files")

    stats = None
    for idx, path in enumerate(paths, start=1):
        log(f"Reading BED {idx}/{len(paths)}: {path}")
        contribution = bin_contribution(path, args.value_col)
        log(f"Merging BED {idx}/{len(paths)}: {contribution.height} usable rows")
        stats = merge_bin_stats(stats, contribution)
        log(f"Running bin summary after BED {idx}/{len(paths)}: {stats.height} coordinate rows")
        del contribution
        gc.collect()

    write_reduced_stats(
        stats,
        args.output,
        ["chrom", "start", "end", "n", "sum_cpgs", "sum_val", "sumsq_val", "min_val", "max_val"],
    )
    log(f"Wrote reduced bin batch statistics with {stats.height} rows to {args.output}")


if __name__ == "__main__":
    main()
