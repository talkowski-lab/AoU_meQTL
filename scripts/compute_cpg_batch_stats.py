#!/usr/bin/env python3
"""Compute reduced CpG statistics for one batch of BED files."""

from __future__ import annotations

import argparse
import gc

import polars as pl

from summary_common import merge_cpg_stats, write_reduced_stats
from table_io import log, numeric, read_bed, read_manifest


def cpg_contribution(path: str, value_col: int) -> pl.DataFrame:
    value = numeric(pl.col(f"column_{value_col}")).alias("value")
    frame = read_bed(path).select(
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").alias("start"),
        pl.col("column_3").alias("end"),
        value,
    )
    return frame.filter(pl.col("value").is_not_null()).select(
        "chrom",
        "start",
        "end",
        pl.lit(1).alias("n"),
        pl.col("value").alias("sum"),
        (pl.col("value") * pl.col("value")).alias("sumsq"),
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

    log(f"Computing CpG batch statistics from {len(paths)} BED files")

    stats = None
    for idx, path in enumerate(paths, start=1):
        log(f"Reading BED {idx}/{len(paths)}: {path}")
        contribution = cpg_contribution(path, args.value_col)
        log(f"Merging BED {idx}/{len(paths)}: {contribution.height} usable rows")
        stats = merge_cpg_stats(stats, contribution)
        log(f"Running CpG summary after BED {idx}/{len(paths)}: {stats.height} coordinate rows")
        del contribution
        gc.collect()

    write_reduced_stats(stats, args.output, ["chrom", "start", "end", "n", "sum", "sumsq", "min_val", "max_val"])
    log(f"Wrote reduced CpG batch statistics with {stats.height} rows to {args.output}")


if __name__ == "__main__":
    main()
