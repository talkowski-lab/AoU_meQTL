#!/usr/bin/env python3
"""Combine reduced CpG batch statistics and write the final BED summary."""

from __future__ import annotations

import argparse
import gc

import polars as pl

from summary_common import coord_name, empty_cpg_stats, merge_cpg_stats
from table_io import NA, coordinate_sort, log, read_manifest


def read_cpg_reduced(path: str) -> pl.DataFrame:
    return (
        pl.read_csv(path, separator="\t", has_header=False, infer_schema_length=0)
        .rename(
            {
                "column_1": "chrom",
                "column_2": "start",
                "column_3": "end",
                "column_4": "n",
                "column_5": "sum",
                "column_6": "sumsq",
                "column_7": "min_val",
                "column_8": "max_val",
            }
        )
        .with_columns(
            pl.col("n").cast(pl.Int64),
            pl.col("sum").cast(pl.Float64),
            pl.col("sumsq").cast(pl.Float64),
            pl.col("min_val").cast(pl.Float64),
            pl.col("max_val").cast(pl.Float64),
        )
    )


def write_cpg_summary(stats: pl.DataFrame, output: str) -> None:
    stats = coordinate_sort(stats)
    with open(output, "wt", encoding="utf-8") as out:
        out.write("#chrom\tstart\tend\tname\tn\tmean\tvariance\tmin\tmax\n")
        for chrom, start, end, n, total, sumsq, min_val, max_val in stats.iter_rows():
            mean = total / n
            variance = (sumsq - (total * total) / n) / (n - 1) if n > 1 else NA
            out.write(
                f"{chrom}\t{start}\t{end}\t{coord_name(chrom, start, end)}\t{n}\t"
                f"{mean}\t{variance}\t{min_val}\t{max_val}\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    paths = [row[0] for row in read_manifest(args.manifest)]
    if not paths:
        raise ValueError("manifest is empty")

    log(f"Combining CpG statistics from {len(paths)} batch files")

    stats = empty_cpg_stats()
    for idx, path in enumerate(paths, start=1):
        log(f"Reading CpG batch stats {idx}/{len(paths)}: {path}")
        batch = read_cpg_reduced(path)
        log(f"Merging CpG batch stats {idx}/{len(paths)}: {batch.height} rows")
        stats = merge_cpg_stats(stats, batch)
        log(f"Running CpG combined stats after batch {idx}/{len(paths)}: {stats.height} coordinate rows")
        del batch
        gc.collect()

    output = f"{args.output_prefix}.bed"
    write_cpg_summary(stats, output)
    log(f"Wrote CpG summary BED with {stats.height} rows to {output}")


if __name__ == "__main__":
    main()
