#!/usr/bin/env python3
"""Build one batch of a BED-style sample matrix."""

from __future__ import annotations

import argparse

import polars as pl

from matrix_common import first_features
from table_io import coordinate_sort, log, read_bed, read_manifest, write_bed


def sample_frame(path: str, sample_id: str, sample_idx: int, feature_col: int, value_col: int) -> pl.DataFrame:
    frame = read_bed(path)
    selected = frame.select(
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").alias("start"),
        pl.col("column_3").alias("end"),
        pl.col(f"column_{feature_col}").alias("feature"),
        pl.col(f"column_{value_col}").alias(sample_id),
    )
    return selected.with_columns(pl.lit(sample_idx).alias("_source_idx"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--feature-col", type=int, required=True)
    parser.add_argument("--value-col", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    entries = read_manifest(args.manifest)
    if not entries:
        raise ValueError("manifest is empty")

    log(f"Building matrix batch from {len(entries)} BED files")

    frames = []
    for idx, (sample_id, path) in enumerate(entries, start=1):
        log(f"Reading BED {idx}/{len(entries)}: sample={sample_id} path={path}")
        frame = sample_frame(path, sample_id, idx - 1, args.feature_col, args.value_col)
        frames.append(frame)

    log("Collecting coordinate union for batch matrix")
    matrix = first_features(frames)

    sample_ids = []
    for idx, (frame, (sample_id, _path)) in enumerate(zip(frames, entries), start=1):
        log(f"Joining sample {idx}/{len(entries)}: {sample_id}")
        sample_ids.append(sample_id)
        matrix = matrix.join(
            frame.select(["chrom", "start", "end", sample_id]),
            on=["chrom", "start", "end"],
            how="left",
        )

    matrix = coordinate_sort(matrix.select(["chrom", "start", "end", "feature", *sample_ids]))
    write_bed(matrix, args.output, ["#chr", "start", "end", "feature", *sample_ids])
    log(f"Wrote matrix batch with {matrix.height} rows and {len(sample_ids)} samples to {args.output}")


if __name__ == "__main__":
    main()
