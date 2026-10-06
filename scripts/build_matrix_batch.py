#!/usr/bin/env python3
"""Build one batch of a BED-style sample matrix."""

from __future__ import annotations

import argparse
import gc
import os
import tempfile

import polars as pl

from table_io import NA, coordinate_sort, log, read_bed, read_manifest


COORDS = ["chrom", "start", "end"]


def empty_features() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            "chrom": pl.Utf8,
            "start": pl.Utf8,
            "end": pl.Utf8,
            "feature": pl.Utf8,
        }
    )


def sample_features(path: str, feature_col: int) -> pl.DataFrame:
    return read_bed(path).select(
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").alias("start"),
        pl.col("column_3").alias("end"),
        pl.col(f"column_{feature_col}").alias("feature"),
    )


def sample_values(path: str, sample_id: str, value_col: int) -> pl.DataFrame:
    return read_bed(path).select(
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").alias("start"),
        pl.col("column_3").alias("end"),
        pl.col(f"column_{value_col}").alias(sample_id),
    )


def merge_features(features: pl.DataFrame, update: pl.DataFrame) -> pl.DataFrame:
    joined = features.join(update, on=COORDS, how="full", coalesce=True, suffix="_new", validate="1:1")
    return joined.select(
        *COORDS,
        pl.when(pl.col("feature").is_null())
        .then(pl.col("feature_new"))
        .otherwise(pl.col("feature"))
        .alias("feature"),
    )


def write_base_columns(features: pl.DataFrame, output: str) -> None:
    features.select(["chrom", "start", "end", "feature"]).write_csv(output, separator="\t", include_header=False)


def write_sample_column(features: pl.DataFrame, sample_id: str, path: str, value_col: int, output: str) -> None:
    values = sample_values(path, sample_id, value_col)
    aligned = features.join(values, on=COORDS, how="left", validate="1:1").sort("_row_idx")
    aligned.select(pl.col(sample_id).cast(pl.Utf8).fill_null(NA)).write_csv(
        output,
        separator="\t",
        include_header=False,
    )


def write_final_matrix(base_path: str, sample_paths: list[str], output: str, sample_ids: list[str]) -> int:
    rows_written = 0
    with open(base_path, "rt", encoding="utf-8") as base, open(output, "wt", encoding="utf-8") as out:
        sample_handles = [open(path, "rt", encoding="utf-8") for path in sample_paths]
        try:
            out.write("\t".join(["#chr", "start", "end", "feature", *sample_ids]) + "\n")
            for base_line in base:
                values = [handle.readline().rstrip("\n") for handle in sample_handles]
                if any(value == "" for value in values):
                    raise ValueError("sample column file ended before coordinate rows")
                out.write(base_line.rstrip("\n") + "\t" + "\t".join(values) + "\n")
                rows_written += 1

            extra = [handle.readline() for handle in sample_handles]
            if any(line for line in extra):
                raise ValueError("sample column file has more rows than coordinate rows")
        finally:
            for handle in sample_handles:
                handle.close()

    return rows_written


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

    features = empty_features()
    for idx, (sample_id, path) in enumerate(entries, start=1):
        log(f"Discovering rows from BED {idx}/{len(entries)}: sample={sample_id} path={path}")
        update = sample_features(path, args.feature_col)
        features = merge_features(features, update)
        log(f"Running coordinate union after BED {idx}/{len(entries)}: {features.height} rows")
        del update
        gc.collect()

    features = coordinate_sort(features).with_row_index("_row_idx")
    sample_ids = [sample_id for sample_id, _path in entries]

    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        base_path = os.path.join(tmpdir, "base.tsv")
        sample_paths = []

        log(f"Writing sorted coordinate/feature columns: {features.height} rows")
        write_base_columns(features, base_path)

        for idx, (sample_id, path) in enumerate(entries, start=1):
            column_path = os.path.join(tmpdir, f"sample_{idx}.tsv")
            log(f"Writing aligned sample column {idx}/{len(entries)}: {sample_id}")
            write_sample_column(features, sample_id, path, args.value_col, column_path)
            sample_paths.append(column_path)
            gc.collect()

        log("Concatenating matrix columns on disk")
        rows_written = write_final_matrix(base_path, sample_paths, args.output, sample_ids)

    log(f"Wrote matrix batch with {rows_written} rows and {len(sample_ids)} samples to {args.output}")


if __name__ == "__main__":
    main()
