#!/usr/bin/env python3
"""Join batch matrices into one BED-style sample matrix."""

from __future__ import annotations

import argparse
import gc
import os
import tempfile

import polars as pl

from table_io import NA, coordinate_sort, log, read_manifest


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


def read_header(path: str) -> list[str]:
    with open(path, "rt", encoding="utf-8") as handle:
        return handle.readline().rstrip("\n").split("\t")


def read_matrix(path: str, columns: list[str]) -> pl.DataFrame:
    frame = pl.read_csv(
        path,
        separator="\t",
        has_header=True,
        columns=columns,
        infer_schema_length=0,
        null_values=[],
    )
    return frame


def matrix_features(path: str) -> pl.DataFrame:
    header = read_header(path)
    if len(header) < 5:
        raise ValueError(f"{path} does not look like a matrix bed with sample columns")

    frame = read_matrix(path, [header[0], "start", "end", "feature"])
    frame = frame.rename({header[0]: "chrom"})
    return frame.select(["chrom", "start", "end", "feature"])


def matrix_sample_columns(path: str) -> list[str]:
    header = read_header(path)
    if len(header) < 5:
        raise ValueError(f"{path} does not look like a matrix bed with sample columns")
    return header[4:]


def matrix_values(path: str, sample_cols: list[str]) -> pl.DataFrame:
    header = read_header(path)
    frame = read_matrix(path, [header[0], "start", "end", *sample_cols])
    frame = frame.rename({header[0]: "chrom"})
    return frame.select(["chrom", "start", "end", *sample_cols])


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


def write_aligned_sample_block(
    features: pl.DataFrame,
    path: str,
    sample_cols: list[str],
    output: str,
) -> None:
    values = matrix_values(path, sample_cols)
    aligned = features.join(values, on=COORDS, how="left", validate="1:1").sort("_row_idx")
    aligned.select([pl.col(sample).cast(pl.Utf8).fill_null(NA) for sample in sample_cols]).write_csv(
        output,
        separator="\t",
        include_header=False,
    )


def write_final_matrix(base_path: str, block_paths: list[str], output: str, sample_cols: list[str]) -> int:
    rows_written = 0
    with open(base_path, "rt", encoding="utf-8") as base, open(output, "wt", encoding="utf-8") as out:
        block_handles = [open(path, "rt", encoding="utf-8") for path in block_paths]
        try:
            out.write("\t".join(["#chr", "start", "end", "feature", *sample_cols]) + "\n")
            for base_line in base:
                blocks = [handle.readline().rstrip("\n") for handle in block_handles]
                if any(block == "" for block in blocks):
                    raise ValueError("aligned sample block ended before coordinate rows")
                out.write(base_line.rstrip("\n") + "\t" + "\t".join(blocks) + "\n")
                rows_written += 1

            extra = [handle.readline() for handle in block_handles]
            if any(line for line in extra):
                raise ValueError("aligned sample block has more rows than coordinate rows")
        finally:
            for handle in block_handles:
                handle.close()

    return rows_written


def check_unique_sample_columns(sample_cols: list[str]) -> None:
    seen = set()
    duplicates = []
    for sample in sample_cols:
        if sample in seen:
            duplicates.append(sample)
        seen.add(sample)
    if duplicates:
        duplicate_text = ", ".join(sorted(set(duplicates)))
        raise ValueError(f"duplicate sample columns across batch matrices: {duplicate_text}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    paths = [row[0] for row in read_manifest(args.manifest)]
    if not paths:
        raise ValueError("manifest is empty")

    log(f"Joining {len(paths)} batch matrix files")

    sample_cols = []
    matrix_sample_cols = []
    features = empty_features()

    for idx, path in enumerate(paths, start=1):
        log(f"Discovering rows from batch matrix {idx}/{len(paths)}: {path}")
        batch_samples = matrix_sample_columns(path)
        update = matrix_features(path)
        features = merge_features(features, update)
        sample_cols.extend(batch_samples)
        matrix_sample_cols.append(batch_samples)
        log(f"Running coordinate union after batch matrix {idx}/{len(paths)}: {features.height} rows")
        del update
        gc.collect()

    check_unique_sample_columns(sample_cols)

    features = coordinate_sort(features).with_row_index("_row_idx")

    output = f"{args.output_prefix}.bed"

    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        base_path = os.path.join(tmpdir, "base.tsv")
        block_paths = []

        log(f"Writing sorted coordinate/feature columns: {features.height} rows")
        write_base_columns(features, base_path)

        for idx, (path, batch_samples) in enumerate(zip(paths, matrix_sample_cols), start=1):
            block_path = os.path.join(tmpdir, f"batch_{idx}.tsv")
            log(f"Writing aligned sample block {idx}/{len(paths)} with {len(batch_samples)} sample columns")
            write_aligned_sample_block(features, path, batch_samples, block_path)
            block_paths.append(block_path)
            gc.collect()

        log("Concatenating aligned sample blocks on disk")
        rows_written = write_final_matrix(base_path, block_paths, output, sample_cols)

    log(f"Wrote joined matrix with {rows_written} rows and {len(sample_cols)} samples to {output}")


if __name__ == "__main__":
    main()
