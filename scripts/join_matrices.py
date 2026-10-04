#!/usr/bin/env python3
"""Join batch matrices into one BED-style sample matrix."""

from __future__ import annotations

import argparse

import polars as pl

from matrix_common import first_features
from table_io import coordinate_sort, log, read_manifest, write_bed


def read_matrix(path: str) -> tuple[list[str], pl.DataFrame]:
    with open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    frame = pl.read_csv(
        path,
        separator="\t",
        has_header=True,
        infer_schema_length=0,
        null_values=[],
    )
    return header, frame


def batch_matrix_frame(path: str, batch_idx: int) -> tuple[list[str], pl.DataFrame]:
    header, frame = read_matrix(path)
    if len(header) < 5:
        raise ValueError(f"{path} does not look like a matrix bed with sample columns")

    frame = frame.rename({header[0]: "chrom"})
    sample_cols = header[4:]
    selected = frame.select(["chrom", "start", "end", "feature", *sample_cols])
    return sample_cols, selected.with_columns(pl.lit(batch_idx).alias("_source_idx"))


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
    frames = []
    for idx, path in enumerate(paths, start=1):
        log(f"Reading batch matrix {idx}/{len(paths)}: {path}")
        batch_samples, frame = batch_matrix_frame(path, idx - 1)
        sample_cols.extend(batch_samples)
        frames.append(frame)

    log("Collecting coordinate union across batch matrices")
    matrix = first_features(frames)

    for idx, frame in enumerate(frames, start=1):
        batch_samples = [c for c in frame.columns if c not in {"chrom", "start", "end", "feature", "_source_idx"}]
        log(f"Joining batch matrix {idx}/{len(frames)} with {len(batch_samples)} sample columns")
        matrix = matrix.join(
            frame.select(["chrom", "start", "end", *batch_samples]),
            on=["chrom", "start", "end"],
            how="left",
        )

    output = f"{args.output_prefix}.bed"
    matrix = coordinate_sort(matrix.select(["chrom", "start", "end", "feature", *sample_cols]))
    write_bed(matrix, output, ["#chr", "start", "end", "feature", *sample_cols])
    log(f"Wrote joined matrix with {matrix.height} rows and {len(sample_cols)} samples to {output}")


if __name__ == "__main__":
    main()
