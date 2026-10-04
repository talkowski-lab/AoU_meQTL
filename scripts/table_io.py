"""Shared IO helpers for BED-like TSV files used by the workflow scripts."""

from __future__ import annotations

import sys

import polars as pl


NA = "NA"


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def read_manifest(path: str) -> list[list[str]]:
    rows: list[list[str]] = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line:
                rows.append(line.split("\t"))
    return rows


def read_bed(path: str) -> pl.DataFrame:
    return pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        comment_prefix="#",
        infer_schema_length=0,
        null_values=[],
    )


def coordinate_sort(frame: pl.DataFrame) -> pl.DataFrame:
    return (
        frame.with_columns(
            pl.col("start").cast(pl.Int64).alias("_start_sort"),
            pl.col("end").cast(pl.Int64).alias("_end_sort"),
        )
        .sort(["chrom", "_start_sort", "_end_sort"])
        .drop(["_start_sort", "_end_sort"])
    )


def write_bed(frame: pl.DataFrame, output: str, header: list[str]) -> None:
    frame = frame.with_columns(pl.all().cast(pl.Utf8)).fill_null(NA)
    with open(output, "wt", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        frame.write_csv(handle, separator="\t", include_header=False)


def numeric(expr: pl.Expr) -> pl.Expr:
    return expr.cast(pl.Float64, strict=False)
