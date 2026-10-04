"""Shared helpers for the one-file-at-a-time summary-stat scripts."""

from __future__ import annotations

import polars as pl

from table_io import coordinate_sort


COORDS = ["chrom", "start", "end"]


def merge_cpg_stats(stats: pl.DataFrame | None, update: pl.DataFrame) -> pl.DataFrame:
    if stats is None:
        return update

    joined = stats.join(update, on=COORDS, how="full", coalesce=True, suffix="_new")
    return joined.select(
        *COORDS,
        (pl.col("n").fill_null(0) + pl.col("n_new").fill_null(0)).alias("n"),
        (pl.col("sum").fill_null(0.0) + pl.col("sum_new").fill_null(0.0)).alias("sum"),
        (pl.col("sumsq").fill_null(0.0) + pl.col("sumsq_new").fill_null(0.0)).alias("sumsq"),
        pl.when(pl.col("min_val").is_null())
        .then(pl.col("min_val_new"))
        .when(pl.col("min_val_new").is_null())
        .then(pl.col("min_val"))
        .otherwise(pl.min_horizontal("min_val", "min_val_new"))
        .alias("min_val"),
        pl.when(pl.col("max_val").is_null())
        .then(pl.col("max_val_new"))
        .when(pl.col("max_val_new").is_null())
        .then(pl.col("max_val"))
        .otherwise(pl.max_horizontal("max_val", "max_val_new"))
        .alias("max_val"),
    )


def merge_bin_stats(stats: pl.DataFrame | None, update: pl.DataFrame) -> pl.DataFrame:
    if stats is None:
        return update

    joined = stats.join(update, on=COORDS, how="full", coalesce=True, suffix="_new")
    return joined.select(
        *COORDS,
        (pl.col("n").fill_null(0) + pl.col("n_new").fill_null(0)).alias("n"),
        (pl.col("sum_cpgs").fill_null(0.0) + pl.col("sum_cpgs_new").fill_null(0.0)).alias("sum_cpgs"),
        (pl.col("sum_val").fill_null(0.0) + pl.col("sum_val_new").fill_null(0.0)).alias("sum_val"),
        (pl.col("sumsq_val").fill_null(0.0) + pl.col("sumsq_val_new").fill_null(0.0)).alias("sumsq_val"),
        pl.when(pl.col("min_val").is_null())
        .then(pl.col("min_val_new"))
        .when(pl.col("min_val_new").is_null())
        .then(pl.col("min_val"))
        .otherwise(pl.min_horizontal("min_val", "min_val_new"))
        .alias("min_val"),
        pl.when(pl.col("max_val").is_null())
        .then(pl.col("max_val_new"))
        .when(pl.col("max_val_new").is_null())
        .then(pl.col("max_val"))
        .otherwise(pl.max_horizontal("max_val", "max_val_new"))
        .alias("max_val"),
    )


def write_reduced_stats(stats: pl.DataFrame, output: str, columns: list[str]) -> None:
    coordinate_sort(stats.select(columns)).write_csv(output, separator="\t", include_header=False)


def coord_name(chrom: str, start: str, end: str) -> str:
    return f"{chrom}:{int(start) + 1}-{end}"

