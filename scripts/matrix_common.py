"""Shared Polars helpers for matrix-building scripts."""

from __future__ import annotations

import polars as pl


def first_features(frames: list[pl.DataFrame]) -> pl.DataFrame:
    return (
        pl.concat([frame.select(["chrom", "start", "end", "feature", "_source_idx"]) for frame in frames])
        .sort("_source_idx")
        .unique(subset=["chrom", "start", "end"], keep="first")
        .drop("_source_idx")
    )

