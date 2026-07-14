"""Detect and repair duplicate rows left behind by GeoPackage writers that
used GDAL's "append" mode to add a layer that already existed.

Background: ``update_gpkg_to_disk`` in ``LEAFs/Name_Cleaning.ipynb`` used to
write the first layer with ``mode="w"`` and every other layer with
``mode="a"``/``append=True``. For GPKG, appending to a layer that already
exists adds rows to it instead of replacing its contents. Re-running that
notebook against an already-populated ``.gpkg`` therefore duplicated every
row in the data tables on each save -- including rows that had just been
set to ``np.nan`` -- leaving the original, un-fixed rows in the file
permanently. This module finds and removes those duplicates in any
already-corrupted ``.gpkg`` file.

Usage:
    python -m sbtn_leaf.gpkg_repair [path] [--dry-run]

``path`` may be a single ``.gpkg`` file or a directory to scan recursively
(default: ``LEAFs``).
"""

from __future__ import annotations

import argparse
import os
from glob import glob
from typing import Iterable

import geopandas as gpd
import pandas as pd
import pyogrio

# Columns holding measured values rather than identity. Two rows that match
# on every other column are the same logical record; if they only differ in
# these columns (e.g. one is NaN and one has the original number), they are
# duplicates from a corrupted append rather than genuinely distinct rows.
VALUE_COLUMNS = {"leaf", "leaf_median", "leaf_std", "cf", "cf_median", "cf_std"}


def load_gpkg(filepath: str) -> dict[str, gpd.GeoDataFrame]:
    layer_info = pyogrio.list_layers(filepath)
    data = {}
    for name, geom_type in layer_info:
        data[name] = gpd.read_file(filepath, layer=name, ignore_geometry=(geom_type is None))
    return data


def write_gpkg(data: dict, gpkg_path: str) -> None:
    """Write every layer fresh. Deletes the file first so no driver-level
    append semantics can silently duplicate rows on the next save."""
    if os.path.exists(gpkg_path):
        os.remove(gpkg_path)
    for name, df in data.items():
        if isinstance(df, gpd.GeoDataFrame):
            df.to_file(gpkg_path, layer=name, driver="GPKG", mode="w")
        else:
            pyogrio.write_dataframe(df, gpkg_path, layer=name, driver="GPKG", append=False)


def dedupe_key(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in VALUE_COLUMNS]


def find_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Rows that share a dedupe key with at least one other row, for inspection."""
    key = dedupe_key(df)
    if not key:
        return df.iloc[0:0]
    return df[df.duplicated(subset=key, keep=False)].sort_values(key)


def repair_gpkg(gpkg_path: str, dry_run: bool = False) -> dict[str, int]:
    """Remove duplicate rows from every non-spatial table in ``gpkg_path``.

    Geometry layers are left untouched -- they were always written with
    real overwrite mode, so they are not the target of this bug, and
    attribute-only deduplication is not safe for geometries (e.g. a single
    ecoregion legitimately split into multiple polygon rows).

    When two rows share the same identity (every column except the value
    columns), the LAST occurrence on disk is kept: each corrupted re-save
    appended the most recently computed in-memory state after the stale
    rows, so the last duplicate is the most up to date one.

    Returns ``{layer_name: n_duplicate_rows_removed}``. Rewrites the file
    in place unless ``dry_run`` is True.
    """
    data = load_gpkg(gpkg_path)
    removed: dict[str, int] = {}
    changed = False

    for name, df in data.items():
        if isinstance(df, gpd.GeoDataFrame):
            continue
        key = dedupe_key(df)
        if not key:
            continue
        before = len(df)
        deduped = df.drop_duplicates(subset=key, keep="last").reset_index(drop=True)
        after = len(deduped)
        if after != before:
            changed = True
            data[name] = deduped
        removed[name] = before - after

    if changed and not dry_run:
        write_gpkg(data, gpkg_path)

    return removed


def iter_gpkgs(root: str) -> Iterable[str]:
    for path in sorted(glob(os.path.join(root, "**", "*.gpkg"), recursive=True)):
        if os.path.getsize(path) > 0:
            yield path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "path", nargs="?", default="LEAFs",
        help="A .gpkg file, or a directory to scan recursively (default: LEAFs)",
    )
    parser.add_argument("--dry-run", action="store_true", help="Report duplicates without rewriting files")
    args = parser.parse_args()

    if os.path.isfile(args.path):
        paths = [args.path]
    else:
        paths = list(iter_gpkgs(args.path))

    if not paths:
        print(f"No non-empty .gpkg files found under {args.path}")
        return

    for path in paths:
        try:
            removed = repair_gpkg(path, dry_run=args.dry_run)
        except Exception as exc:
            print(f"[SKIP] {path}: could not read ({exc})")
            continue
        total = sum(removed.values())
        if total == 0:
            print(f"[OK]   {path}: no duplicates found")
        else:
            verb = "would remove" if args.dry_run else "removed"
            detail = ", ".join(f"{k}={v}" for k, v in removed.items() if v)
            print(f"[FIX]  {path}: {verb} {total} duplicate rows ({detail})")


if __name__ == "__main__":
    main()
