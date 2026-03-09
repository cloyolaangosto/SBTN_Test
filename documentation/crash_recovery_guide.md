# Crash Recovery Guide for `build_cfs_gpkg_from_rasters`

`build_cfs_gpkg_from_rasters` now supports **per-raster checkpointing**, so if the
process crashes or is interrupted, it can resume from where it left off instead of
reprocessing everything from scratch.

## How It Works

Each time a raster is successfully processed, its results are saved as a small
Parquet file inside a `_checkpoints_{layer_name}/` directory in the output folder.
On the next run, if a checkpoint file exists for a given raster, the function loads
the cached result instantly instead of re-running the expensive computation.

Checkpoints are written **atomically** (write to a `.tmp` file, then rename), so a
crash mid-write cannot produce a corrupt checkpoint.

Additionally, the pre-reprojected shapefile is now passed to the inner calculator
function, eliminating redundant shapefile loading and CRS reprojection on every
raster (a significant speedup for subcountry-level runs).

---

## Usage Scenarios

### Scenario 1: Fresh Run (default behavior, fully backward-compatible)

```python
from sbtn_leaf.map_calculations import build_cfs_gpkg_from_rasters

gpkg_path, results_df = build_cfs_gpkg_from_rasters(
    input_folder="/path/to/rasters/",
    output_folder="/path/to/output/",
    layer_name="soil_erosion",
    master_gdf=my_master_gdf,
    master_key="ADM0_NAME",
    result_key="country",
    cf_name="Soil Erosion",
    cf_unit="t/ha/yr",
    area_type="country",
    reset_gpkg=True,          # default: starts fresh, deletes old gpkg + checkpoints
)
```

**What happens:**
- Old GeoPackage and checkpoint directory are deleted
- All rasters are processed from scratch
- Checkpoint files are saved as each raster completes
- If interrupted, the checkpoints survive for the next run (with `reset_gpkg=False`)

---

### Scenario 2: Resume After a Crash

If a run was interrupted (crash, timeout, manual stop), call the function again
with `reset_gpkg=False` to resume:

```python
gpkg_path, results_df = build_cfs_gpkg_from_rasters(
    input_folder="/path/to/rasters/",
    output_folder="/path/to/output/",
    layer_name="soil_erosion",
    master_gdf=my_master_gdf,
    master_key="ADM0_NAME",
    result_key="country",
    cf_name="Soil Erosion",
    cf_unit="t/ha/yr",
    area_type="country",
    reset_gpkg=False,         # <-- keep checkpoints, resume where we left off
)
```

**What happens:**
- Existing checkpoints are preserved
- The GeoPackage is rebuilt cleanly (deleted and re-created from all results)
- Rasters that already have a checkpoint file are loaded instantly (milliseconds)
- Only rasters without checkpoints are computed (the expensive step)
- The final GeoPackage and CSV are identical to a full clean run

---

### Scenario 3: Force Full Recomputation

To discard all cached results and start completely fresh:

```python
gpkg_path, results_df = build_cfs_gpkg_from_rasters(
    input_folder="/path/to/rasters/",
    output_folder="/path/to/output/",
    layer_name="soil_erosion",
    master_gdf=my_master_gdf,
    master_key="ADM0_NAME",
    result_key="country",
    cf_name="Soil Erosion",
    cf_unit="t/ha/yr",
    area_type="country",
    reset_gpkg=True,          # deletes gpkg AND checkpoint directory
)
```

This is the default behavior and is identical to how the function worked before
the crash recovery feature was added.

---

### Scenario 4: Large Subcountry Run with Resume Safety

For subcountry runs with thousands of regions and many rasters, the risk of
interruption is highest. Set up the initial run, then if interrupted, just
re-run with `reset_gpkg=False`:

```python
# First run (may take hours)
gpkg_path, results_df = build_cfs_gpkg_from_rasters(
    input_folder="/data/acidification_rasters/",
    output_folder="/output/acid_subcountry/",
    layer_name="acidification_subcountry",
    master_gdf=subcountry_gdf,
    master_key="ADM1_CODE",
    result_key="ADM1_CODE",
    cf_name="Acidification",
    cf_unit="mol H+/ha/yr",
    area_type="subcountry",
    reset_gpkg=True,          # fresh start
)

# If interrupted at raster 45 of 100, re-run:
gpkg_path, results_df = build_cfs_gpkg_from_rasters(
    input_folder="/data/acidification_rasters/",
    output_folder="/output/acid_subcountry/",
    layer_name="acidification_subcountry",
    master_gdf=subcountry_gdf,
    master_key="ADM1_CODE",
    result_key="ADM1_CODE",
    cf_name="Acidification",
    cf_unit="mol H+/ha/yr",
    area_type="subcountry",
    reset_gpkg=False,         # resume from raster 46
)
```

**What happens on resume:**
- Rasters 1-45: loaded from checkpoint in milliseconds each
- Rasters 46-100: computed normally, checkpointed as they complete
- Final output is identical to a successful uninterrupted run

---

### Scenario 5: CSV-Only Mode with Crash Recovery

Crash recovery also works when `write_gpkg=False` (CSV-only mode):

```python
_, results_df = build_cfs_gpkg_from_rasters(
    input_folder="/path/to/rasters/",
    output_folder="/path/to/output/",
    layer_name="erosion_csv",
    master_gdf=my_master_gdf,
    master_key="ADM0_NAME",
    result_key="country",
    cf_name="Soil Erosion",
    cf_unit="t/ha/yr",
    area_type="country",
    write_gpkg=False,         # no GeoPackage, CSV only
    reset_gpkg=False,         # resume from checkpoints
)
```

Checkpoints are saved regardless of whether a GeoPackage is written.

---

## Checkpoint File Details

| Item | Details |
|------|---------|
| **Location** | `{output_folder}/_checkpoints_{layer_name}/` |
| **Format** | Apache Parquet (one file per raster) |
| **Naming** | `{flow_name}.parquet` |
| **Size** | ~50-100 KB each (tabular data only, no geometry) |
| **Integrity** | Atomic writes via `os.replace` prevent corruption |

## Performance Improvements

| Optimization | Impact |
|---|---|
| **Checkpoint resume** | Skips expensive per-raster computation entirely (~seconds to minutes saved per cached raster) |
| **Shapefile reuse** | Eliminates redundant shapefile load + CRS reprojection per raster call (~5-15s saved per raster for subcountry) |
| **CSV refactor** | Builds long-format CSV data directly from `flow_values` instead of re-merging from raw results |

## FAQ

**Q: Do I need to change my existing code?**
No. The default behavior (`reset_gpkg=True`) is fully backward-compatible. To use
crash recovery, the only change is passing `reset_gpkg=False` on re-runs.

**Q: Can I manually delete checkpoints?**
Yes. Delete the `_checkpoints_{layer_name}/` directory to force recomputation of
all rasters on the next run.

**Q: What if my input rasters changed since the last run?**
Use `reset_gpkg=True` to force a clean recomputation. Checkpoints assume the input
rasters have not changed between runs.

**Q: Do checkpoints accumulate across runs?**
Only when `reset_gpkg=False`. With the default `reset_gpkg=True`, the checkpoint
directory is wiped at the start of each run, then rebuilt as rasters are processed.
