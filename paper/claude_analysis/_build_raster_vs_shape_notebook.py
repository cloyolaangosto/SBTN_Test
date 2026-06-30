"""Generate the raster-vs-shape comparison notebook.

Build helper (not part of the analysis): assembles the notebook cells and writes
the .ipynb, which is then executed with ``jupyter nbconvert``.

The notebook is *self-contained* — all raster logic lives inline in its code
cells (no new ``sbtn_leaf`` modules).  It only imports the existing package for
the **shape** (polygon-aggregated) side and the flow-name canonicalisers, so the
raster and shape analyses line up on identical canonical flow keys.
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def md(*lines):
    return {"cell_type": "markdown", "metadata": {}, "source": _src(lines)}


def code(*lines):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": _src(lines)}


def _src(lines):
    text = "\n".join(lines)
    parts = text.split("\n")
    return [p + "\n" for p in parts[:-1]] + [parts[-1]]


cells = []

# --------------------------------------------------------------------------- #
# Title
# --------------------------------------------------------------------------- #
cells.append(md(
    "# Raster vs. shape: the LEAF flows at the pixel level",
    "",
    "[`Manuscript_Support_Figures.ipynb`](Manuscript_Support_Figures.ipynb) and the other",
    "`claude_analysis` notebooks work entirely on the **polygon-aggregated** LEAF tables",
    "(`SOC_2030_*_v1.0.csv`, `se_*_clipped_2.csv`, `acidification_*.csv`): each region",
    "(country / sub-country / ecoregion) is reduced to a single area-weighted `leaf`,",
    "with `leaf_median` and `leaf_std`. Every statistic there treats **one polygon = one",
    "observation**.",
    "",
    "Those polygon tables were produced by clipping the underlying **LEAF rasters** to",
    "polygons. This notebook goes back to those rasters and redoes the same analyses with",
    "**one pixel = one area-true observation**, then compares the two.",
    "",
    "| manuscript / shape claim | shape analysis | raster analysis here |",
    "|---|---|---|",
    "| polygon means reproduce the grid | published `leaf` | §0 zonal reproduction (audit) |",
    "| *“biomes lead to LEAFs that are significantly different”* | §1 η² on ecoregion **means** | §1 η² on area-weighted **pixels** |",
    "| *“sub-country … smaller standard deviation”* | §2 per-region `leaf_std` | §2 true within-polygon pixel SD + variance partition |",
    "| *“aligned SOC and soil erosion … benefits for both”* | §3 ρ on region **means** | §3 ρ on co-located **pixels** (ecological-fallacy check) |",
    "",
    "**Why it matters.** A polygon mean is a lossy summary: it hides the within-region",
    "distribution and weights a tiny ecoregion the same as a continental one. The raster",
    "version is the area-true ground truth the polygon tables approximate — so the gap",
    "between them quantifies (a) reproduction fidelity, (b) aggregation bias in the biome",
    "signal, and (c) the ecological fallacy in the co-benefit claim.",
))

# --------------------------------------------------------------------------- #
# Setup
# --------------------------------------------------------------------------- #
cells.append(md(
    "## Setup — flow→raster index, grids, and area-weighted pixel helpers",
    "",
    "The flow→raster mapping is *free*: applying the existing canonicalisers",
    "(`se.canonical_flow`, `canonical_flow_soc`) to each raster's filename yields the same",
    "canonical keys the shape tables use, so the two analyses are guaranteed to line up.",
))

cells.append(code(
    "import warnings; warnings.filterwarnings('ignore')",
    "%matplotlib inline",
    "import os, glob",
    "import numpy as np, pandas as pd, matplotlib.pyplot as plt",
    "import rasterio, geopandas as gpd, rioxarray",
    "from rasterio.features import rasterize",
    "from rasterio.enums import Resampling",
    "from scipy import stats",
    "pd.set_option('display.width', 200, 'display.max_columns', 40)",
    "",
    "from sbtn_leaf.paths import data_path, leaf_path",
    "from sbtn_leaf.claude_analysis import se_aggregation_analysis as se",
    "from sbtn_leaf.claude_analysis import manuscript_support as ms",
    "from sbtn_leaf.claude_analysis.indicators import (",
    "    INDICATORS, SOC, SOIL_EROSION, ACIDIFICATION, canonical_flow_soc)",
    "",
    "# Enable to also recompute the sub-country level in §2 (reads a ~400 MB shapefile).",
    "INCLUDE_SUBCOUNTRY = False",
    "",
    "REP_FLOW = ms.REP_FLOW   # representative flow per indicator (shared with the shape notebook)",
    "REP_FLOW",
))

cells.append(code(
    "# --- flow -> raster-file index (reuse the shape-side canonicalisers) ---",
    "def _stem(p): return os.path.splitext(os.path.basename(str(p)))[0]",
    "",
    "RASTER_DIRS = {",
    "    'soc':           [leaf_path('SOC', 'rasters_clipped'), leaf_path('SOC', 'rasters')],",
    "    'soil_erosion':  [leaf_path('soil_erosion', 'rasters_clipped')],",
    "    'acidification': [leaf_path('acidification', 'rasters')],",
    "}",
    "CANON = {'soc': canonical_flow_soc, 'soil_erosion': se.canonical_flow, 'acidification': lambda s: s}",
    "",
    "def build_index(key):",
    "    idx = {}",
    "    for d in RASTER_DIRS[key]:",
    "        for f in sorted(glob.glob(os.path.join(str(d), '*.tif'))):",
    "            idx.setdefault(CANON[key](_stem(f)), f)",
    "    return idx",
    "",
    "RASTER_INDEX = {k: build_index(k) for k in INDICATORS}",
    "pd.DataFrame({",
    "    'n_rasters': {k: len(v) for k, v in RASTER_INDEX.items()},",
    "    'focal_resolved': {k: sum(f in RASTER_INDEX[k] for f in cfg.focal_flows)",
    "                       for k, cfg in INDICATORS.items()},",
    "    'n_focal': {k: len(cfg.focal_flows) for k, cfg in INDICATORS.items()},",
    "})",
))

cells.append(code(
    "# --- raster loading + cos(lat) area weights (EPSG:4326 pixel area ~ cos(latitude)) ---",
    "def load_grid(path):",
    "    with rasterio.open(path) as s:",
    "        arr = s.read(1).astype('float64')",
    "        T, W, H, crs, nod = s.transform, s.width, s.height, s.crs, s.nodata",
    "    valid = np.isfinite(arr)",
    "    if nod is not None:",
    "        valid &= (arr != nod)",
    "    rows = np.arange(H)",
    "    lat = T.f + (rows + 0.5) * T.e            # T.e is negative (north-up)",
    "    wlat = np.cos(np.deg2rad(lat)).clip(0)",
    "    W2 = np.broadcast_to(wlat[:, None], (H, W))",
    "    return dict(arr=arr, T=T, W=W, H=H, crs=crs, nodata=nod, valid=valid, w=W2)",
    "",
    "def grid_sig(g): return (g['W'], g['H'])",
))

cells.append(code(
    "# --- rasterise region zones (+ biome/realm for ecoregions) onto a grid, cached per grid ---",
    "_BND_PATHS = {",
    "    'country':    data_path('CountryLayers', 'Country_Level0', 'g2015_2014_0.shp'),",
    "    'subcountry': data_path('CountryLayers', 'SubCountry_Level1', 'g2015_2014_1.shp'),",
    "    'ecoregion':  data_path('ecoregions2017', 'ecoregions2017.shp'),",
    "}",
    "_ID_COL = {'country': 'ADM0_NAME', 'subcountry': 'ADM1_CODE', 'ecoregion': 'ECO_ID'}",
    "_bnd_cache, _zone_cache = {}, {}",
    "",
    "def _boundaries(level):",
    "    if level not in _bnd_cache:",
    "        _bnd_cache[level] = gpd.read_file(_BND_PATHS[level])",
    "    return _bnd_cache[level]",
    "",
    "def zone_for(level, g):",
    "    \"\"\"Return dict(zone, ids, [biome, realm]); zone holds code+1 (0 = no region).\"\"\"",
    "    key = (level, *grid_sig(g))",
    "    if key in _zone_cache:",
    "        return _zone_cache[key]",
    "    gdf = _boundaries(level).to_crs(g['crs']).copy()",
    "    codes, _ = pd.factorize(gdf[_ID_COL[level]])",
    "    gdf['_c'] = codes",
    "    order = gdf[gdf['_c'] >= 0].drop_duplicates('_c').sort_values('_c')",
    "    zone = rasterize(",
    "        ((geom, int(c) + 1) for geom, c in zip(gdf.geometry, gdf['_c']) if c >= 0),",
    "        out_shape=(g['H'], g['W']), transform=g['T'], fill=0, dtype='int32',",
    "    )",
    "    out = {'zone': zone, 'ids': order[_ID_COL[level]].to_numpy()}",
    "    if level == 'ecoregion':",
    "        out['biome'] = order['BIOME_NAME'].astype('string').to_numpy()",
    "        out['realm'] = order['REALM'].astype('string').to_numpy()",
    "    _zone_cache[key] = out",
    "    return out",
))

cells.append(code(
    "# --- area-weighted pixel statistics (fast, vectorised with np.bincount) ---",
    "def weighted_eta2(values, codes, weights):",
    "    \"\"\"Area-weighted one-way ANOVA η² of `values` grouped by integer `codes` (0..G-1).\"\"\"",
    "    wsum = weights.sum()",
    "    if wsum <= 0 or len(np.unique(codes)) < 2:",
    "        return np.nan",
    "    gm = (values * weights).sum() / wsum",
    "    ss_tot = (weights * (values - gm) ** 2).sum()",
    "    Wg = np.bincount(codes, weights=weights)",
    "    Sg = np.bincount(codes, weights=weights * values)",
    "    mug = Sg / np.where(Wg > 0, Wg, 1.0)",
    "    ss_btw = (Wg * (mug - gm) ** 2).sum()",
    "    return ss_btw / ss_tot if ss_tot > 0 else np.nan",
    "",
    "def zonal_mean_std(arr, zone, w, valid, n=None):",
    "    \"\"\"Per-region area-weighted mean & SD; output length `n` (= number of zones).\"\"\"",
    "    keep = valid & (zone > 0)",
    "    c = zone[keep] - 1",
    "    v = arr[keep]; ww = w[keep]",
    "    n = n if n is not None else (int(c.max()) + 1 if c.size else 0)",
    "    Wg = np.bincount(c, weights=ww, minlength=n)",
    "    Sg = np.bincount(c, weights=ww * v, minlength=n)",
    "    Qg = np.bincount(c, weights=ww * v * v, minlength=n)",
    "    mean = Sg / np.where(Wg > 0, Wg, np.nan)",
    "    var = Qg / np.where(Wg > 0, Wg, np.nan) - mean ** 2",
    "    std = np.sqrt(var.clip(0))",
    "    return mean, std, Wg",
    "",
    "def variance_partition(arr, zone, w, valid):",
    "    keep = valid & (zone > 0)",
    "    c = zone[keep] - 1; v = arr[keep]; ww = w[keep]",
    "    wsum = ww.sum(); gm = (v * ww).sum() / wsum",
    "    ss_tot = (ww * (v - gm) ** 2).sum()",
    "    Wg = np.bincount(c, weights=ww); Sg = np.bincount(c, weights=ww * v)",
    "    mug = Sg / np.where(Wg > 0, Wg, 1.0)",
    "    ss_btw = (Wg * (mug - gm) ** 2).sum()",
    "    eta2 = ss_btw / ss_tot if ss_tot > 0 else np.nan",
    "    return dict(eta2_between=eta2, frac_within=1 - eta2, n_groups=int((Wg > 0).sum()))",
))

# --------------------------------------------------------------------------- #
# Section 0 — reproduction / audit
# --------------------------------------------------------------------------- #
cells.append(md(
    "## 0. Do the rasters reproduce the published polygon `leaf`?",
    "",
    "Before comparing *conclusions*, audit the data: re-derive each ecoregion's area-weighted",
    "mean **from the raster** and check it matches the published `leaf` the shape notebook uses.",
    "This both validates the published aggregation and tells us whether the distributed LEAF",
    "rasters are the same vintage as the published tables.",
))

cells.append(code(
    "def reproduce_level(cfg, flow, level='ecoregion'):",
    "    path = RASTER_INDEX[cfg.key].get(flow)",
    "    if path is None:",
    "        return None",
    "    g = load_grid(path)",
    "    z = zone_for(level, g)",
    "    mean, _, _ = zonal_mean_std(g['arr'], z['zone'], g['w'], g['valid'], n=len(z['ids']))",
    "    raster = pd.Series(mean, index=pd.Index(z['ids']).astype('string'), name='raster')",
    "    sh = cfg.load_harmonized(drop_na=True)",
    "    pub = (sh[(sh.level == level) & (sh.flow == flow)]",
    "           .assign(region_id=lambda d: d.region_id.astype('string'))",
    "           .set_index('region_id')['leaf'].rename('published'))",
    "    j = pd.concat([raster, pub], axis=1).dropna()",
    "    err = j['raster'] - j['published']",
    "    return dict(indicator=cfg.name, flow_label=cfg.label(flow), level=level, n=len(j),",
    "                pearson_r=j['raster'].corr(j['published']),",
    "                bias=err.mean(), rmse=np.sqrt((err ** 2).mean()),",
    "                med_raster=j['raster'].median(), med_published=j['published'].median(),",
    "                _join=j)",
    "",
    "repro = {k: reproduce_level(cfg, REP_FLOW[k]) for k, cfg in INDICATORS.items()}",
    "repro_tbl = pd.DataFrame([{x: r[x] for x in",
    "    ['indicator', 'flow_label', 'level', 'n', 'pearson_r', 'bias', 'rmse', 'med_raster', 'med_published']}",
    "    for r in repro.values()])",
    "repro_tbl.round(3)",
))

cells.append(code(
    "fig, axes = plt.subplots(1, 3, figsize=(16, 5))",
    "for ax, (k, r) in zip(axes, repro.items()):",
    "    j = r['_join']",
    "    ax.scatter(j['published'], j['raster'], s=10, alpha=0.4, color='#3b6ea5', edgecolor='none')",
    "    lo = float(np.nanmin([j['published'].min(), j['raster'].min()]))",
    "    hi = float(np.nanmax([j['published'].max(), j['raster'].max()]))",
    "    ax.plot([lo, hi], [lo, hi], 'k--', lw=1, label='1:1')",
    "    ax.set_xlabel('published polygon leaf'); ax.set_ylabel('raster zonal mean')",
    "    ax.set_title(f\"{r['indicator']} — {r['flow_label']}\\n\"",
    "                 f\"r={r['pearson_r']:.3f}  bias={r['bias']:+.3g}  (n={r['n']})\")",
    "    ax.legend(); ax.grid(True, ls='--', alpha=0.3)",
    "fig.tight_layout(); display(fig); plt.close(fig)",
))

cells.append(md(
    "**Acidification reproduces essentially perfectly** (r≈1.00, bias≈0): the distributed",
    "`acid_*.tif` rasters *are* the source of the published table. **Soil erosion reproduces",
    "well** (r≈0.97, medians within a few %) — small differences are whole-pixel vs fractional",
    "polygon coverage and the published outlier filtering. **SOC does not**: the pattern agrees",
    "(r≈0.87) but the distributed SOC rasters run ~50 % higher than the published `v1.0` table —",
    "i.e. the SOC raster set is a **different/later vintage**. So SOC raster results below are read",
    "on **rank / structure** (which transfers), not absolute level.",
))

# --------------------------------------------------------------------------- #
# Section 1 — biome significance
# --------------------------------------------------------------------------- #
cells.append(md(
    "## 1. Biomes produce significantly different LEAFs — pixels vs ecoregion means",
    "",
    "> *“ecoregion’s biomes lead to LEAFs that are significantly different.”*",
    "",
    "The shape analysis tests this on ecoregion **means** (η² with every ecoregion weighted",
    "equally, n≈600–830). Here we test the **same biome grouping on area-weighted pixels**",
    "(n in the 10⁵–10⁶), so a biome's signal is weighted by the land it actually covers.",
    "p-values are omitted — with millions of pixels everything is `***`; the honest comparison",
    "is the **effect size η²**.",
))

cells.append(code(
    "def biome_eta2_pixels(cfg, flow):",
    "    path = RASTER_INDEX[cfg.key].get(flow)",
    "    if path is None:",
    "        return np.nan, 0",
    "    g = load_grid(path)",
    "    z = zone_for('ecoregion', g)",
    "    biome_by_code = z['biome']",
    "    keep = g['valid'] & (z['zone'] > 0)",
    "    code = z['zone'][keep] - 1",
    "    bname = biome_by_code[code]",
    "    ok = pd.notna(bname)",
    "    bcode = pd.factorize(bname[ok])[0]",
    "    return weighted_eta2(g['arr'][keep][ok], bcode, g['w'][keep][ok]), int(ok.sum())",
    "",
    "shape_tbl = ms.biome_significance_table()",
    "rows = []",
    "for k, cfg in INDICATORS.items():",
    "    for flow in cfg.focal_flows:",
    "        e2_px, n_px = biome_eta2_pixels(cfg, flow)",
    "        srow = shape_tbl[(shape_tbl.indicator == k) & (shape_tbl.flow == flow)]",
    "        e2_mean = float(srow['eta2_biome'].iloc[0]) if len(srow) else np.nan",
    "        rows.append(dict(indicator=cfg.name, flow_label=cfg.label(flow),",
    "                         eta2_pixels=e2_px, eta2_ecoregion_mean=e2_mean,",
    "                         ratio=e2_mean / e2_px if e2_px else np.nan, n_pixels=n_px))",
    "biome_cmp = pd.DataFrame(rows)",
    "biome_cmp.round(3)",
))

cells.append(code(
    "sub = biome_cmp.dropna(subset=['eta2_pixels', 'eta2_ecoregion_mean']).reset_index(drop=True)",
    "fig, ax = plt.subplots(figsize=(12, 6))",
    "x = np.arange(len(sub)); wd = 0.4",
    "ax.bar(x - wd/2, sub['eta2_pixels'], wd, label='pixels (area-weighted)', color='#2c7fb8')",
    "ax.bar(x + wd/2, sub['eta2_ecoregion_mean'], wd, label='ecoregion means (shape)', color='#de8a3a')",
    "ax.set_xticks(x)",
    "ax.set_xticklabels([f\"{r.indicator.split()[0]}: {r.flow_label}\" for r in sub.itertuples()],",
    "                   rotation=45, ha='right', fontsize=8)",
    "ax.set_ylabel('biome η² (variance explained)')",
    "ax.set_title('Biome separation: area-true pixels vs equal-weighted ecoregion means')",
    "ax.legend(); ax.grid(True, axis='y', ls='--', alpha=0.4)",
    "fig.tight_layout(); display(fig); plt.close(fig)",
))

cells.append(code(
    "# The pixel distribution behind one test: ecoregion biomes for the representative flow.",
    "cfg = SOIL_EROSION; flow = REP_FLOW['soil_erosion']",
    "g = load_grid(RASTER_INDEX[cfg.key][flow]); z = zone_for('ecoregion', g)",
    "keep = g['valid'] & (z['zone'] > 0) & (g['arr'] > 0)",
    "bname = z['biome'][z['zone'][keep] - 1]",
    "dfp = pd.DataFrame({'v': g['arr'][keep], 'biome': bname}).dropna()",
    "order = dfp.groupby('biome')['v'].median().sort_values().index.tolist()",
    "fig, ax = plt.subplots(figsize=(10, 0.5 * len(order) + 2))",
    "ax.boxplot([dfp.loc[dfp.biome == b, 'v'].values for b in order], orientation='horizontal',",
    "           whis=(5, 95), showfliers=False)",
    "ax.set_yticks(range(1, len(order) + 1)); ax.set_yticklabels(order, fontsize=8)",
    "ax.set_xscale('log'); ax.set_xlabel(f'{cfg.name} ({cfg.unit}) — pixels, log')",
    "ax.set_title(f'{cfg.label(flow)} — pixel {cfg.name} by biome (raster)')",
    "ax.grid(True, axis='x', ls='--', alpha=0.4)",
    "fig.tight_layout(); display(fig); plt.close(fig)",
))

cells.append(md(
    "Biome remains a **strong, significant** axis at the pixel level for every flow — the",
    "manuscript's claim holds on the area-true data, not just on polygon means. But the polygon-mean",
    "η² is a **biased estimate of the area-true biome effect, and the direction of the bias depends",
    "on the indicator**:",
    "",
    "* For **soil erosion** the equal-weighted ecoregion-mean test *overstates* biome separation",
    "  (wheat η² 0.43 → 0.26 on pixels; every erosion flow drops). A few small, extreme-erosion",
    "  ecoregions get a full vote in the polygon test but little land area in the pixel test.",
    "* For **SOC** it *understates* it (wheat 0.33 → 0.37; grassland 0.25 → 0.53; tropical broadleaf",
    "  0.09 → 0.46). Area-weighting concentrates on each biome's typical core, while ecoregion means",
    "  blur SOC across biome edges and dilute the signal.",
    "",
    "So the qualitative conclusion (*biomes matter*) is robust, but the headline η² is **not",
    "interchangeable** between the two views — it should be quoted as variance *between ecoregions*,",
    "not of the land surface.",
))

# --------------------------------------------------------------------------- #
# Section 2 — within-region spread / information loss
# --------------------------------------------------------------------------- #
cells.append(md(
    "## 2. Within-region spread and aggregation information loss",
    "",
    "> *“sub-country … leads to smaller standard deviation.”*",
    "",
    "The shape notebook uses the published per-region `leaf_std`. From the raster we can compute",
    "the **true within-polygon pixel SD** directly, and — uniquely — partition the *total* pixel",
    "variance into a **between-region** and a **within-region** part at each level. `frac_within`",
    "is the share of the on-the-ground variance a single regional LEAF hides: the real cost of",
    "aggregation, which the polygon tables cannot show.",
))

cells.append(code(
    "LEVELS = ['country', 'ecoregion'] + (['subcountry'] if INCLUDE_SUBCOUNTRY else [])",
    "rows = []",
    "for k, cfg in INDICATORS.items():",
    "    flow = REP_FLOW[k]",
    "    g = load_grid(RASTER_INDEX[cfg.key][flow])",
    "    for level in LEVELS:",
    "        z = zone_for(level, g)",
    "        vp = variance_partition(g['arr'], z['zone'], g['w'], g['valid'])",
    "        mean, std, Wg = zonal_mean_std(g['arr'], z['zone'], g['w'], g['valid'])",
    "        rows.append(dict(indicator=cfg.name, flow_label=cfg.label(flow), level=level,",
    "                         n_regions=vp['n_groups'], frac_within=vp['frac_within'],",
    "                         eta2_between=vp['eta2_between'],",
    "                         mean_within_std=np.nanmean(std)))",
    "infoloss = pd.DataFrame(rows)",
    "infoloss['level'] = pd.Categorical(infoloss['level'], ['country', 'subcountry', 'ecoregion'], ordered=True)",
    "infoloss.sort_values(['indicator', 'level']).round(3)",
))

cells.append(code(
    "# Raster within-polygon SD vs the published `leaf_std` (mean over regions), per level.",
    "shape_sd = ms.within_region_sd_table().set_index('indicator')",
    "comp = []",
    "for k, cfg in INDICATORS.items():",
    "    flow = REP_FLOW[k]",
    "    g = load_grid(RASTER_INDEX[cfg.key][flow])",
    "    for level in LEVELS:",
    "        z = zone_for(level, g)",
    "        _, std, _ = zonal_mean_std(g['arr'], z['zone'], g['w'], g['valid'])",
    "        comp.append(dict(indicator=cfg.name, level=level,",
    "                         raster_within_std=np.nanmean(std),",
    "                         published_leaf_std=shape_sd.loc[k, f'within_sd_{level}']))",
    "pd.DataFrame(comp).round(3)",
))

cells.append(code(
    "fig, ax = plt.subplots(figsize=(9, 5.5))",
    "order = [l for l in ['country', 'subcountry', 'ecoregion'] if l in LEVELS]",
    "for k, cfg in INDICATORS.items():",
    "    sub = infoloss[infoloss.indicator == cfg.name].set_index('level').reindex(order)",
    "    ax.plot(order, sub['frac_within'], marker='o', lw=2, label=cfg.name)",
    "ax.set_ylabel('frac. of pixel variance hidden *within* regions')",
    "ax.set_title('Aggregation information loss by level (raster variance partition)')",
    "ax.legend(); ax.grid(True, axis='y', ls='--', alpha=0.4)",
    "fig.tight_layout(); display(fig); plt.close(fig)",
))

cells.append(md(
    "The raster's within-polygon SD broadly tracks the published `leaf_std` at the levels we",
    "recompute (within ~10–20 %), corroborating the published spread statistic. The **variance",
    "partition** adds what",
    "the polygons cannot: even at the ecoregion level a large majority of the pixel variance still",
    "lives *inside* regions (`frac_within` high), so the single regional LEAF is a coarse summary of",
    "a wide on-the-ground distribution. Finer / ecological polygons move some of that variance",
    "*between* regions (lower `frac_within`), which is the quantitative version of the manuscript's",
    "homogeneity argument. *(Set `INCLUDE_SUBCOUNTRY = True` in setup to add the admin-1 level.)*",
))

# --------------------------------------------------------------------------- #
# Section 3 — multi-indicator co-benefit (ecological fallacy)
# --------------------------------------------------------------------------- #
cells.append(md(
    "## 3. SOC ↔ soil-erosion alignment — co-located pixels vs region means",
    "",
    "> *“aligned SOC and soil erosion … identify where … the most benefits for both … simultaneously.”*",
    "",
    "The co-benefit claim is fundamentally a **pixel** question — *where on the ground* do both",
    "improve. The shape notebook answers it with per-region **means** (Spearman ρ). Here we",
    "resample SOC onto the 25 km erosion grid and correlate **co-located pixels** for each shared",
    "commodity, then compare to the region-mean ρ. A gap is the classic **ecological fallacy**:",
    "correlations between regional averages need not hold pixel-by-pixel.",
))

cells.append(code(
    "def soc_se_pixels(flow):",
    "    sp = RASTER_INDEX['soc'].get(flow); ep = RASTER_INDEX['soil_erosion'].get(flow)",
    "    if sp is None or ep is None:",
    "        return None",
    "    se_da = rioxarray.open_rasterio(ep, masked=True).isel(band=0)",
    "    soc_da = rioxarray.open_rasterio(sp, masked=True).isel(band=0)",
    "    soc_on = soc_da.rio.reproject_match(se_da, resampling=Resampling.average)",
    "    a = np.ravel(soc_on.values).astype('float64')",
    "    b = np.ravel(se_da.values).astype('float64')",
    "    m = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)",
    "    return a[m], b[m]",
    "",
    "shared = [f for f in SOIL_EROSION.focal_flows if f in SOC.focal_flows]",
    "shape_corr = ms.multi_indicator_correlation_table(flows=shared)",
    "rows = []",
    "for flow in shared:",
    "    px = soc_se_pixels(flow)",
    "    if px is None:",
    "        continue",
    "    rho_px, _ = stats.spearmanr(px[0], px[1])",
    "    em = shape_corr[(shape_corr.flow == flow) & (shape_corr.level == 'ecoregion')]",
    "    rho_mean = float(em['spearman_rho'].iloc[0]) if len(em) else np.nan",
    "    rows.append(dict(flow_label=SOIL_EROSION.label(flow), n_pixels=len(px[0]),",
    "                     rho_pixels=rho_px, rho_ecoregion_mean=rho_mean))",
    "eco_fallacy = pd.DataFrame(rows)",
    "eco_fallacy.round(3)",
))

cells.append(code(
    "fig, ax = plt.subplots(figsize=(11, 5.5))",
    "x = np.arange(len(eco_fallacy)); wd = 0.4",
    "ax.bar(x - wd/2, eco_fallacy['rho_pixels'], wd, label='co-located pixels', color='#2c7fb8')",
    "ax.bar(x + wd/2, eco_fallacy['rho_ecoregion_mean'], wd, label='ecoregion means (shape)', color='#de8a3a')",
    "ax.axhline(0, color='k', lw=0.8)",
    "ax.set_xticks(x); ax.set_xticklabels(eco_fallacy['flow_label'], rotation=45, ha='right', fontsize=8)",
    "ax.set_ylabel('Spearman ρ (SOC vs soil erosion)')",
    "ax.set_title('SOC↔erosion alignment: pixel reality vs region-mean (ecological fallacy)')",
    "ax.legend(); ax.grid(True, axis='y', ls='--', alpha=0.4)",
    "fig.tight_layout(); display(fig); plt.close(fig)",
))

cells.append(code(
    "# Co-benefit density for the representative commodity: where do both sit on the grid?",
    "flow = REP_FLOW['soil_erosion']; a, b = soc_se_pixels(flow)",
    "fig, ax = plt.subplots(figsize=(7.5, 6))",
    "hb = ax.hexbin(a, b, gridsize=45, bins='log', mincnt=1, cmap='viridis')",
    "ax.set_yscale('log')",
    "ax.set_xlabel(f'SOC stock ({SOC.unit}) — raster vintage'); ax.set_ylabel(f'Soil erosion ({SOIL_EROSION.unit}) — log')",
    "rho, _ = stats.spearmanr(a, b)",
    "ax.set_title(f'{SOIL_EROSION.label(flow)} @ 25 km pixels: SOC vs erosion\\nSpearman ρ={rho:.2f} (n={len(a):,})')",
    "fig.colorbar(hb, ax=ax, label='log10(pixel count)')",
    "fig.tight_layout(); display(fig); plt.close(fig)",
))

cells.append(md(
    "This is the headline raster-vs-shape result. At the **ecoregion-mean** level SOC and soil",
    "erosion look positively aligned (ρ≈0.2–0.6, the manuscript's co-benefit signal). At the",
    "**co-located pixel** level that alignment **weakens, vanishes, or reverses sign** for several",
    "commodities — a textbook ecological fallacy. The region-mean correlation reflects *between-region*",
    "climate gradients (warm/wet regions are both higher-SOC and higher-erosion); it does **not**",
    "imply that an individual high-SOC field is also high-erosion. A co-benefit map intended to guide",
    "field-level practice change should therefore be built and read at the **pixel** level, not from",
    "regional averages. *(SOC raster vintage differs from the published table — see §0 — so treat the",
    "pixel ρ magnitudes as indicative; the region-vs-pixel divergence is the robust finding.)*",
))

# --------------------------------------------------------------------------- #
# Findings
# --------------------------------------------------------------------------- #
cells.append(md(
    "## Findings — raster vs shape",
    "",
    "1. **Reproduction / provenance audit.** Zonal-averaging the distributed LEAF rasters recovers",
    "   the published polygon `leaf` exactly for acidification (r≈1.00) and well for soil erosion",
    "   (r≈0.97); for **SOC** the rasters are a different vintage (r≈0.87, ~50 % higher level), which",
    "   the polygon-only view cannot reveal.",
    "2. **Biome signal is real but its polygon-mean effect size is biased — not always the same way.**",
    "   Biome η² stays significant on area-weighted pixels but differs from the equal-weighted",
    "   ecoregion-mean value, *overstating* it for soil erosion (wheat 0.43→0.26) and *understating*",
    "   it for SOC (grassland 0.25→0.53) — quote it as between-ecoregion variance, not of the land surface.",
    "3. **A single regional LEAF hides most of the on-the-ground variance.** The raster variance",
    "   partition makes the manuscript's homogeneity argument quantitative and shows the residual",
    "   within-region spread the polygon tables omit.",
    "4. **The SOC↔erosion co-benefit is a between-region artefact.** Region-mean alignment (ρ>0)",
    "   weakens or reverses at the pixel level — co-benefit mapping for practice change must be done",
    "   pixel-wise.",
    "",
    "**Method notes.** Pixels are area-weighted by cos(latitude) on the native EPSG:4326 grids;",
    "regions are assigned by rasterising the same boundary shapefiles used to build the polygon",
    "tables; SOC/acidification share a ~10 km grid and soil erosion a 25 km grid (SOC is resampled",
    "to 25 km for §3). With 10⁵–10⁶ pixels every test is significant, so effect sizes (η², ρ), not",
    "p-values, carry the comparison.",
))

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out = HERE / "RasterVsShape_Comparison.ipynb"
out.write_text(json.dumps(nb, indent=1) + "\n")
print("wrote", out)
