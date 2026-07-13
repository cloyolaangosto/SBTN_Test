# Creating SBTN Land LEAFs: Challenges and Learnings Across Three Soil-Quality Indicators

*A retrospective on the design, data, and engineering decisions behind the
Soil Organic Carbon, Soil Erosion, and Terrestrial Acidification Land
Environmental Assessment Factors (LEAFs).*

---

## Executive Summary

Land Environmental Assessment Factors (LEAFs) are default, location-specific
factors that let organizations compare their land use and management practices
against ecoregional thresholds and set Science Based Targets Network (SBTN) Land
v2 targets. This repository delivers LEAFs for three soil-quality indicators —
**Soil Organic Carbon (SOC)**, **Soil Erosion**, and **Terrestrial
Acidification** — together with the open data, documentation, and reproducible
workflows needed to extend them.

Building these three indicators surfaced a consistent set of lessons that are
worth carrying into future SBTN work:

1. **Match modeling ambition to data availability and decision relevance.** The
   three indicators required very different levels of effort — from a full
   dynamic simulation to simply adopting and correcting published factors — and
   resisting the urge to over-engineer the simpler ones was a strategic win.
2. **The dominant challenge was data harmonization, not modeling.** No single
   global dataset is complete or shares a common resolution, so most of the
   effort went into stitching fragmented sources onto a single reference grid
   and gracefully filling gaps.
3. **Translating plot-scale science to a global grid is iterative.** Several
   scientific inputs had to be revisited and recomputed mid-project; the work is
   never "run once."
4. **A research pipeline at global scale must be engineered like production
   software** — with fault tolerance, testing, and validation across every
   geographic scale — because multi-hour runs and scale-dependent bugs are the
   norm.
5. **The "last mile" — aggregating raw maps into clean, decision-ready factors —
   is as important as the modeling itself** for the LEAFs to be usable and
   trustworthy.

The remainder of this report develops these themes and closes with consolidated
recommendations and an honest account of the work that remains.

---

## Background: What We Built and Why

A LEAF answers a deceptively simple question for a company sourcing a commodity:
*given where this is grown and how it is managed, what is the expected impact on
this soil-quality indicator, and how does that compare to a safe operating
threshold?* To answer it consistently for any commodity, anywhere in the world,
each indicator passes through four broad stages:

1. **Gather** global input data (climate, soils, crop yields, land use).
2. **Harmonize** every input onto a common spatial grid and projection.
3. **Model or compute** the indicator value for each location and management
   option.
4. **Aggregate** the resulting maps into representative factors for ecoregions,
   countries, and subcountry regions that companies can actually look up.

The three indicators differ most in stage 3:

| Indicator | What it measures | Unit | Underlying method | Native resolution |
|---|---|---|---|---|
| **Soil Organic Carbon** | Estimated SOC stock on the land | t C/ha | RothC dynamic soil-carbon model, projected to 2030 from a 2016 baseline | ~9 km (1/12°) |
| **Soil Erosion** | Annual soil loss | t soil/ha/yr | RUSLE (Revised Universal Soil Loss Equation) | 25 km |
| **Terrestrial Acidification** | Soil acidification potential of emitted gases | kg SO₂-eq./kg | Published characterization factors (Roy et al., 2014) | ~2° × 2.5° |

All three were ultimately delivered for 42 land-use classes (28 agricultural,
15 forest, 1 grassland) and aggregated to ecoregion, country, and subcountry
levels.

---

## A Spectrum of Modeling Ambition

The single most useful strategic learning was that **the three indicators sit on
a spectrum of modeling effort, and that this is a feature, not an
inconsistency.**

- **Soil Organic Carbon** sits at the demanding end. It uses RothC, a *dynamic
  process model* that simulates monthly carbon turnover across five soil pools,
  driven by temperature, moisture, soil texture, plant cover, and organic-matter
  inputs. Producing a single SOC LEAF means assembling a dozen distinct inputs
  and running a time-stepped simulation forward to 2030.

- **Soil Erosion** sits in the middle. RUSLE reduces erosion to a *product of
  five factors* (rainfall erosivity, soil erodibility, slope, ground cover, and
  protection). Because three of those factors depend only on soil and weather —
  not on the crop — they could be pre-combined into a single reusable base
  layer, leaving only a cover-and-management factor to vary by commodity. This
  turned a potentially heavy computation into a lightweight multiplication.

- **Terrestrial Acidification** sits at the light end. Robust, peer-reviewed
  characterization factors already existed (Roy et al., 2014, as used in the
  IMPACT World+ method), so *no new modeling was performed at all*. The work was
  limited to correcting calculation errors in the source, normalizing the
  factors to a common SO₂-equivalent basis, and re-gridding them.

The learning: **invest scientific and computational effort where it changes
decisions and where good factors do not already exist.** Spending the same
modeling budget on acidification that SOC required would have added cost without
adding accuracy or insight. Conversely, the deliberate pre-computation of the
erosion base layer shows how a small amount of up-front design dramatically
simplified everything downstream.

---

## Challenge 1 — Stitching Together a Fragmented Global Data Landscape

By far the largest share of effort went not into the science of any single model
but into **making fundamentally incompatible datasets work together.** The
necessary inputs came from many institutions, each with its own resolution,
projection, and coverage:

| Input | Source | Native resolution |
|---|---|---|
| Precipitation | NASA GPM (IMERG) | ~0.1° |
| Temperature | NASA GLDAS | ~1.0° |
| Soil carbon, clay, sand | ISRIC SoilGrids | ~250 m |
| Erosion base factors | GloSEM (JRC) | 25 km |
| Acidification factors | Roy et al. (2014) | ~2° × 2.5° |
| Crop yields | FAOSTAT (country) + SPAM (gridded) | mixed |

Reconciling resolutions that span four orders of magnitude was a recurring tax
on the project. Two design choices proved decisive:

**A single reference grid as the backbone.** Every input was resampled onto a
common base — the UHTH ecoregion zoning derived from earlier published work
(Morais, Teixeira & Domingos, 2019, following FAO GAEZ suitability maps). Having
one canonical grid that all indicators share meant harmonization logic could be
written once and reused, and that results across indicators remain spatially
comparable.

**A layered "best-available" strategy for crop yields.** Yields are a critical
SOC input (they drive how much plant residue returns carbon to the soil), yet no
single dataset covers every crop in every location. The solution was a
prioritized fallback: use gridded SPAM yields where available; otherwise fall
back to country-level FAOSTAT averages adjusted for irrigation practice; then to
ecoregion averages; then to broader biome averages; and finally to a
nearest-neighbor fill. This cascade guarantees a complete global surface while
always preferring the most specific data available — a pragmatic answer to
pervasive data gaps that recurs throughout environmental modeling.

**Provenance and openness tension.** Some inputs (notably the acidification
factors) were obtained directly from the original authors rather than a public
portal, and the largest datasets — global rasters and intermediate outputs — are
simply too big to live inside a public Git repository. This created an ongoing
tension between full reproducibility and practical distribution that the project
is still resolving (several "where will we host this data" questions remain open
in the documentation). The learning is to **budget for data hosting and
provenance from the start**, treating them as first-class deliverables rather
than afterthoughts.

---

## Challenge 2 — Translating Scientific Models to a Global, Operational Scale

The published science behind these indicators was largely developed and
validated at the scale of a plot, a field, or a research station. Operationalizing
it for *every commodity, everywhere* introduced challenges that the original
science never had to confront.

**From a point to a planet.** Models like RothC and the Thornthwaite
evapotranspiration method were conceived as single-location calculations. They
had to be re-expressed to run across global raster grids — for example,
evapotranspiration depends on day length, which varies continuously with
latitude and month and therefore had to be computed cell by cell across the
entire map. The project ultimately maintains both a simple single-point version
(useful for understanding one farm) and a full gridded version (for global
coverage), which is itself a useful pattern for balancing interpretability
against scale.

**A combinatorial explosion of management options.** A central purpose of LEAFs
is to let companies compare *practices*, not just locations — rainfed versus
irrigated, conventional versus conservation tillage, residues left on or removed
from the field, with or without cover crops. These options multiply: for cereals
alone, SOC required up to eight distinct scenario combinations per crop, and
each new option (such as conservation tillage, which had to be adapted from a
separate published method) multiplies the number of model runs and outputs to
manage. Designing the pipeline so that management options compose cleanly — and
so that a new option does not require rewriting the core model — was essential to
keeping the work tractable.

**Science is iterative, even after "completion."** Perhaps the most honest
learning is that several methodological inputs had to be discovered to be wrong,
corrected, and re-run well into the project. The evapotranspiration inputs were
recomputed and replaced; the timing of how crop residues are returned to the
soil over the year was revised; the definition of which months count as
"soil-covered" was corrected; and at one point the land-management cover factors
were found not to be flowing correctly into the SOC results, forcing a
regeneration of outputs. None of these were failures of competence — they are
the normal texture of building a complex scientific pipeline. The practical
consequence is that **the system must be built to be re-run cheaply and often**,
because "run once and publish" is not how this work actually proceeds.

---

## Challenge 3 — Engineering for Scale, Reliability, and Reproducibility

Producing LEAFs globally is a genuinely large computation, and treating the code
as throwaway research scripts would not have survived contact with that scale.
Several engineering lessons emerged.

**Long runs demand fault tolerance.** Generating factors across thousands of
subcountry regions and dozens of rasters can take hours. A single crash, timeout,
or power interruption partway through was expensive enough that the team built a
**checkpoint-and-resume system**: each unit of work is saved as it completes, so
an interrupted run can pick up where it left off instead of starting over.
Checkpoints are written atomically so an interruption can never leave a corrupted
result. This fault-tolerance work was not "extra" — at this scale it was the
difference between a pipeline that finishes and one that perpetually restarts.

**Bugs hide at scale.** A recurring and somewhat counterintuitive lesson was that
code which worked perfectly at one geographic level failed at another — logic
validated on a few hundred countries broke on thousands of ecoregions or
subcountry units, where memory and performance behave differently. On more than
one occasion the team deliberately **reverted a faster implementation in favor of
a slower one that was correct at every scale.** The takeaway is to validate
across the full range of real inputs, not a convenient subset, and to treat
correctness as non-negotiable relative to speed.

**Reproducibility involves real trade-offs.** Beyond the data-hosting issue
already noted, the project repeatedly chose to *lock* input dataset versions —
keeping an older version of a reference dataset even after a newer one was
released — specifically so that new outputs stay consistent and comparable with
the established methodology. Reproducibility is not free: it sometimes means
forgoing the latest data to preserve a coherent baseline, and that trade-off
should be made deliberately and documented.

**Testing as a safety net.** The presence of an automated test suite covering the
core scientific calculations and the parallel processing paths reflects a mature
recognition that, in a pipeline this interconnected, a small change in one place
can silently corrupt distant results. Tests are how a project of this complexity
keeps iterating without regressing.

---

## Challenge 4 — From Raw Maps to Decision-Ready Factors

A global raster of model outputs is not yet something a company can use. The
final stage — turning pixels into a small, trustworthy set of lookup values — carried
its own challenges and lessons.

**Aggregation and representativeness.** Results are summarized to ecoregion,
country, and subcountry levels. A clear lesson, surfaced directly in the guidance
to users, is that **the coarsest aggregation is the least trustworthy**: country
averages can blend wildly different growing conditions and are recommended only
as a last resort, with ecoregion-level factors being far more representative.
Communicating these representativeness caveats is part of delivering the factors
responsibly.

**Outlier handling.** A handful of extreme pixels — often artifacts of sparse
input data — can badly distort an aggregated factor. Considerable effort went
into filtering strategies so that anomalous values do not propagate into the
published numbers, while genuine variation is preserved. Notably, this was
*not* applied uniformly: acidification, whose factors are well-behaved across
their coverage, needed no such filtering. Knowing when robust statistics are
necessary and when they are unnecessary overhead is itself a useful judgment.

**The last mile matters.** It would be easy to view aggregation and quality
control as clean-up after the "real" modeling. In practice this stage determines
whether the LEAFs are usable, comparable, and credible — and it absorbed a
meaningful share of the project's debugging effort (including a late discovery
that erosion averages had been computed incorrectly and had to be regenerated).
The last mile deserves the same rigor as the science.

---

## Key Learnings and Recommendations

**Strategy**
- Calibrate modeling effort to data availability and decision relevance; adopt
  and correct existing factors where they are sound rather than rebuilding them.
- A small amount of up-front design (e.g., pre-combining invariant factors into a
  reusable base layer) pays large dividends downstream.

**Data**
- Standardize on a single reference grid early and write harmonization logic
  once; it is the backbone everything else depends on.
- Expect incomplete data and design explicit, prioritized fallback strategies
  rather than ad-hoc patches.
- Treat data hosting and provenance as first-class deliverables, planned from the
  outset.

**Modeling**
- Build for re-runnability: methodological inputs *will* need correction and
  recomputation, so make that cheap.
- Design management options to compose cleanly so new practices do not require
  reworking the core models.

**Engineering**
- Engineer global pipelines like production software: checkpointing, atomic
  writes, automated tests, and validation across *every* geographic scale.
- Prefer slower-but-correct over faster-but-fragile, and make version-locking
  decisions deliberately and transparently.

**Communication**
- Ship representativeness and uncertainty guidance alongside the factors; the
  appropriate level of aggregation is part of the product.

**Honest account of remaining work.** This is a living project, and several items
are openly outstanding: the documentation still contains placeholders and "to do"
notes (including unconfirmed resolutions and citations, and the public-data
hosting locations); managed grassland is not yet fully implemented; and the
aggregation methodology sections are not finished. None of these undermine the
delivered LEAFs, but closing them — and completing the public-data and
documentation story — is the most valuable next step to make the work fully
transparent and reproducible for the organizations that will rely on it.

---

*Sources: this report synthesizes the repository's indicator documentation (SOC,
Soil Erosion, and Terrestrial Acidification), the crash-recovery engineering
guide, the project's open task notes, the LEAF output guidance, and the
implementation history of the processing code.*
