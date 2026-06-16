"""Command-line entry point for the cross-indicator aggregation analysis.

Run as::

    python -m sbtn_leaf.claude_analysis.run_cross_indicator [OUTDIR] [--no-figures]

Writes the tables, figures and ``README.md`` (see
:func:`sbtn_leaf.claude_analysis.cross_indicator.run_cross_indicator_analysis`).
With no ``OUTDIR`` the default ``paper/claude_analysis/outputs`` is used.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from sbtn_leaf.claude_analysis.cross_indicator import (
    default_output_dir,
    run_cross_indicator_analysis,
)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "outdir",
        nargs="?",
        default=None,
        help="Output directory (default: paper/claude_analysis/outputs).",
    )
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Write only the CSV tables and README, skip the (regenerable) PNGs.",
    )
    args = parser.parse_args(argv)

    outdir = Path(args.outdir) if args.outdir else default_output_dir()
    result = run_cross_indicator_analysis(outdir=outdir, make_figures=not args.no_figures)

    sig = result["significance"]
    print(f"Cross-indicator analysis written to: {result['outdir']}")
    print("\nEcoregion significance (median across focal flows):")
    print(
        sig[["indicator_name", "eco_biome_eta2", "eco_realm_eta2", "subcty_within_country_frac"]]
        .round(3)
        .to_string(index=False)
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
