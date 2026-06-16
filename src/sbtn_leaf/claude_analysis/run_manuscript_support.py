"""Command-line entry point for the manuscript-support statistics & figures.

Run as::

    python -m sbtn_leaf.claude_analysis.run_manuscript_support [OUTDIR] [--no-figures]

Writes the tables, figures and ``README.md`` (see
:func:`sbtn_leaf.claude_analysis.manuscript_support.run_manuscript_support`).
With no ``OUTDIR`` the default ``paper/claude_analysis/outputs/manuscript_support``
is used.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from sbtn_leaf.claude_analysis.manuscript_support import (
    default_output_dir,
    run_manuscript_support,
)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", nargs="?", default=None, help="Output directory.")
    parser.add_argument("--no-figures", action="store_true", help="Write only CSV tables + README.")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir) if args.outdir else default_output_dir()
    result = run_manuscript_support(outdir=outdir, make_figures=not args.no_figures)

    print(f"Manuscript-support artifacts written to: {result['outdir']}")
    print("\nBiome significance (median over focal flows; p = worst-case Kruskal–Wallis):")
    print(result["biome_summary"].round(4).to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
