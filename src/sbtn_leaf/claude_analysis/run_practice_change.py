"""Command-line entry point for the practice-change co-benefit analysis.

Run as::

    python -m sbtn_leaf.claude_analysis.run_practice_change [OUTDIR] [--no-figures] [--no-maps] [--no-download]

Writes tables, figures, maps and ``README.md`` (see
:func:`sbtn_leaf.claude_analysis.practice_change.run_practice_change_analysis`).
Maps fetch Natural Earth / RESOLVE Ecoregions geometry once (cached); pass
``--no-download`` to use only locally available geometry.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from sbtn_leaf.claude_analysis.practice_change import (
    default_output_dir,
    run_practice_change_analysis,
)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", nargs="?", default=None, help="Output directory.")
    parser.add_argument("--no-figures", action="store_true", help="Skip all figures.")
    parser.add_argument("--no-maps", action="store_true", help="Skip the choropleth maps only.")
    parser.add_argument("--no-download", action="store_true", help="Use only locally available geometry.")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir) if args.outdir else default_output_dir()
    result = run_practice_change_analysis(
        outdir=outdir,
        make_figures=not args.no_figures,
        make_maps=not args.no_maps,
        allow_download=not args.no_download,
    )
    print(f"Practice-change artifacts written to: {result['outdir']}  (maps rendered: {result['n_maps']})")
    print("\nExtent of benefits (ecoregion, median over regions):")
    cols = ["commodity", "switch", "n_regions", "med_d_soc", "med_d_se_red", "d_se_red_pct", "pct_win_win"]
    print(result["summary"][cols].round(2).to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
