#!/usr/bin/env python3
"""Draw the branch / merge topology of this repository as a static image.

Reconstructs, from the git history, the main integration line (the ``paper``
branch) and every feature branch that was merged into it, then renders a clean
swim-lane "git graph" to:

    documentation/figures/branch_topology.png

Each merged feature branch is shown as an arc that forks off the main line at
its branch point and merges back at the merge commit, labelled with the branch
name and its lifespan. Everything is derived from ``git`` at run time, so the
picture stays correct as history grows.

Usage:
    python scripts/branch_graph.py
"""
from __future__ import annotations

import re
import subprocess
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT = REPO_ROOT / "documentation" / "figures" / "branch_topology.png"

MAIN_COLOR = "#3A3A3A"
BRANCH_COLORS = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756",
    "#72B7B2", "#B279A2", "#EECA3B", "#9D755D",
]


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL
    ).strip()


def dt(iso: str) -> datetime:
    return datetime.fromisoformat(iso)


def branch_name(subject: str) -> str:
    """Pull a readable branch name out of a merge commit subject."""
    m = re.search(r"(claude/[\w./-]+?)(?:-[a-z0-9]{6})?(?:'| of | into |$)", subject)
    if m:
        return m.group(1)
    m = re.search(r"claude/[\w./-]+", subject)
    return m.group(0) if m else subject[:40]


def collect_merges() -> list[dict]:
    shas = git("rev-list", "--all", "--merges").split()
    merges = []
    for sha in shas:
        subject = git("log", "-1", "--format=%s", sha)
        parents = git("log", "-1", "--format=%P", sha).split()
        if len(parents) < 2:
            continue
        p1, p2 = parents[0], parents[1]
        try:
            base = git("merge-base", p1, p2)
        except subprocess.CalledProcessError:
            base = p1
        merges.append({
            "merge_date": dt(git("log", "-1", "--format=%aI", sha)),
            "start_date": dt(git("log", "-1", "--format=%aI", base)),
            "name": branch_name(subject),
        })
    merges.sort(key=lambda m: m["merge_date"])
    return merges


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    first = dt(git("log", "--all", "--reverse", "--format=%aI").splitlines()[0])
    last = dt(git("log", "--all", "--format=%aI").splitlines()[0])
    merges = collect_merges()

    fig, ax = plt.subplots(figsize=(14, 6.5))
    fig.suptitle(
        "SBTN Land LEAFs — branch & merge history",
        fontsize=17, fontweight="bold",
    )

    x0, x1 = mdates.date2num(first), mdates.date2num(last)

    # Main integration line (paper branch).
    ax.plot([x0, x1], [0, 0], color=MAIN_COLOR, lw=3, zorder=2)
    ax.scatter([x0, x1], [0, 0], s=90, color=MAIN_COLOR, zorder=3)
    ax.annotate("paper  (integration branch)", (x0, 0),
                textcoords="offset points", xytext=(4, 12),
                fontweight="bold", color=MAIN_COLOR)

    # Feature branches, staggered above/below to avoid overlap.
    for i, m in enumerate(merges):
        color = BRANCH_COLORS[i % len(BRANCH_COLORS)]
        y = 1.0 + (i % 3) * 0.9
        if i % 2:
            y = -y
        bx = mdates.date2num(m["start_date"])
        mx = mdates.date2num(m["merge_date"])

        # Branch line.
        ax.plot([bx, mx], [y, y], color=color, lw=2.4, zorder=2)
        ax.scatter([mx], [y], s=55, color=color, zorder=3)

        # Fork off the main line, and merge back into it (curved arrows).
        ax.add_patch(FancyArrowPatch(
            (bx, 0), (bx, y), connectionstyle="arc3,rad=0.25",
            arrowstyle="-", color=color, lw=1.6, alpha=0.8, zorder=1))
        ax.add_patch(FancyArrowPatch(
            (mx, y), (mx, 0), connectionstyle="arc3,rad=0.25",
            arrowstyle="-|>", mutation_scale=14, color=color, lw=1.8, zorder=1))

        span = (m["merge_date"] - m["start_date"]).days
        va = "bottom" if y > 0 else "top"
        ax.annotate(
            f"{m['name']}\n{m['start_date']:%d %b} → {m['merge_date']:%d %b}"
            f"  ({span}d)",
            ((bx + mx) / 2, y), textcoords="offset points",
            xytext=(0, 8 if y > 0 else -8), ha="center", va=va,
            fontsize=8.5, color=color, fontweight="bold")

    ax.set_ylim(-4.2, 4.2)
    ax.set_xlim(x0 - 3, x1 + 3)
    ax.get_yaxis().set_visible(False)
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.grid(axis="x", alpha=0.25)
    ax.set_title(
        f"{len(merges)} feature branches merged · "
        f"{first:%d %b %Y} – {last:%d %b %Y}",
        fontsize=11, color="#555",
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT, dpi=140)
    plt.close(fig)
    print(f"Wrote {OUT.relative_to(REPO_ROOT)}  ({len(merges)} branches)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
