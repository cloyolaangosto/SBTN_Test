#!/usr/bin/env python3
"""Visualize the evolution of this repository's work over time.

Reads the git history (all branches) and produces:

  1. documentation/figures/activity_dashboard.png
       A poster-style static summary: cumulative commits, commits per month,
       commits per contributor, and lines added/removed over time.

  2. documentation/figures/history.gif  (or .mp4 if ffmpeg is available)
       An animated commit-by-commit timeline showing the project growing:
       cumulative commits and total files tracked, coloured by contributor.

No third-party git bindings are required -- the git CLI is called via
subprocess. Only matplotlib and Pillow are needed:

    pip install matplotlib pillow

Usage:
    python scripts/visualize_history.py [--no-gif] [--fps 20]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = REPO_ROOT / "documentation" / "figures"

# A colour-blind-friendly qualitative palette.
PALETTE = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756",
    "#72B7B2", "#B279A2", "#EECA3B", "#9D755D",
]


def _git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL
    )


class Commit:
    __slots__ = ("sha", "author", "when", "adds", "dels", "files")

    def __init__(self, sha, author, when, adds, dels, files):
        self.sha = sha
        self.author = author
        self.when = when
        self.adds = adds
        self.dels = dels
        self.files = files


def load_commits() -> list[Commit]:
    """Parse `git log` (all branches, chronological) with per-commit numstat."""
    sep = "\x1f"  # unit separator, unlikely to appear in author names
    fmt = f"__C__%H{sep}%an{sep}%aI"
    raw = _git(
        "log", "--all", "--no-merges", "--reverse",
        "--numstat", f"--pretty=format:{fmt}",
    )
    commits: list[Commit] = []
    cur = None
    adds = dels = 0
    files: set[str] = set()

    def flush():
        nonlocal cur, adds, dels, files
        if cur is not None:
            sha, author, when = cur
            commits.append(
                Commit(sha, author, when, adds, dels, set(files))
            )

    for line in raw.splitlines():
        if line.startswith("__C__"):
            flush()
            body = line[len("__C__"):]
            sha, author, when = body.split(sep)
            cur = (sha, author, datetime.fromisoformat(when))
            adds = dels = 0
            files = set()
        elif line.strip():
            parts = line.split("\t")
            if len(parts) == 3:
                a, d, path = parts
                adds += int(a) if a.isdigit() else 0
                dels += int(d) if d.isdigit() else 0
                files.add(path)
    flush()
    commits.sort(key=lambda c: c.when)
    return commits


def build_dashboard(commits: list[Commit]) -> Path:
    authors = list(dict.fromkeys(c.author for c in commits))
    color_for = {a: PALETTE[i % len(PALETTE)] for i, a in enumerate(authors)}

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(
        "SBTN Land LEAFs — Project Activity",
        fontsize=18, fontweight="bold",
    )

    # 1. Cumulative commits over time.
    ax = axes[0][0]
    dates = [c.when for c in commits]
    ax.plot(dates, range(1, len(commits) + 1), color=PALETTE[0], lw=2.2)
    ax.fill_between(dates, range(1, len(commits) + 1), color=PALETTE[0], alpha=0.15)
    ax.set_title("Cumulative commits", fontweight="bold")
    ax.set_ylabel("commits")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.grid(alpha=0.3)

    # 2. Commits per month.
    ax = axes[0][1]
    per_month: Counter = Counter(c.when.strftime("%Y-%m") for c in commits)
    months = sorted(per_month)
    ax.bar(months, [per_month[m] for m in months], color=PALETTE[2])
    ax.set_title("Commits per month", fontweight="bold")
    ax.set_ylabel("commits")
    ax.tick_params(axis="x", rotation=45)
    ax.grid(alpha=0.3, axis="y")

    # 3. Commits per contributor.
    ax = axes[1][0]
    per_author: Counter = Counter(c.author for c in commits)
    names = [a for a, _ in per_author.most_common()]
    ax.barh(
        names, [per_author[a] for a in names],
        color=[color_for[a] for a in names],
    )
    ax.set_title("Commits per contributor", fontweight="bold")
    ax.set_xlabel("commits")
    ax.invert_yaxis()
    ax.grid(alpha=0.3, axis="x")

    # 4. Lines added / removed over time (cumulative).
    ax = axes[1][1]
    cum_add, cum_del = [], []
    a = d = 0
    for c in commits:
        a += c.adds
        d += c.dels
        cum_add.append(a)
        cum_del.append(d)
    ax.plot(dates, cum_add, color=PALETTE[2], lw=2, label="added")
    ax.plot(dates, cum_del, color=PALETTE[3], lw=2, label="removed")
    ax.set_title("Cumulative lines changed", fontweight="bold")
    ax.set_ylabel("lines")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = OUT_DIR / "activity_dashboard.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def build_animation(commits: list[Commit], fps: int) -> Path:
    """Animated timeline of cumulative commits + files tracked."""
    authors = list(dict.fromkeys(c.author for c in commits))
    color_for = {a: PALETTE[i % len(PALETTE)] for i, a in enumerate(authors)}

    dates = [c.when for c in commits]
    n = len(commits)

    # Running count of distinct files touched so far (proxy for project size).
    seen: set[str] = set()
    files_curve = []
    for c in commits:
        seen |= c.files
        files_curve.append(len(seen))

    fig, ax1 = plt.subplots(figsize=(12, 6.75))
    ax2 = ax1.twinx()
    fig.suptitle(
        "SBTN Land LEAFs — repository growth", fontsize=16, fontweight="bold"
    )

    ax1.set_xlim(dates[0], dates[-1])
    ax1.set_ylim(0, n * 1.05)
    ax2.set_ylim(0, max(files_curve) * 1.1)
    ax1.set_ylabel("cumulative commits", color=PALETTE[0])
    ax2.set_ylabel("distinct files touched", color=PALETTE[1])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax1.grid(alpha=0.3)

    (commits_line,) = ax1.plot([], [], color=PALETTE[0], lw=2.4)
    (files_line,) = ax2.plot([], [], color=PALETTE[1], lw=2.0, ls="--")
    scat = ax1.scatter([], [], s=[], c=[], zorder=5, edgecolors="white", lw=0.4)
    caption = ax1.text(
        0.02, 0.95, "", transform=ax1.transAxes, va="top",
        fontsize=11, fontweight="bold",
    )

    # Legend for contributors.
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=color_for[a], label=a)
        for a in authors
    ]
    ax1.legend(handles=handles, loc="lower right", fontsize=9)

    # Hold on the final frame for ~1.5s.
    hold = int(fps * 1.5)
    total_frames = n + hold

    def update(frame):
        i = min(frame, n)
        xs = dates[:i]
        commits_line.set_data(xs, range(1, i + 1))
        files_line.set_data(xs, files_curve[:i])
        if i:
            offsets = [[mdates.date2num(d), y] for d, y in zip(xs, range(1, i + 1))]
            scat.set_offsets(offsets)
            scat.set_color([color_for[c.author] for c in commits[:i]])
            scat.set_sizes([18] * i)
            caption.set_text(
                f"{commits[i-1].when:%d %b %Y}   "
                f"{i}/{n} commits   {files_curve[i-1]} files"
            )
        return commits_line, files_line, scat, caption

    anim = FuncAnimation(
        fig, update, frames=total_frames, interval=1000 / fps, blit=False
    )

    out = OUT_DIR / "history.gif"
    anim.save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-gif", action="store_true", help="skip the animation")
    ap.add_argument("--fps", type=int, default=20, help="animation frames/sec")
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    commits = load_commits()
    if not commits:
        print("No commits found.", file=sys.stderr)
        return 1
    print(f"Loaded {len(commits)} commits "
          f"({commits[0].when:%Y-%m-%d} → {commits[-1].when:%Y-%m-%d}).")

    dash = build_dashboard(commits)
    print(f"Wrote {dash.relative_to(REPO_ROOT)}")

    if not args.no_gif:
        gif = build_animation(commits, args.fps)
        print(f"Wrote {gif.relative_to(REPO_ROOT)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
