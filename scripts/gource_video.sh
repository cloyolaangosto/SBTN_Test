#!/usr/bin/env bash
#
# Render a Gource "growing file-tree" animation of this repository's history
# to an MP4. Gource shows every file as a node in a tree, with contributors
# moving around and building/editing files as their commits are replayed --
# the most vivid way to "show all the work done" over the life of the project.
#
# This needs OpenGL and is best run on your own machine (macOS/Linux desktop),
# NOT in a headless CI container.
#
# Install the tools first:
#   macOS:         brew install gource ffmpeg
#   Debian/Ubuntu: sudo apt-get install gource ffmpeg
#
# Then, from the repo root:
#   ./scripts/gource_video.sh                # -> documentation/figures/gource.mp4
#   ./scripts/gource_video.sh my_video.mp4   # custom output path
#
set -euo pipefail

OUT="${1:-documentation/figures/gource.mp4}"
mkdir -p "$(dirname "$OUT")"

# --seconds-per-day  : playback speed (lower = faster)
# --auto-skip-seconds: skip idle gaps so the video stays lively
# --key              : show the file-extension colour key
gource \
  --seconds-per-day 0.75 \
  --auto-skip-seconds 0.5 \
  --max-file-lag 0.2 \
  --hide filenames,mouse,progress \
  --key \
  --title "SBTN Land LEAFs — project history" \
  --highlight-users \
  --font-size 20 \
  -1280x720 \
  --output-ppm-stream - \
  | ffmpeg -y -r 30 -f image2pipe -vcodec ppm -i - \
      -vcodec libx264 -preset medium -pix_fmt yuv420p -crf 21 \
      "$OUT"

echo "Wrote $OUT"
