#!/usr/bin/env bash
# Before file searches: remind Claude the graph exists so it reads it once
# and navigates by structure instead of crawling raw files blindly.
# Prints a one-line reminder only — not the full report (that would add
# hundreds of tokens on every search and negate the savings).
REPORT="$(pwd)/graphify-out/GRAPH_REPORT.md"
[[ -f "$REPORT" ]] && echo "graphify: graph available — read graphify-out/GRAPH_REPORT.md for god nodes and community structure before searching files."
exit 0
