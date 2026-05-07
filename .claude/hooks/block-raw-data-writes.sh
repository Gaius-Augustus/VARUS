#!/usr/bin/env bash
# PreToolUse hook for Write|Edit.
# Blocks any write to paths under raw/, data/raw/, or *_raw/.
# Reads tool input JSON on stdin; exit 2 = block, exit 0 = allow.

set -euo pipefail

path=$(python3 -c '
import json, sys
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
ti = d.get("tool_input", {}) or {}
print(ti.get("file_path") or ti.get("path") or "")
')

if [[ -z "$path" ]]; then
  exit 0
fi

case "$path" in
  */raw/*|raw/*|*/data/raw/*|data/raw/*|*_raw/*|*/_raw/*)
    echo "BLOCKED: writes under raw/ are forbidden (path: $path). Raw data is read-only. If you truly need to modify this file, do it manually outside Claude." >&2
    exit 2
    ;;
esac

exit 0
