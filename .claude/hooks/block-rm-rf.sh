#!/usr/bin/env bash
# PreToolUse hook for Bash.
# Blocks broad destructive rm / find -delete / dd / mkfs / shred / wipefs.
# Reads tool input JSON on stdin; exit 2 = block, exit 0 = allow.

set -euo pipefail

cmd=$(python3 -c '
import json, sys
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
print((d.get("tool_input") or {}).get("command") or "")
')

if [[ -z "$cmd" ]]; then
  exit 0
fi

block() {
  echo "BLOCKED: destructive command refused by block-rm-rf.sh." >&2
  echo "Command: $cmd" >&2
  echo "$1" >&2
  echo "If you truly mean to do this, run it yourself in a terminal." >&2
  exit 2
}

if [[ "$cmd" =~ (^|[[:space:]\&\;\|])rm[[:space:]]+(-[a-zA-Z]*r[a-zA-Z]*f|-[a-zA-Z]*f[a-zA-Z]*r|-rf|-fr) ]]; then
  block "rm -rf / -fr is not allowed from Claude."
fi

if [[ "$cmd" =~ find[[:space:]].*(-delete|-exec[[:space:]]+rm) ]]; then
  block "find -delete / find -exec rm is not allowed from Claude."
fi

if [[ "$cmd" =~ (^|[[:space:]\&\;\|])(dd|mkfs[.a-zA-Z0-9]*|shred|wipefs)[[:space:]] ]]; then
  block "Disk-level commands (dd/mkfs/shred/wipefs) are not allowed from Claude."
fi

exit 0
