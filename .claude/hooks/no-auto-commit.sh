#!/usr/bin/env bash
# PreToolUse hook for Bash.
# Blocks git commit / push / reset --hard / rebase / clean.
# Inspects staged + modified files, flags risky ones, and prints the
# exact git add + git commit commands the user should run.
# Reads tool input JSON on stdin; exit 2 = block, exit 0 = allow.

set -uo pipefail

cmd=$(python3 -c '
import json, sys
try:
    d = json.load(sys.stdin)
except Exception:
    sys.exit(0)
print((d.get("tool_input") or {}).get("command") or "")
')

[[ -z "$cmd" ]] && exit 0

if [[ "$cmd" =~ (^|[[:space:]\&\;\|])git[[:space:]]+(commit|push|reset[[:space:]]+--hard|rebase|clean) ]]; then

  python3 - "$cmd" >&2 <<'PYEOF'
import sys, subprocess, re, os
from pathlib import Path

cmd = sys.argv[1]

# --- Extract commit message Claude was about to use -------------------------
msg_match = re.search(r'-m\s+["\'](.+?)["\']', cmd) \
         or re.search(r'-m\s+(\S+)', cmd)
proposed_msg = msg_match.group(1) if msg_match else "<describe the change: what and why>"

# --- File patterns considered risky -----------------------------------------
RISKY_SUFFIXES = {
    '.log', '.tmp', '.pyc', '.o', '.a', '.so',
    '.h5', '.hdf5', '.pkl', '.pickle', '.npy', '.npz', '.pt', '.pth',
    '.gz', '.bz2', '.xz', '.zip', '.tar', '.7z',
    '.bam', '.sam', '.cram', '.fastq', '.fq', '.fasta', '.fa', '.vcf',
    '.bed', '.bigwig', '.bw', '.bigbed',
}
RISKY_DIRS = {
    'results', 'output', 'outputs', 'tmp', 'temp', '__pycache__',
    '.pytest_cache', 'graphify-out', 'data', 'raw',
}
RISKY_NAMES = {'.DS_Store', 'Thumbs.db'}

def is_risky(path_str):
    p = Path(path_str)
    if p.name in RISKY_NAMES:
        return True
    if p.suffix.lower() in RISKY_SUFFIXES:
        return True
    # any path component matches a risky directory name
    if any(part in RISKY_DIRS for part in p.parts):
        return True
    # large file threshold: 1 MB
    try:
        if p.exists() and p.stat().st_size > 1_000_000:
            return True
    except OSError:
        pass
    return False

# --- Inspect git status -----------------------------------------------------
try:
    result = subprocess.run(
        ['git', 'status', '--porcelain'],
        capture_output=True, text=True, check=True
    )
    lines = result.stdout.splitlines()
except Exception:
    lines = []

staged_safe, staged_risky, unstaged_safe, unstaged_risky = [], [], [], []

for line in lines:
    if len(line) < 4:
        continue
    xy = line[:2]
    # handle renames: "R  old -> new" — take the new name
    raw_path = line[3:].split(' -> ')[-1].strip().strip('"')
    index_status = xy[0]
    worktree_status = xy[1]

    if index_status not in (' ', '?') and index_status != '!':
        # staged
        if is_risky(raw_path):
            staged_risky.append(raw_path)
        else:
            staged_safe.append(raw_path)

    if worktree_status not in (' ', '?') and worktree_status != '!':
        # unstaged modification (not untracked)
        if is_risky(raw_path):
            unstaged_risky.append(raw_path)
        else:
            unstaged_safe.append(raw_path)

    if xy == '??':
        # untracked — flag risky ones, suggest adding safe ones
        if is_risky(raw_path):
            unstaged_risky.append(raw_path)
        else:
            unstaged_safe.append(raw_path)

# --- Print the block message -------------------------------------------------
print("BLOCKED: Claude must not run git commit / push / reset --hard / rebase / clean.")
print(f"Command attempted: {cmd}")
print()

if staged_risky:
    print("WARNING — these staged files look like temporary outputs or large data.")
    print("Unstage them before committing:")
    for f in staged_risky:
        print(f"  git restore --staged {f}")
    print()

if staged_risky and staged_safe:
    print("Safe files already staged:")
    for f in staged_safe:
        print(f"  {f}")
    print()

# Suggest what to add
files_to_add = staged_safe[:]  # already staged safe ones stay
if unstaged_safe:
    print("Suggested git add (safe files only):")
    add_targets = unstaged_safe
    print("  git add " + " ".join(add_targets))
    files_to_add.extend(unstaged_safe)
    print()

if unstaged_risky:
    print("Skipped (risky — do NOT add these):")
    for f in unstaged_risky:
        print(f"  {f}")
    print()

# Suggest commit command
print("Suggested commit:")
print(f'  git commit -m "{proposed_msg}"')
print()
print("You commit. Claude only proposes.")
PYEOF

  exit 2
fi

exit 0
