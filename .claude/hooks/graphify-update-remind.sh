#!/usr/bin/env bash
# PostToolUse(Write|Edit): after a code file edit, remind Claude to keep the
# graphify graph current if one exists in this project.
REPORT="$(pwd)/graphify-out/GRAPH_REPORT.md"
[[ -f "$REPORT" ]] || exit 0

# Read the file path from the tool-input JSON on stdin.
FILE=$(python3 -c "
import json, sys
try:
    d = json.load(sys.stdin)
    print(d.get('tool_input', {}).get('file_path', ''))
except Exception:
    print('')
" 2>/dev/null)

[[ -z "$FILE" ]] && exit 0

# Only remind for code / script files — not markdown, JSON, YAML, data files.
case "$FILE" in
  *.py|*.ts|*.js|*.mjs|*.go|*.rs|*.java|*.cpp|*.c|*.cc|*.cxx|\
  *.h|*.hpp|*.rb|*.swift|*.kt|*.cs|*.scala|*.php|*.lua|\
  *.sh|*.bash|*.nf|*.smk|*.R|*.Rmd)
    echo "graphify: graph exists. If this edit adds, removes, or renames a public symbol or module, run /graphify --update when done."
    ;;
esac
exit 0
