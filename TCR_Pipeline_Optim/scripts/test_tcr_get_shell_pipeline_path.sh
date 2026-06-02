#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_SCRIPT="$SCRIPT_DIR/tcr_get_shell_fixed_primers_v3.py"

fail() {
    echo "FAIL: $1" >&2
    exit 1
}

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

copied_scripts="$tmpdir/copied_scripts"
work_dir="$tmpdir/work"
raw_dir="$tmpdir/raw"
mkdir -p "$copied_scripts" "$work_dir" "$raw_dir"

cp "$SOURCE_SCRIPT" "$copied_scripts/tcr_get_shell_fixed_primers_v3.py"
touch "$copied_scripts/TCR_analysis_pipeline.v7.sh"

sample_info="$tmpdir/fq_matched.tsv"
r1="$raw_dir/sample1_R1.fq.gz"
r2="$raw_dir/sample1_R2.fq.gz"
touch "$r1" "$r2"

printf "Sample_ID\tChain\tRaw_Path_R1\tRaw_Path_R2\nsample1\tTRB\t%s\t%s\n" "$r1" "$r2" > "$sample_info"

python3 "$copied_scripts/tcr_get_shell_fixed_primers_v3.py" "$sample_info" "$work_dir" 8 generic-amplicon >/dev/null

expected_pipeline="$copied_scripts/TCR_analysis_pipeline.v7.sh"
first_command="$(head -n 1 "$work_dir/get_shell.sh")"

if [[ "$first_command" != "$expected_pipeline "* ]]; then
    echo "Expected command to start with: $expected_pipeline" >&2
    echo "Actual command: $first_command" >&2
    fail "get_shell.sh should use TCR_analysis_pipeline.v7.sh from the same directory as tcr_get_shell_fixed_primers_v3.py"
fi

echo "PASS: tcr_get_shell pipeline path test"
