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

copied_scripts="$tmpdir/copied scripts"
work_dir="$tmpdir/work dir"
raw_dir="$tmpdir/raw dir"
mkdir -p "$copied_scripts" "$work_dir" "$raw_dir"

cp "$SOURCE_SCRIPT" "$copied_scripts/tcr_get_shell_fixed_primers_v3.py"
cat > "$copied_scripts/TCR_analysis_pipeline.v7.sh" <<'SH'
#!/bin/bash
printf '%s\n' "$@" > "$(dirname "$0")/received_args.txt"
SH
chmod 0644 "$copied_scripts/TCR_analysis_pipeline.v7.sh"

sample_info="$tmpdir/fq_matched.tsv"
r1="$raw_dir/sample 1_R1.fq.gz"
r2="$raw_dir/sample 1_R2.fq.gz"
touch "$r1" "$r2"

sample_id="sample 1"
printf "Sample_ID\tChain\tRaw_Path_R1\tRaw_Path_R2\n%s\tTRB\t%s\t%s\n" "$sample_id" "$r1" "$r2" > "$sample_info"

python3 "$copied_scripts/tcr_get_shell_fixed_primers_v3.py" "$sample_info" "$work_dir" 8 generic-amplicon >/dev/null

expected_pipeline="$copied_scripts/TCR_analysis_pipeline.v7.sh"
first_command="$(head -n 1 "$work_dir/get_shell.sh")"
expected_prefix="bash '$expected_pipeline' "

if [[ "$first_command" != "$expected_prefix"* ]]; then
    echo "Expected command to start with: $expected_prefix" >&2
    echo "Actual command: $first_command" >&2
    fail "get_shell.sh should invoke the quoted sibling pipeline through bash"
fi

if [ -x "$expected_pipeline" ]; then
    fail "pipeline fixture should remain non-executable"
fi

bash "$work_dir/get_shell.sh"

received_args="$copied_scripts/received_args.txt"
if [ ! -f "$received_args" ]; then
    fail "generated command should execute the non-executable pipeline through bash"
fi

sample_dir="$work_dir/$sample_id"
expected_args="$(printf "%s\n" "$r1" "$r2" "$sample_dir" "$sample_id" 8 TRB generic-amplicon TCR)"
actual_args="$(cat "$received_args")"
if [ "$actual_args" != "$expected_args" ]; then
    echo "Expected arguments:" >&2
    printf "%s\n" "$expected_args" >&2
    echo "Actual arguments:" >&2
    printf "%s\n" "$actual_args" >&2
    fail "generated command should preserve shell argument boundaries"
fi

echo "PASS: tcr_get_shell pipeline path test"
