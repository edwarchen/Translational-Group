#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_SCRIPT="$SCRIPT_DIR/tcr_process_v3.sh"
SOURCE_FILTER="$SCRIPT_DIR/filter_info.py"

fail() {
    echo "FAIL: $1" >&2
    exit 1
}

assert_contains() {
    local haystack="$1"
    local needle="$2"
    local message="$3"
    if [[ "$haystack" != *"$needle"* ]]; then
        echo "$haystack" >&2
        fail "$message"
    fi
}

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

copied_scripts="$tmpdir/scripts"
rawfq_root="$tmpdir/rawfq"
work_root="$tmpdir/work"
info_csv="$tmpdir/info.csv"
mkdir -p "$copied_scripts" "$rawfq_root" "$work_root"
cp "$SOURCE_SCRIPT" "$copied_scripts/tcr_process_v3.sh"
cp "$SOURCE_FILTER" "$copied_scripts/filter_info.py"

generator_empty="$tmpdir/generator_empty.py"
generator_noop="$tmpdir/generator_noop.py"

cat > "$generator_empty" <<'PY'
import os
import sys

open(os.path.join(sys.argv[2], "get_shell.sh"), "w").close()
PY

cat > "$generator_noop" <<'PY'
PY

cp "$generator_empty" "$copied_scripts/tcr_get_shell_fixed_primers_v3.py"

cat > "$tmpdir/matcher_header_only.py" <<'PY'
import sys

with open(sys.argv[3], "w") as output_fh:
    output_fh.write("Sample_ID\tChain\tRaw_Path_R1\tRaw_Path_R2\n")
PY

cat > "$tmpdir/matcher_with_row.py" <<'PY'
import sys

with open(sys.argv[3], "w") as output_fh:
    output_fh.write("Sample_ID\tChain\tRaw_Path_R1\tRaw_Path_R2\n")
    output_fh.write("SAMPLE_A\tTRB\t/tmp/sample_R1.fq.gz\t/tmp/sample_R2.fq.gz\n")
PY

cat > "$tmpdir/matcher_noop.py" <<'PY'
PY

printf "RUN001,SAMPLE_A,LIB001,Name001,c5,c6,c7,肺结节,c9,c10,c11,Adapter_TCR_V1,c13,c14,c15,PanelA,c17,c18,c19,c20,c21,c22,c23,cos://bucket/sample_a\n" > "$info_csv"

run_batch() {
    local batch_id=$1
    local match_script=$2
    bash "$copied_scripts/tcr_process_v3.sh" \
        --batch-id "$batch_id" \
        --info-csv "$info_csv" \
        --sed-id RUN001 \
        --chain TRB \
        --rawfq-root "$rawfq_root" \
        --work-root "$work_root" \
        --match-script "$match_script" \
        --python-bin python3 \
        --skip-download
}

set +e
header_only_output="$(run_batch batch-header-only "$tmpdir/matcher_header_only.py" 2>&1)"
header_only_status=$?
set -e
if [ "$header_only_status" -eq 0 ]; then
    fail "header-only fq_matched.tsv should fail"
fi
assert_contains "$header_only_output" "ERROR: match_sample.py produced no matched FASTQ rows" "header-only fq_matched.tsv should fail before task generation"

set +e
empty_shell_output="$(run_batch batch-empty-shell "$tmpdir/matcher_with_row.py" 2>&1)"
empty_shell_status=$?
set -e
if [ "$empty_shell_status" -eq 0 ]; then
    fail "empty get_shell.sh should fail"
fi
assert_contains "$empty_shell_output" "ERROR: tcr_get_shell_fixed_primers_v3.py produced no analysis commands" "empty get_shell.sh should fail before add_wait"

stale_match_batch="$work_root/batch-stale-match"
mkdir -p "$stale_match_batch"
printf "Sample_ID\tChain\tRaw_Path_R1\tRaw_Path_R2\nSTALE\tTRB\t/tmp/stale_R1.fq.gz\t/tmp/stale_R2.fq.gz\n" > "$stale_match_batch/fq_matched.tsv"

set +e
stale_match_output="$(run_batch batch-stale-match "$tmpdir/matcher_noop.py" 2>&1)"
stale_match_status=$?
set -e
if [ "$stale_match_status" -eq 0 ]; then
    fail "stale fq_matched.tsv should not satisfy a no-op matcher run"
fi
assert_contains "$stale_match_output" "ERROR: match_sample.py produced no matched FASTQ rows" "stale fq_matched.tsv should be removed before matching"

cp "$generator_noop" "$copied_scripts/tcr_get_shell_fixed_primers_v3.py"
stale_shell_batch="$work_root/batch-stale-shell"
mkdir -p "$stale_shell_batch"
printf "echo stale command\n" > "$stale_shell_batch/get_shell.sh"

set +e
stale_shell_output="$(run_batch batch-stale-shell "$tmpdir/matcher_with_row.py" 2>&1)"
stale_shell_status=$?
set -e
if [ "$stale_shell_status" -eq 0 ]; then
    fail "stale get_shell.sh should not satisfy a no-op generator run"
fi
assert_contains "$stale_shell_output" "ERROR: tcr_get_shell_fixed_primers_v3.py produced no analysis commands" "stale get_shell.sh should be removed before task generation"

echo "PASS: tcr_process_v3 generation guard tests"
