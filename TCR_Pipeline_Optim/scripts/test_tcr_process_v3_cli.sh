#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="$SCRIPT_DIR/tcr_process_v3.sh"

fail() {
    echo "FAIL: $1" >&2
    exit 1
}

assert_contains() {
    local haystack="$1"
    local needle="$2"
    local message="$3"
    if [[ "$haystack" != *"$needle"* ]]; then
        fail "$message"
    fi
}

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

info_csv="$tmpdir/info.csv"
rawfq_root="$tmpdir/rawfq"
work_root="$tmpdir/work"
download_script="$tmpdir/download_fq.py"
match_script="$tmpdir/match_sample.py"
mkdir -p "$rawfq_root" "$work_root"
touch "$download_script" "$match_script"

make_row() {
    local sed_id=$1
    local sample_id=$2
    local project=$3
    local adapter=$4
    local panel=$5
    local data_id=$6
    printf "%s,%s,LIB001,Name001,c5,c6,c7,%s,c9,c10,c11,%s,c13,c14,c15,%s,c17,c18,c19,c20,c21,c22,c23,%s\n" \
        "$sed_id" "$sample_id" "$project" "$adapter" "$panel" "$data_id"
}

make_row_with_name() {
    local sed_id=$1
    local sample_id=$2
    local name=$3
    local project=$4
    local adapter=$5
    local panel=$6
    local data_id=$7
    printf '%s,%s,LIB001,"%s",c5,c6,c7,%s,c9,c10,c11,%s,c13,c14,c15,%s,c17,c18,c19,c20,c21,c22,c23,%s\n' \
        "$sed_id" "$sample_id" "$name" "$project" "$adapter" "$panel" "$data_id"
}

{
    make_row "RUN001" "SAMPLE_A" "肺结节" "Adapter_TCR_V1" "PanelA" "cos://bucket/sample_a"
    make_row "RUN001" "SAMPLE_B" "其他项目" "Adapter_TCR_V1" "PanelB" "cos://bucket/sample_b"
    make_row "RUN002" "SAMPLE_C" "肺结节" "Adapter_TCR_V1" "PanelC" "cos://bucket/sample_c"
    make_row "RUN003" "SAMPLE_D" "肺结节" "Adapter_BCR_V1" "PanelD" "cos://bucket/sample_d"
    make_row_with_name "RUN004" "SAMPLE_E" "Name, With Comma" "肺结节" "Adapter_TCR_V1" "PanelE" "cos://bucket/sample_e"
    make_row "RUN005" "SAMPLE_LITERAL" "项目.精确" "Adapter.TCR" "PanelLiteral" "cos://bucket/sample_literal"
    make_row "RUN005" "SAMPLE_OTHER" "项目X精确" "AdapterXTCR" "PanelOther" "cos://bucket/sample_other"
} > "$info_csv"

help_output="$(bash "$SCRIPT" --help)"
assert_contains "$help_output" "analysis_operation_manual.md" "--help should point to the operation manual"
if [[ "$help_output" == *"USAGE:"* ]]; then
    fail "--help should not print a long usage block"
fi

set +e
old_mode_output="$(bash "$SCRIPT" --batch-id batch01 --sample-info "$info_csv" --dry-run 2>&1)"
old_mode_status=$?
set -e
if [[ "$old_mode_status" -eq 0 ]]; then
    fail "--sample-info mode should no longer be accepted"
fi
assert_contains "$old_mode_output" "unknown argument: --sample-info" "old sample-info mode should be rejected"
if [[ "$old_mode_output" == *"USAGE:"* ]]; then
    fail "argument errors should not print a long usage block"
fi

set +e
missing_output="$(bash "$SCRIPT" --info-csv "$info_csv" --chain TRB --prepare-only 2>&1)"
missing_status=$?
set -e
if [[ "$missing_status" -eq 0 ]]; then
    fail "missing --batch-id should exit non-zero"
fi
assert_contains "$missing_output" "ERROR: --batch-id is required" "missing --batch-id should explain the error"
if [[ "$missing_output" == *"USAGE:"* ]]; then
    fail "validation errors should not print a long usage block"
fi

bash "$SCRIPT" \
    --batch-id batch01 \
    --info-csv "$info_csv" \
    --sed-id RUN001 \
    --chain TRB \
    --rawfq-root "$rawfq_root" \
    --work-root "$work_root" \
    --download-script "$download_script" \
    --match-script "$match_script" \
    --python-bin python3 \
    --prepare-only >/dev/null

filtered_info="$work_root/batch01/filtered_info.csv"
if [ ! -f "$filtered_info" ]; then
    fail "prepare-only should write filtered_info.csv"
fi

line_count="$(wc -l < "$filtered_info" | tr -d ' ')"
if [ "$line_count" != "2" ]; then
    echo "filtered_info.csv:" >&2
    cat "$filtered_info" >&2
    fail "filtered_info.csv should contain header plus one matching sample"
fi

header="$(head -n 1 "$filtered_info")"
assert_contains "$header" "Chain" "filtered_info.csv should include Chain column"
assert_contains "$header" "Receptor" "filtered_info.csv should include Receptor column"

sample_line="$(tail -n 1 "$filtered_info")"
assert_contains "$sample_line" "SAMPLE_A" "filtered_info.csv should keep matching sample"
assert_contains "$sample_line" "TRB" "filtered_info.csv should append requested chain"
assert_contains "$sample_line" "TCR" "legacy runs should append default TCR receptor"

bash "$SCRIPT" \
    --batch-id batch-quoted-comma \
    --info-csv "$info_csv" \
    --sed-id RUN004 \
    --chain TRB \
    --rawfq-root "$rawfq_root" \
    --work-root "$work_root" \
    --download-script "$download_script" \
    --match-script "$match_script" \
    --python-bin python3 \
    --prepare-only >/dev/null

quoted_comma_line="$(tail -n 1 "$work_root/batch-quoted-comma/filtered_info.csv")"
assert_contains "$quoted_comma_line" "SAMPLE_E" "CSV filtering should keep rows with quoted commas"
assert_contains "$quoted_comma_line" '"Name, With Comma"' "CSV filtering should preserve quoted comma fields"

bash "$SCRIPT" \
    --batch-id batch-literal-match \
    --info-csv "$info_csv" \
    --sed-id RUN005 \
    --chain TRB \
    --project-keyword . \
    --adapter-pattern . \
    --rawfq-root "$rawfq_root" \
    --work-root "$work_root" \
    --download-script "$download_script" \
    --match-script "$match_script" \
    --python-bin python3 \
    --prepare-only >/dev/null

literal_filtered_info="$work_root/batch-literal-match/filtered_info.csv"
literal_line_count="$(wc -l < "$literal_filtered_info" | tr -d ' ')"
if [ "$literal_line_count" != "2" ]; then
    echo "filtered_info.csv:" >&2
    cat "$literal_filtered_info" >&2
    fail "adapter and project matching should treat CLI values as literal substrings"
fi
literal_sample_line="$(tail -n 1 "$literal_filtered_info")"
assert_contains "$literal_sample_line" "SAMPLE_LITERAL" "literal matching should keep rows containing a literal dot"
if [[ "$literal_sample_line" == *"SAMPLE_OTHER"* ]]; then
    fail "literal matching should exclude rows without a literal dot"
fi

set +e
missing_bcr_adapter_output="$(
    bash "$SCRIPT" \
        --batch-id batch-bcr-missing-adapter \
        --info-csv "$info_csv" \
        --sed-id RUN003 \
        --receptor BCR \
        --chain BCR_ALL \
        --rawfq-root "$rawfq_root" \
        --work-root "$work_root" \
        --download-script "$download_script" \
        --match-script "$match_script" \
        --python-bin python3 \
        --prepare-only 2>&1
)"
missing_bcr_adapter_status=$?
set -e
if [ "$missing_bcr_adapter_status" -eq 0 ]; then
    fail "BCR runs without --adapter-pattern should fail"
fi
assert_contains "$missing_bcr_adapter_output" "ERROR: --adapter-pattern is required for BCR runs" "BCR adapter validation should explain the error"

bash "$SCRIPT" \
    --batch-id batch-bcr \
    --info-csv "$info_csv" \
    --sed-id RUN003 \
    --receptor bcr \
    --chain bcr_all \
    --adapter-pattern Adapter_BCR_ \
    --rawfq-root "$rawfq_root" \
    --work-root "$work_root" \
    --download-script "$download_script" \
    --match-script "$match_script" \
    --python-bin python3 \
    --prepare-only >/dev/null

bcr_filtered_info="$work_root/batch-bcr/filtered_info.csv"
bcr_sample_line="$(tail -n 1 "$bcr_filtered_info")"
assert_contains "$bcr_sample_line" "SAMPLE_D" "BCR prepare-only should keep matching sample"
assert_contains "$bcr_sample_line" "BCR_ALL" "BCR prepare-only should append requested chain selection"
assert_contains "$bcr_sample_line" "BCR" "BCR prepare-only should append normalized receptor"

bcr_igh_dry_run_output="$(
    bash "$SCRIPT" \
        --batch-id batch-bcr-igh \
        --info-csv "$info_csv" \
        --sed-id RUN003 \
        --receptor BCR \
        --chain IGH \
        --adapter-pattern Adapter_BCR_ \
        --rawfq-root "$rawfq_root" \
        --work-root "$work_root" \
        --download-script "$download_script" \
        --match-script "$match_script" \
        --python-bin python3 \
        --dry-run
)"
assert_contains "$bcr_igh_dry_run_output" "Receptor: BCR" "BCR IGH dry-run should show BCR receptor"
assert_contains "$bcr_igh_dry_run_output" "Chain: IGH" "BCR IGH dry-run should accept a single BCR chain"

set +e
invalid_receptor_output="$(
    bash "$SCRIPT" \
        --batch-id batch-invalid-receptor \
        --info-csv "$info_csv" \
        --sed-id RUN003 \
        --receptor XYZ \
        --chain IGH \
        --adapter-pattern Adapter_BCR_ \
        --rawfq-root "$rawfq_root" \
        --work-root "$work_root" \
        --download-script "$download_script" \
        --match-script "$match_script" \
        --python-bin python3 \
        --prepare-only 2>&1
)"
invalid_receptor_status=$?
set -e
if [ "$invalid_receptor_status" -eq 0 ]; then
    fail "invalid receptor values should fail"
fi
assert_contains "$invalid_receptor_output" "ERROR: --receptor must be TCR or BCR, got 'XYZ'" "invalid receptor should explain the error"

set +e
invalid_combo_output="$(
    bash "$SCRIPT" \
        --batch-id batch-invalid \
        --info-csv "$info_csv" \
        --sed-id RUN003 \
        --receptor TCR \
        --chain IGH \
        --rawfq-root "$rawfq_root" \
        --work-root "$work_root" \
        --download-script "$download_script" \
        --match-script "$match_script" \
        --python-bin python3 \
        --prepare-only 2>&1
)"
invalid_combo_status=$?
set -e
if [ "$invalid_combo_status" -eq 0 ]; then
    fail "invalid receptor/chain combinations should fail"
fi
assert_contains "$invalid_combo_output" "ERROR: chain 'IGH' is not valid for receptor 'TCR'" "invalid combination should explain the error"

dry_run_output="$(
    bash "$SCRIPT" \
        -b batch01 \
        --info-csv "$info_csv" \
        --sed-id RUN001 \
        --chain BOTH \
        --rawfq-root "$rawfq_root" \
        --work-root "$work_root" \
        --download-script "$download_script" \
        --match-script "$match_script" \
        --python-bin python3 \
        --skip-download \
        --dry-run
)"

assert_contains "$dry_run_output" "Dry run: info-mode configuration validated." "dry-run should report success"
assert_contains "$dry_run_output" "Filtered info: $work_root/batch01/filtered_info.csv" "dry-run should show filtered_info path"
assert_contains "$dry_run_output" "Chain: BOTH" "dry-run should show chain"
assert_contains "$dry_run_output" "Receptor: TCR" "legacy dry-run should show default TCR receptor"

echo "PASS: tcr_process_v3 CLI tests"
