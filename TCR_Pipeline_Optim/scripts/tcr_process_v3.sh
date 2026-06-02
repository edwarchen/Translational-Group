#!/bin/bash
set -euo pipefail

# TCR batch pipeline entry point.
# Input starts from the latest sequencing info.csv, not from a pre-matched sample table.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

batch_id=""
info_csv=""
sed_id=""
receptor="TCR"
chain=""
project_keyword="肺结节"
adapter_pattern="Adapter_TCR_"
adapter_pattern_explicit=0
rawfq_dir_total=/haplox/rawfq/TCR
working_dir=/haplox/users/chenya/TCR_TEST
threads=5
PRESET=generic-amplicon
PYTHON_BIN=/x03_haplox/users/donglf/miniconda3/bin/python
download_script=/x03_haplox/users/donglf/common_tools/download_scripts/download_fq.py
match_sample_script=/x03_haplox/users/donglf/haima_pipeline_tools/match_sample.py
skip_download=0
prepare_only=0
dry_run=0

TCR_GET_SHELL_SCRIPT="$SCRIPT_DIR/tcr_get_shell_fixed_primers_v3.py"
FILTER_INFO_SCRIPT="$SCRIPT_DIR/filter_info.py"
ADD_WAIT_SCRIPT=/x03_haplox/users/donglf/common_tools/add_wait.py
TCR_QC_SCRIPT=/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/tcr_qc_allchain.py
FASTP_QC_SCRIPT=/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/fastp_qc.v2.py

function help_message(){
    echo "See analysis_operation_manual.md for usage and troubleshooting."
}

function require_value(){
    local option=$1
    local value=${2:-}
    if [ -z "$value" ] || [[ "$value" == -* ]]; then
        echo "ERROR: $option requires a value" >&2
        exit 1
    fi
}

while [ $# -gt 0 ]; do
    case "$1" in
        -b|--batch-id)
            require_value "$1" "${2:-}"
            batch_id=$2
            shift 2
            ;;
        --info-csv)
            require_value "$1" "${2:-}"
            info_csv=$2
            shift 2
            ;;
        --sed-id)
            require_value "$1" "${2:-}"
            sed_id=$2
            shift 2
            ;;
        --chain)
            require_value "$1" "${2:-}"
            chain=$2
            shift 2
            ;;
        --receptor)
            require_value "$1" "${2:-}"
            receptor=$2
            shift 2
            ;;
        --project-keyword)
            require_value "$1" "${2:-}"
            project_keyword=$2
            shift 2
            ;;
        --adapter-pattern)
            require_value "$1" "${2:-}"
            adapter_pattern=$2
            adapter_pattern_explicit=1
            shift 2
            ;;
        -t|--threads)
            require_value "$1" "${2:-}"
            threads=$2
            shift 2
            ;;
        -p|--preset)
            require_value "$1" "${2:-}"
            PRESET=$2
            shift 2
            ;;
        --rawfq-root)
            require_value "$1" "${2:-}"
            rawfq_dir_total=$2
            shift 2
            ;;
        --work-root)
            require_value "$1" "${2:-}"
            working_dir=$2
            shift 2
            ;;
        --download-script)
            require_value "$1" "${2:-}"
            download_script=$2
            shift 2
            ;;
        --match-script)
            require_value "$1" "${2:-}"
            match_sample_script=$2
            shift 2
            ;;
        --python-bin)
            require_value "$1" "${2:-}"
            PYTHON_BIN=$2
            shift 2
            ;;
        --skip-download)
            skip_download=1
            shift
            ;;
        --prepare-only)
            prepare_only=1
            shift
            ;;
        --dry-run)
            dry_run=1
            shift
            ;;
        -h|--help)
            help_message
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [ -z "$batch_id" ]; then
    echo "ERROR: --batch-id is required" >&2
    exit 1
fi

if [ -z "$info_csv" ]; then
    echo "ERROR: --info-csv is required" >&2
    exit 1
fi

if [ -z "$sed_id" ]; then
    echo "ERROR: --sed-id is required" >&2
    exit 1
fi

receptor="$(echo "$receptor" | tr '[:lower:]' '[:upper:]')"
chain="$(echo "$chain" | tr '[:lower:]' '[:upper:]')"

if [ "$receptor" != "TCR" ] && [ "$receptor" != "BCR" ]; then
    echo "ERROR: --receptor must be TCR or BCR, got '$receptor'" >&2
    exit 1
fi

case "$receptor:$chain" in
    TCR:TRA|TCR:TRB|TCR:BOTH|BCR:IGH|BCR:IGK|BCR:IGL|BCR:BCR_ALL)
        ;;
    *)
        echo "ERROR: chain '$chain' is not valid for receptor '$receptor'" >&2
        exit 1
        ;;
esac

if [ "$receptor" = "BCR" ] && [ "$adapter_pattern_explicit" -eq 0 ]; then
    echo "ERROR: --adapter-pattern is required for BCR runs" >&2
    exit 1
fi

if [ "$PRESET" != "generic-amplicon" ] && [ "$PRESET" != "rna-seq" ]; then
    echo "ERROR: --preset must be generic-amplicon or rna-seq, got '$PRESET'" >&2
    exit 1
fi

if ! [[ "$threads" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: --threads must be a positive integer, got '$threads'" >&2
    exit 1
fi

if [ ! -f "$info_csv" ]; then
    echo "ERROR: info.csv not found: $info_csv" >&2
    exit 1
fi

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo "ERROR: Python interpreter not found: $PYTHON_BIN" >&2
    exit 1
fi

if [ "$skip_download" -eq 0 ] && [ "$prepare_only" -eq 0 ] && [ ! -f "$download_script" ]; then
    echo "ERROR: download_fq.py not found: $download_script" >&2
    exit 1
fi

if [ ! -f "$match_sample_script" ]; then
    echo "ERROR: match_sample.py not found: $match_sample_script" >&2
    exit 1
fi

if [ ! -f "$TCR_GET_SHELL_SCRIPT" ]; then
    echo "ERROR: tcr_get_shell_fixed_primers_v3.py not found: $TCR_GET_SHELL_SCRIPT" >&2
    exit 1
fi

if [ ! -f "$FILTER_INFO_SCRIPT" ]; then
    echo "ERROR: filter_info.py not found: $FILTER_INFO_SCRIPT" >&2
    exit 1
fi

rawfq_dir_batch=$rawfq_dir_total/$batch_id
working_dir_batch=$working_dir/$batch_id
filtered_info_csv=$working_dir_batch/filtered_info.csv

if [ "$dry_run" -eq 1 ]; then
    echo "Dry run: info-mode configuration validated."
    echo "Batch ID: $batch_id"
    echo "Info CSV: $info_csv"
    echo "Sed ID: $sed_id"
    echo "Project keyword: $project_keyword"
    echo "Adapter pattern: $adapter_pattern"
    echo "Receptor: $receptor"
    echo "Chain: $chain"
    echo "Preset: $PRESET"
    echo "Threads: $threads"
    echo "FASTQ batch dir: $rawfq_dir_batch"
    echo "Working batch dir: $working_dir_batch"
    echo "Filtered info: $filtered_info_csv"
    echo "Download: $([ "$skip_download" -eq 1 ] && echo "skip" || echo "enabled")"
    exit 0
fi

mkdir -p "$working_dir_batch" "$working_dir_batch/log" "$rawfq_dir_batch"

"$PYTHON_BIN" "$FILTER_INFO_SCRIPT" \
    --input "$info_csv" \
    --output "$filtered_info_csv" \
    --sed-id "$sed_id" \
    --adapter-pattern "$adapter_pattern" \
    --project-keyword "$project_keyword" \
    --chain "$chain" \
    --receptor "$receptor"

filtered_count=$(( $(wc -l < "$filtered_info_csv" | tr -d ' ') - 1 ))
if [ "$filtered_count" -le 0 ]; then
    echo "ERROR: no samples matched sed_id=$sed_id adapter=$adapter_pattern project=$project_keyword" >&2
    echo "Filtered info: $filtered_info_csv" >&2
    exit 1
fi

echo "Filtered samples: $filtered_count"
echo "Filtered info: $filtered_info_csv"

if [ "$prepare_only" -eq 1 ]; then
    echo "Prepare-only mode: stop before download and analysis."
    exit 0
fi

if [ "$skip_download" -eq 0 ]; then
    echo "Downloading FASTQ files to $rawfq_dir_batch"
    (
        cd "$rawfq_dir_batch"
        "$PYTHON_BIN" "$download_script" "$filtered_info_csv" "$rawfq_dir_batch"
        bash download_fq.sh
    )
else
    echo "Skip download: using existing FASTQ files in $rawfq_dir_batch"
fi

cd "$working_dir_batch"
mkdir -p log

rm -f ./fq_matched.tsv
if ! "$PYTHON_BIN" "$match_sample_script" "$filtered_info_csv" "$rawfq_dir_batch" ./fq_matched.tsv; then
    echo "ERROR: match_sample.py failed - cannot match FASTQ files from $rawfq_dir_batch"
    exit 1
fi

matched_count=0
if [ -f ./fq_matched.tsv ]; then
    matched_count=$(( $(wc -l < ./fq_matched.tsv | tr -d ' ') - 1 ))
fi
if [ "$matched_count" -le 0 ]; then
    echo "ERROR: match_sample.py produced no matched FASTQ rows: $working_dir_batch/fq_matched.tsv" >&2
    exit 1
fi

rm -f get_shell.sh
"$PYTHON_BIN" "$TCR_GET_SHELL_SCRIPT" "$working_dir_batch/fq_matched.tsv" "$working_dir_batch" "$threads" "$PRESET"
if [ ! -s get_shell.sh ]; then
    echo "ERROR: tcr_get_shell_fixed_primers_v3.py produced no analysis commands: $working_dir_batch/get_shell.sh" >&2
    exit 1
fi

"$ADD_WAIT_SCRIPT" get_shell.sh 8 run2.sh
bash run2.sh > log/run2.log 2>&1
"$TCR_QC_SCRIPT" . total_result.tsv
"$FASTP_QC_SCRIPT" -i . -o fastp_qc.csv -t .json

echo "combined_stat: ${working_dir_batch}/all_combined_stat.csv"
