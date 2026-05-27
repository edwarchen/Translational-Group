#!/bin/bash
set -e

# ==========================
# batch_vdjtools_presets.sh — VDJTools 批量分析 (含 OverlapPair / Rarefaction / ClusterSamples)
#
# USAGE:
#   bash batch_vdjtools_presets.sh                 # 默认过滤 singletons (readCount>1)
#   MIN_READCOUNT=1 bash batch_vdjtools_presets.sh # 不过滤（保留所有克隆）
#   MIN_READCOUNT=3 bash batch_vdjtools_presets.sh # 过滤 readCount<3
#
# 分析模块：
#   Phase 1: CalcDiversityStats / CalcSegmentUsage / CalcSpectratype (per sample)
#   Phase 2: OverlapPair (per biological sample, amp vs rna)
#   Phase 3: Rarefaction (all 10 samples combined)
#   Phase 4: ClusterSamples (all 10 samples combined)
# ==========================

# ---------- 过滤配置 ----------
MIN_READCOUNT="${MIN_READCOUNT:-2}"

# ---------- 服务器路径配置 ----------
JAVA11="/haplox/users/xuliu/software/java/jdk-11.0.30/bin/java"
VDJTOOLS_JAR="/x03_haplox/users/donglf/TCR_chenyr/shell/vdjtools-1.2.1/vdjtools-1.2.1.jar"
PYTHON3="${PYTHON3:-python3}"

# MiXCR 结果根目录
BASE_V3="/haplox/users/chenya/TCR_TEST/test_v3"
BASE_V3_RNA="/haplox/users/chenya/TCR_TEST/test_v3_rnaseq"

# VDJTools 输出根目录
if [ "$MIN_READCOUNT" -eq 1 ]; then
    VDJTOOLS_BASE="/haplox/users/chenya/TCR_TEST/vdjtools_comparison"
    FILTER_LABEL="unfiltered"
else
    VDJTOOLS_BASE="/haplox/users/chenya/TCR_TEST/vdjtools_comparison_minCount${MIN_READCOUNT}"
    FILTER_LABEL="minReadCount=${MIN_READCOUNT}"
fi

# 样本列表
SAMPLES=(
    "S044_SZ20250508032WHB-2_gdna_genome_2537265"
    "S045_SZ20250522049WHB-1_gdna_genome_2537266"
    "S046_SZ20250522052WHB-9_gdna_genome_2537267"
    "S047_SZ20250522068WHB-0_gdna_genome_2537268"
    "S048_SZ20250621039WHB-2_gdna_genome_2537269"
)

declare -A PRESET_PATHS
PRESET_PATHS["generic-amplicon"]="$BASE_V3"
PRESET_PATHS["rna-seq"]="$BASE_V3_RNA"

CHAIN="TRB"
LOGDIR="$VDJTOOLS_BASE/logs"

# ---------- 前置检查 ----------
if [ ! -f "$VDJTOOLS_JAR" ]; then
    echo "ERROR: VDJTools jar not found at $VDJTOOLS_JAR"
    exit 1
fi

if ! command -v "$PYTHON3" &>/dev/null; then
    echo "ERROR: python3 not found. Set PYTHON3 env var."
    exit 1
fi

echo "====================================================================="
echo " VDJTools Batch Analysis"
echo " Filter: $FILTER_LABEL"
echo " Output: $VDJTOOLS_BASE"
echo " Phases: CalcDiversityStats | CalcSegmentUsage | CalcSpectratype"
echo "         OverlapPair | Rarefaction | ClusterSamples"
echo "====================================================================="

# =====================================================================
# PHASE 1: Per-sample analysis
# =====================================================================
TOTAL=$(( ${#SAMPLES[@]} * 2 ))
CURRENT=0

# Track convert file paths for Phase 2-4
declare -a CONVERT_PATHS_AMP
declare -a CONVERT_PATHS_RNA

for SAMPLE in "${SAMPLES[@]}"; do
    for PRESET in "generic-amplicon" "rna-seq"; do

        CURRENT=$((CURRENT + 1))
        SRC_DIR="${PRESET_PATHS[$PRESET]}"
        RAW_FILE="$SRC_DIR/$SAMPLE/Map_Clone_Analysis/$SAMPLE.clonotypes.TRB.raw.txt"
        VDJTOOLS_OUT="$VDJTOOLS_BASE/$PRESET/$SAMPLE"
        CONVERT_FILE="$VDJTOOLS_BASE/$PRESET/$SAMPLE/$SAMPLE.TRB.convert.txt"
        METADATA="$VDJTOOLS_BASE/$PRESET/$SAMPLE/$SAMPLE.TRB.metadata.txt"

        echo ""
        echo "====================================================================="
        echo " [Phase 1: $CURRENT/$TOTAL] $SAMPLE  |  $PRESET  |  $FILTER_LABEL"
        echo "====================================================================="

        if [ ! -f "$RAW_FILE" ]; then
            echo "  SKIP: raw file not found: $RAW_FILE"
            continue
        fi

        mkdir -p "$VDJTOOLS_OUT" "$LOGDIR"

        # ---- Convert ----
        echo "  [0/3] Convert MiXCR -> VDJTools format (minReadCount=$MIN_READCOUNT) ..."
        $PYTHON3 -c "
import csv, re

raw = '$RAW_FILE'
out = '$CONVERT_FILE'
min_count = $MIN_READCOUNT

with open(raw) as f:
    reader = csv.DictReader(f, delimiter='\t')
    rows = list(reader)

kept = 0
dropped = 0
with open(out, 'w') as f:
    f.write('count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\n')
    for r in rows:
        count = int(float(r['readCount']))
        if count < min_count:
            dropped += 1
            continue
        kept += 1
        v_raw = r.get('allVHitsWithScore', '')
        v = re.sub(r'\([^)]*\)', '', v_raw.split(',')[0]).strip() if v_raw else ''
        d_raw = r.get('allDHitsWithScore', '')
        d = re.sub(r'\([^)]*\)', '', d_raw.split(',')[0]).strip() if d_raw else ''
        j_raw = r.get('allJHitsWithScore', '')
        j = re.sub(r'\([^)]*\)', '', j_raw.split(',')[0]).strip() if j_raw else ''
        f.write('\t'.join([
            str(count),
            r['readFraction'],
            r['nFeatureSequences(CDR3)'],
            r['aaFeatureSequences(CDR3)'],
            v, d, j
        ]) + '\n')

print(f'Total clones: {len(rows)}  |  Kept: {kept}  |  Dropped (readCount<{min_count}): {dropped}')
" > "$LOGDIR/$SAMPLE.$PRESET.convert.log" 2>&1
        cat "$LOGDIR/$SAMPLE.$PRESET.convert.log"

        # Store convert path for later phases
        if [ "$PRESET" == "generic-amplicon" ]; then
            CONVERT_PATHS_AMP+=("$CONVERT_FILE")
        else
            CONVERT_PATHS_RNA+=("$CONVERT_FILE")
        fi

        # ---- CalcDiversityStats ----
        printf "#filepath\tsample_id\n%s\t%s\n" "$CONVERT_FILE" "$SAMPLE" > "$METADATA"
        echo "  [1/3] CalcDiversityStats ..."
        $JAVA11 -Xmx8g -jar $VDJTOOLS_JAR CalcDiversityStats \
            -m "$METADATA" \
            "$VDJTOOLS_OUT/$SAMPLE.TRB.diversity" \
            > "$LOGDIR/$SAMPLE.$PRESET.diversity.log" 2>&1
        ls "$VDJTOOLS_OUT/$SAMPLE.TRB.diversity"* 2>/dev/null && echo "    -> OK" || echo "    -> FAILED"

        # ---- CalcSegmentUsage ----
        echo "  [2/3] CalcSegmentUsage ..."
        $JAVA11 -Xmx8g -jar $VDJTOOLS_JAR CalcSegmentUsage \
            -m "$METADATA" \
            -p "$VDJTOOLS_OUT/$SAMPLE.TRB.segment" \
            > "$LOGDIR/$SAMPLE.$PRESET.segment.log" 2>&1
        ls "$VDJTOOLS_OUT/$SAMPLE.TRB.segment"* 2>/dev/null && echo "    -> OK" || echo "    -> FAILED"

        # ---- CalcSpectratype ----
        echo "  [3/3] CalcSpectratype ..."
        $JAVA11 -Xmx8g -jar $VDJTOOLS_JAR CalcSpectratype \
            -m "$METADATA" \
            "$VDJTOOLS_OUT/$SAMPLE.TRB.spectratype" \
            > "$LOGDIR/$SAMPLE.$PRESET.spectratype.log" 2>&1
        ls "$VDJTOOLS_OUT/$SAMPLE.TRB.spectratype"* 2>/dev/null && echo "    -> OK" || echo "    -> FAILED"

        rm -f "$METADATA"

    done
done

echo ""
echo "====================================================================="
echo " Phase 1 Complete ($TOTAL samples processed)"
echo "====================================================================="


# =====================================================================
# PHASE 2: OverlapPair (amp vs rna per biological sample)
#   语法: OverlapPair [options] sample1_file sample2_file output_prefix
# =====================================================================
echo ""
echo "====================================================================="
echo " Phase 2: OverlapPair (amp vs rna per sample)"
echo "====================================================================="

OVERLAP_DIR="$VDJTOOLS_BASE/overlap"
mkdir -p "$OVERLAP_DIR"

for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    CONVERT_AMP="${CONVERT_PATHS_AMP[$i]}"
    CONVERT_RNA="${CONVERT_PATHS_RNA[$i]}"

    echo ""
    echo "  [$((i + 1))/${#SAMPLES[@]}] $SAMPLE: amp vs rna"

    # OverlapPair uses FILENAME as sample ID — need distinct names
    LINK_AMP="$OVERLAP_DIR/${SAMPLE}_amp.TRB.convert.txt"
    LINK_RNA="$OVERLAP_DIR/${SAMPLE}_rna.TRB.convert.txt"
    ln -sf "$CONVERT_AMP" "$LINK_AMP"
    ln -sf "$CONVERT_RNA" "$LINK_RNA"

    echo "  Running OverlapPair (strict intersection) ..."
    $JAVA11 -Xmx8g -jar $VDJTOOLS_JAR OverlapPair \
        -i strict \
        "$LINK_AMP" "$LINK_RNA" \
        "$OVERLAP_DIR/$SAMPLE.TRB.overlap" \
        > "$LOGDIR/$SAMPLE.overlap.log" 2>&1
    ls "$OVERLAP_DIR/$SAMPLE.TRB.overlap"* 2>/dev/null && echo "    -> OK" || echo "    -> FAILED"

    rm -f "$LINK_AMP" "$LINK_RNA"
done

echo ""
echo " Phase 2 Complete"
echo " Output: $OVERLAP_DIR/"


# =====================================================================
# PHASE 3: RarefactionPlot (all samples combined)
#   语法: RarefactionPlot [options] -m metadata output_prefix
# =====================================================================
echo ""
echo "====================================================================="
echo " Phase 3: RarefactionPlot (all 10 samples)"
echo "====================================================================="

RAREFACTION_DIR="$VDJTOOLS_BASE/rarefaction"
mkdir -p "$RAREFACTION_DIR"
RAREFACTION_META="$RAREFACTION_DIR/all.TRB.rarefaction.metadata.txt"

true > "$RAREFACTION_META"
for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    printf "#filepath\tsample_id\n%s\t%s_amp\n%s\t%s_rna\n" \
        "${CONVERT_PATHS_AMP[$i]}" "$SAMPLE" \
        "${CONVERT_PATHS_RNA[$i]}" "$SAMPLE" \
        >> "$RAREFACTION_META"
done

echo "  Running RarefactionPlot on 10 samples ..."
$JAVA11 -Xmx8g -jar $VDJTOOLS_JAR RarefactionPlot \
    -m "$RAREFACTION_META" \
    "$RAREFACTION_DIR/all.TRB.rarefaction" \
    > "$LOGDIR/rarefaction.log" 2>&1
ls "$RAREFACTION_DIR/all.TRB.rarefaction"* 2>/dev/null && echo "    -> OK" || echo "    -> FAILED"

rm -f "$RAREFACTION_META"

echo ""
echo " Phase 3 Complete"
echo " Output: $RAREFACTION_DIR/"


# =====================================================================
# PHASE 4: CalcPairwiseDistances + ClusterSamples
#   步骤 4a: CalcPairwiseDistances -m metadata output_prefix
#   步骤 4b: ClusterSamples -p input_prefix output_prefix
# =====================================================================
echo ""
echo "====================================================================="
echo " Phase 4: ClusterSamples (all 10 samples)"
echo "====================================================================="

CLUSTER_DIR="$VDJTOOLS_BASE/cluster"
mkdir -p "$CLUSTER_DIR"
CLUSTER_META="$CLUSTER_DIR/all.TRB.cluster.metadata.txt"

true > "$CLUSTER_META"
for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"
    printf "#filepath\tsample_id\n%s\t%s_amp\n%s\t%s_rna\n" \
        "${CONVERT_PATHS_AMP[$i]}" "$SAMPLE" \
        "${CONVERT_PATHS_RNA[$i]}" "$SAMPLE" \
        >> "$CLUSTER_META"
done

# Step 4a: CalcPairwiseDistances
echo "  [4a] CalcPairwiseDistances on 10 samples ..."
$JAVA11 -Xmx8g -jar $VDJTOOLS_JAR CalcPairwiseDistances \
    -m "$CLUSTER_META" \
    -i aa \
    "$CLUSTER_DIR/all.TRB.pwdist" \
    > "$LOGDIR/cluster.pwdist.log" 2>&1
ls "$CLUSTER_DIR/all.TRB.pwdist"* 2>/dev/null && echo "    -> OK" || echo "    -> FAILED"

# Step 4b: ClusterSamples (uses CalcPairwiseDistances output as input)
echo "  [4b] ClusterSamples ..."
$JAVA11 -Xmx8g -jar $VDJTOOLS_JAR ClusterSamples \
    -i aa \
    -p \
    "$CLUSTER_DIR/all.TRB.pwdist" \
    "$CLUSTER_DIR/all.TRB.cluster" \
    > "$LOGDIR/cluster.log" 2>&1
ls "$CLUSTER_DIR/all.TRB.cluster"* 2>/dev/null && echo "    -> OK" || echo "    -> FAILED"

rm -f "$CLUSTER_META"

echo ""
echo " Phase 4 Complete"
echo " Output: $CLUSTER_DIR/"


# =====================================================================
# PHASE 5: PlotFancyVJUsage (per sample, both presets)
#   语法: PlotFancyVJUsage [options] input_file output_prefix
# =====================================================================
echo ""
echo "====================================================================="
echo " Phase 5: PlotFancyVJUsage (V-J pairing heatmaps)"
echo "====================================================================="

VJUSAGE_DIR="$VDJTOOLS_BASE/vjusage"
mkdir -p "$VJUSAGE_DIR"

VJ_CURRENT=0
VJ_TOTAL=$(( ${#SAMPLES[@]} * 2 ))

for PRESET in "generic-amplicon" "rna-seq"; do
    for i in "${!SAMPLES[@]}"; do
        VJ_CURRENT=$((VJ_CURRENT + 1))
        SAMPLE="${SAMPLES[$i]}"
        CONVERT="${CONVERT_PATHS_AMP[$i]}"  # reuse amp array for paths
        if [ "$PRESET" = "rna-seq" ]; then
            CONVERT="${CONVERT_PATHS_RNA[$i]}"
        fi
        OUTPUT_PREFIX="$VJUSAGE_DIR/$SAMPLE.$PRESET.TRB"

        echo ""
        echo "  [$VJ_CURRENT/$VJ_TOTAL] $SAMPLE | $PRESET"

        $JAVA11 -Xmx8g -jar $VDJTOOLS_JAR PlotFancyVJUsage \
            --plot-type pdf \
            "$CONVERT" \
            "$OUTPUT_PREFIX" \
            > "$LOGDIR/$SAMPLE.$PRESET.vjusage.log" 2>&1

        ls "$OUTPUT_PREFIX"* 2>/dev/null && echo "    -> OK" || echo "    -> FAILED"
    done
done

echo ""
echo " Phase 5 Complete"
echo " Output: $VJUSAGE_DIR/"


# =====================================================================
# DONE
# =====================================================================
echo ""
echo "====================================================================="
echo " ALL DONE"
echo " Filter:   $FILTER_LABEL"
echo " Results:  $VDJTOOLS_BASE/"
echo "   diversity/segment/spectratype -> generic-amplicon/ & rna-seq/"
echo "   overlap       -> $OVERLAP_DIR/"
echo "   rarefaction   -> $RAREFACTION_DIR/"
echo "   cluster       -> $CLUSTER_DIR/"
echo "   vjusage       -> $VJUSAGE_DIR/"
echo " Logs:     $LOGDIR/"
echo "====================================================================="
