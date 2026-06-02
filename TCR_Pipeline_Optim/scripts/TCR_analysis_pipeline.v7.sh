#!/bin/bash
set -e

# ==========================
# TCR Analysis Pipeline v7
# v6 基础上新增 preset 参数:
#   generic-amplicon → MiXCR align --preset generic-amplicon --dna + floating boundaries
#   rna-seq          → MiXCR align --preset rna-seq
# ==========================

# ==========================
# 参数检查
# ==========================
if [ $# -lt 7 ]; then
echo ""
echo "USAGE:"
echo " $0 input_R1.fastq.gz input_R2.fastq.gz outputdir sample_key threads chain_type preset [receptor]"
echo ""
echo " chain_type: TRA / TRB / BOTH / IGH / IGK / IGL / BCR_ALL"
echo " preset:     generic-amplicon / rna-seq"
echo " receptor:   TCR (default) / BCR"
echo ""
echo "VDJTools output is written to Stat_Picture/vdjtools/"
echo ""
exit 1
fi

R1=$1
R2=$2
outdir=$3
key=$4
thread=$5
CHAIN=$6
PRESET=$7
RECEPTOR=${8:-TCR}

# ==========================
# 验证 preset 参数
# ==========================
if [ "$PRESET" != "generic-amplicon" ] && [ "$PRESET" != "rna-seq" ]; then
    echo "ERROR: preset must be 'generic-amplicon' or 'rna-seq', got '$PRESET'"
    exit 1
fi

# ==========================
# 验证 receptor 参数
# ==========================
RECEPTOR="$(echo "$RECEPTOR" | tr '[:lower:]' '[:upper:]')"
if [ "$RECEPTOR" != "TCR" ] && [ "$RECEPTOR" != "BCR" ]; then
    echo "ERROR: receptor must be TCR or BCR, got '$RECEPTOR'"
    exit 1
fi

# ==========================
# 根据 preset 设置 MiXCR align 参数
# ==========================
if [ "$PRESET" == "generic-amplicon" ]; then
    ALIGN_EXTRA="--dna --floating-left-alignment-boundary --floating-right-alignment-boundary J"
else
    ALIGN_EXTRA=""
fi

# ==========================
# 基础路径
# ==========================
Bin="/x03_haplox/users/donglf/TCR_chenyr/shell"
Bin2="/x03_haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline"

if [ "$RECEPTOR" == "BCR" ]; then
    V_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/BCR_primers/V_primer.txt"
    J_rc_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/BCR_primers/J_rc_primer.txt"
else
    V_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/V_primer.txt"
    J_rc_primers="/x03_haplox/users/donglf/tcr_scripts/total_primers/V10_primers/J_rc_primer.txt"
fi
Mixcr_Jar="/haplox/users/xuliu/software/mixcr/mixcr.jar"
JAVA11="/haplox/users/xuliu/software/java/jdk-11.0.30/bin/java"

# VDJTools — modify this path to match your installation
VDJTOOLS_JAR="/haplox/users/xuliu/software/vdjtools/vdjtools.jar"

# ==========================
# 创建目录
# ==========================
mkdir -p $outdir/{Cleanfq,Merge_PE,Map_Clone_Analysis,Stat_Picture/vdjtools,log}

# ==========================
# chain控制
# ==========================
CHAINS_TO_RUN=()

if [ "$RECEPTOR" == "BCR" ]; then
    case "$CHAIN" in
        IGH)   CHAINS_TO_RUN=("IGH") ;;
        IGK)   CHAINS_TO_RUN=("IGK") ;;
        IGL)   CHAINS_TO_RUN=("IGL") ;;
        BCR_ALL) CHAINS_TO_RUN=("IGH" "IGK" "IGL") ;;
        *)
            echo "ERROR: chain_type must be IGH / IGK / IGL / BCR_ALL for BCR"
            exit 1
            ;;
    esac
else
    case "$CHAIN" in
        TRB)  CHAINS_TO_RUN=("TRB") ;;
        TRA)  CHAINS_TO_RUN=("TRA") ;;
        BOTH) CHAINS_TO_RUN=("TRA" "TRB") ;;
        *)
            echo "ERROR: chain_type must be TRA / TRB / BOTH for TCR"
            exit 1
            ;;
    esac
fi

echo "====================================================================="
echo " MiXCR version: 4.6.0, preset: $PRESET, receptor: $RECEPTOR"
echo " VDJTools analysis: enabled"
echo " Sample: $key    Chain: $CHAIN    Threads: $thread"
echo "====================================================================="

# ==========================
# fastp
# ==========================
echo "[1/9] fastp — quality trimming & adapter removal"
$Bin/fastp-0.19.7/fastp \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -i $R1 -I $R2 \
    -o $outdir/Cleanfq/$key.good_R1.fq.gz \
    -O $outdir/Cleanfq/$key.good_R2.fq.gz \
    -t 1 -T 1 -w $thread \
    -j $outdir/Cleanfq/$key.json \
    -h $outdir/Cleanfq/$key.html \
    > $outdir/log/$key.fastp.log 2>&1

# ==========================
# FLASH
# ==========================
echo "[2/9] FLASH — paired-end merging"
$Bin/FLASH-1.2.11/FLASH-1.2.11-Linux-x86_64/flash \
    $outdir/Cleanfq/$key.good_R1.fq.gz \
    $outdir/Cleanfq/$key.good_R2.fq.gz \
    -d $outdir/Merge_PE -o $key \
    -m 10 -M 150 -p 33 -r 150 -x 0.1 \
    > $outdir/log/$key.flash.log

# ==========================
# primer筛选
# ==========================
echo "[3/9] seqkit — V-primer filtering"
/x03_haplox/users/donglf/miniconda3/envs/tcr_hapyun/bin/seqkit grep -s -i --degenerate \
    -f $V_primers \
    $outdir/Merge_PE/$key.extendedFrags.fastq \
    -o $outdir/Merge_PE/$key.extendedFrags.containV.fastq

perl $Bin/stat_primer.pl \
    -i $outdir/Merge_PE/$key.extendedFrags.fastq \
    -v $outdir/Merge_PE/$key.extendedFrags.containV.fastq \
    -o $outdir/Merge_PE/$key.primer.stat.txt

echo "[4/9] seqkit — J-primer filtering"
/x03_haplox/users/donglf/miniconda3/envs/tcr_hapyun/bin/seqkit grep -s -i --degenerate \
    -f $J_rc_primers \
    $outdir/Merge_PE/$key.extendedFrags.containV.fastq \
    -o $outdir/Merge_PE/$key.extendedFrags.containVJ.fastq

perl $Bin/stat_primer.pl \
    -i $outdir/Merge_PE/$key.extendedFrags.containV.fastq \
    -v $outdir/Merge_PE/$key.extendedFrags.containVJ.fastq \
    -o $outdir/Merge_PE/$key.primer.stat.txt

# ==========================
# MiXCR
# ==========================
echo "[5/9] MiXCR align — preset=$PRESET"
$JAVA11 -Xmx16g -Xms4g -jar $Mixcr_Jar align \
    --preset $PRESET --species hs $ALIGN_EXTRA \
    --report $outdir/Map_Clone_Analysis/$key.report \
    $outdir/Merge_PE/$key.extendedFrags.containVJ.fastq \
    $outdir/Map_Clone_Analysis/$key.vdjca \
    -t $thread \
    >$outdir/Map_Clone_Analysis/$key.log 2>&1

echo "[6/9] MiXCR assemble — clonotype clustering"
$JAVA11 -Xmx16g -Xms4g -jar $Mixcr_Jar assemble \
    --write-alignments \
    -OassemblingFeatures="[CDR3]" \
    -OseparateByV=false -OseparateByJ=false -OseparateByC=false \
    $outdir/Map_Clone_Analysis/$key.vdjca \
    $outdir/Map_Clone_Analysis/$key.clna \
    >>$outdir/Map_Clone_Analysis/$key.log 2>&1

# ==========================
# export + 下游分析（按链循环）
# ==========================
for c in "${CHAINS_TO_RUN[@]}"; do

echo "====================================================================="
echo " Processing $c ..."
echo "====================================================================="

echo "[7/9] MiXCR exportClones — $c"
$JAVA11 -Xmx16g -Xms4g -jar $Mixcr_Jar exportClones \
    --chains $c \
    $outdir/Map_Clone_Analysis/$key.clna \
    $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.tsv \
    >>$outdir/Map_Clone_Analysis/$key.log 2>&1

# 修复表头并改名为下游脚本需要的 .raw.txt
sed -i '1s/nSeqCDR3/nFeatureSequences(CDR3)/g; 1s/aaSeqCDR3/aaFeatureSequences(CDR3)/g; 1s/minQualCDR3/minFeatureQualities(CDR3)/g' \
    $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.tsv

mv $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.tsv \
   $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.raw.txt

# CDR3长度过滤
perl $Bin/filter_cdr3_length.pl \
    -i $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.raw.txt \
    -o $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.txt \
    > $outdir/log/$key.cdr3.filter.${c}.log 2>&1

# CDR3长度统计
perl $Bin/stat_cdr3_length.pl \
    -i $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.txt \
    -o $outdir/Stat_Picture/$key.clonotypes.${c}.cdr3.stat.txt \
    -p $outdir/Stat_Picture/$key.cdr3.len.${c}.png \
    >$outdir/log/$key.cdr3.len.${c}.log 2>&1

# VDJtools格式转换
/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/TCR_export_standard.py \
    $outdir/Map_Clone_Analysis/$key.clonotypes.${c}.txt \
    $outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt \
    > $outdir/log/$key.convert.${c}.log 2>&1

# VJ配对统计
perl $Bin/stat_VJ_count.pl \
    -i $outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt \
    -o $outdir/Stat_Picture/$key.${c}.VJ.stat.txt \
    -p $outdir/Stat_Picture/$key.${c}.VJ.heatmap.png \
    >$outdir/log/$key.VJ.${c}.log 2>&1

# ==========================
# VDJTools 分析
# ==========================
echo "[8/9] VDJTools — $c"

if [ -f "$VDJTOOLS_JAR" ]; then

VDJTOOLS_OUT="$outdir/Stat_Picture/vdjtools"
METADATA="$VDJTOOLS_OUT/$key.${c}.metadata.txt"

# 生成 VDJTools metadata 文件
printf "#filepath\tsample_id\n%s\t%s\n" \
    "$outdir/Map_Clone_Analysis/convert.$key.clonotypes.${c}.txt" \
    "$key" > "$METADATA"

# 8a. CalcDiversityStats — 多样性指数
echo "      8a: CalcDiversityStats (Shannon, Simpson, clonality ...)"
$JAVA11 -Xmx8g -jar $VDJTOOLS_JAR CalcDiversityStats \
    -m "$METADATA" \
    -o "$VDJTOOLS_OUT/$key.${c}.diversity" \
    > $outdir/log/$key.vdjtools.diversity.${c}.log 2>&1

# 8b. CalcSegmentUsage — V/J 基因片段使用频率
echo "      8b: CalcSegmentUsage (V/J gene segment usage)"
$JAVA11 -Xmx8g -jar $VDJTOOLS_JAR CalcSegmentUsage \
    -m "$METADATA" \
    -p -n "$VDJTOOLS_OUT/$key.${c}.segment" \
    > $outdir/log/$key.vdjtools.segment.${c}.log 2>&1

# 8c. CalcSpectratype — CDR3 长度分布 (by V gene)
echo "      8c: CalcSpectratype (CDR3 length by V gene)"
$JAVA11 -Xmx8g -jar $VDJTOOLS_JAR CalcSpectratype \
    -m "$METADATA" \
    -o "$VDJTOOLS_OUT/$key.${c}.spectratype" \
    > $outdir/log/$key.vdjtools.spectratype.${c}.log 2>&1

# 清理临时 metadata
rm -f "$METADATA"

echo "      VDJTools output -> $VDJTOOLS_OUT/"

else
echo "      [WARNING] VDJTools jar not found at $VDJTOOLS_JAR — skipping"
fi

done

# ==========================
# 总体统计
# ==========================
echo "[9/9] total statistic"
perl $Bin2/total_statistic.v3.pl \
    -i $outdir/ \
    -o $outdir/Stat_Picture \
    -k $key -c $CHAIN \
    >$outdir/log/$key.total_stat.log 2>&1

# ==========================
# 清理中间文件
# ==========================
rm -rf $outdir/Cleanfq/*fq.gz $outdir/Merge_PE/*fastq

echo ""
echo "====================================================================="
echo " Done: $key  (preset: $PRESET, receptor: $RECEPTOR)"
if [ -f "$VDJTOOLS_JAR" ]; then
echo " VDJTools results: $outdir/Stat_Picture/vdjtools/"
echo "   - diversity  : *.diversity.txt"
echo "   - segment    : *.segment.V.txt / *.segment.J.txt (heatmap PNGs)"
echo "   - spectratype: *.spectratype.txt (per V-gene CDR3 length)"
fi
echo "====================================================================="
echo ""
