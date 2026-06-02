# TCR / BCR Pipeline 分析操作文档

本文档说明如何从最新下机信息表 `info.csv` 启动 TCR 或 BCR 分析流程：

```text
tcr_process_v3.sh
```

每批分析都从 `info.csv` 开始。通过 `--receptor` 参数切换 TCR / BCR，流程自动适配 chain 类型和引物文件。

## 1. 流程总览

```text
info.csv
    |
    |  tcr_process_v3.sh 按 Sed_ID / Adapter / 项目关键词筛选
    v
filtered_info.csv
    |
    |  download_fq.py -> download_fq.sh
    v
FASTQ 下载到 <rawfq-root>/<batch_id>
    |
    |  match_sample.py
    v
fq_matched.tsv
    |
    |  tcr_get_shell_fixed_primers_v3.py
    v
get_shell.sh
    |
    |  add_wait.py -> run2.sh
    v
TCR_analysis_pipeline.v7.sh
    |
    |  fastp -> FLASH -> primer filtering -> MiXCR -> VDJTools -> summary
    v
每个样本的分析结果 + 批次汇总结果
```

## 2. 脚本对应关系

| 脚本 | 作用 |
| --- | --- |
| `tcr_process_v3.sh` | 唯一主入口：筛选 `info.csv`、下载 FASTQ、匹配 R1/R2、生成任务、运行分析、汇总 QC |
| `sample_filter_info.sh` | 原始筛选脚本，当前仅保留作历史参考，筛选逻辑已被 `tcr_process_v3.sh` 吸收 |
| `match_sample.py` | 根据 `filtered_info.csv` 和 FASTQ 目录生成 `fq_matched.tsv` |
| `tcr_get_shell_fixed_primers_v3.py` | 读取 `fq_matched.tsv`，生成 `get_shell.sh` |
| `TCR_analysis_pipeline.v7.sh` | 单样本分析脚本 |
| `batch_vdjtools_presets.sh` | 独立 preset 对比脚本，非日常主流程必需 |

## 3. 运行前准备

### 3.1 准备 info.csv

`info.csv` 来自最新下机信息表。原始流程中通常由 xlsx 经 R 脚本转换得到。

`tcr_process_v3.sh` 默认沿用原 `sample_filter_info.sh` 的筛选逻辑：

```text
第 1 列  == --sed-id
第 12 列 包含 --adapter-pattern，TCR 默认 Adapter_TCR_，BCR 需显式指定（如 Adapter_BCR_）
第 8 列  包含 --project-keyword，默认 肺结节
```

筛选后输出列为：

```text
Sed_ID,Sample_ID,Lib_number,Name,Panel,Data_ID,Chain,Receptor
```

其中 `Chain` 由命令行 `--chain` 统一指定，`Receptor` 由 `--receptor` 指定。

**TCR chain 选项：**

```text
TRA / TRB / BOTH
```

**BCR chain 选项：**

```text
IGH / IGK / IGL / BCR_ALL
```

> `BOTH` 展开为 TRA + TRB；`BCR_ALL` 展开为 IGH + IGK + IGL。未来 TCR 侧会扩展 TRD、TRG 等链。

### 3.2 确认脚本目录

建议至少将以下脚本放在同一目录：

```text
tcr_process_v3.sh
tcr_get_shell_fixed_primers_v3.py
TCR_analysis_pipeline.v7.sh
```

`tcr_process_v3.sh` 会自动调用同级目录下的 `tcr_get_shell_fixed_primers_v3.py`。

`tcr_get_shell_fixed_primers_v3.py` 会自动把同级目录下的 `TCR_analysis_pipeline.v7.sh` 写入 `get_shell.sh`，不需要手动修改 `TCR_PIPELINE`。

### 3.3 确认外部依赖

需要确认以下公共路径可访问：

```text
/x03_haplox/users/donglf/common_tools/download_scripts/download_fq.py
/x03_haplox/users/donglf/haima_pipeline_tools/match_sample.py
/x03_haplox/users/donglf/common_tools/add_wait.py
/x03_haplox/users/donglf/miniconda3/bin/python
/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/tcr_qc_allchain.py
/haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/fastp_qc.v2.py
```

`TCR_analysis_pipeline.v7.sh` 还依赖：

```text
fastp
FLASH
seqkit
Java 11
MiXCR jar
VDJTools jar
V/J primer 文件
filter_cdr3_length.pl
stat_cdr3_length.pl
stat_VJ_count.pl
total_statistic.v3.pl
TCR_export_standard.py
```

### 3.4 MiXCR license 检查

可以手动激活MiXCR，license code 路径在 /haplox/users/xuliu/software/mixcr/license/milicense

```bash
mixcr activate-license
# 然后粘贴license code
```

测试 MiXCR：

```bash
JAVA11=/haplox/users/xuliu/software/java/jdk-11.0.30/bin/java
MIXCR_JAR=/haplox/users/xuliu/software/mixcr/mixcr.jar
$JAVA11 -jar $MIXCR_JAR --version
```

如果输出 `MiXCR v4.6.0`，说明 MiXCR 和 license 可用。

## 4. 主流程运行方法

如果自己新建了运行的目录，确保这些要执行的脚本有执行的权限（`chmod +x`）

进入项目目录：

```bash
cd /haplox/users/chenya/TCR_TEST
```

### 4.1 TCR 分析

先只生成 `filtered_info.csv`，不下载、不分析：

```bash
bash optim_scripts/tcr_process_v3.sh \
  -b 250129_TCR_batch01 \
  --info-csv /haplox/users/chenya/TCR_TEST/optim_scripts/info.csv \
  --sed-id 20260304_LH00348_0559_B23JWFJLT4 \
  --chain TRB \
  --prepare-only
```

检查筛选结果：

```bash
cat /haplox/users/chenya/TCR_TEST/250129_TCR_batch01/filtered_info.csv
```

正式运行，包含下载 FASTQ 和后续分析：

```bash
nohup bash optim_scripts/tcr_process_v3.sh \
  -b 250129_TCR_batch01 \
  --info-csv /haplox/users/chenya/TCR_TEST/optim_scripts/info.csv \
  --sed-id 20260304_LH00348_0559_B23JWFJLT4 \
  --chain TRB \
  -t 8 \
  -p generic-amplicon \
  > pipeline_250129.log 2>&1 &
```

如果 FASTQ 已经下载过，跳过下载，直接从本地 FASTQ 开始（参数`skip-download`）：

```bash
nohup bash optim_scripts/tcr_process_v3.sh \
  -b 250129_TCR_batch01 \
  --info-csv /haplox/users/chenya/TCR_TEST/optim_scripts/info.csv \
  --sed-id 20260304_LH00348_0559_B23JWFJLT4 \
  --chain TRB \
  --skip-download \
  -t 8 \
  -p generic-amplicon \
  > pipeline_250129.log 2>&1 &
```

### 4.2 BCR 分析

BCR 分析必须指定 `--receptor BCR` 和 `--adapter-pattern`：

```bash
# 先预览筛选结果
bash optim_scripts/tcr_process_v3.sh \
  -b 250130_BCR_batch01 \
  --info-csv /haplox/users/chenya/TCR_TEST/optim_scripts/info.csv \
  --sed-id 20260304_LH00348_0559_B23JWFJLT4 \
  --receptor BCR \
  --chain BCR_ALL \
  --adapter-pattern Adapter_BCR_ \
  --prepare-only
```

正式运行 BCR：

```bash
nohup bash optim_scripts/tcr_process_v3.sh \
  -b 250130_BCR_batch01 \
  --info-csv /haplox/users/chenya/TCR_TEST/optim_scripts/info.csv \
  --sed-id 20260304_LH00348_0559_B23JWFJLT4 \
  --receptor BCR \
  --chain BCR_ALL \
  --adapter-pattern Adapter_BCR_ \
  -t 8 \
  -p generic-amplicon \
  > pipeline_250130.log 2>&1 &
```

也可以只跑单条链，如只跑重链：

```bash
--receptor BCR --chain IGH --adapter-pattern Adapter_BCR_
```

### 4.3 通用参数

如需修改筛选条件：

```bash
--project-keyword 肺结节
--adapter-pattern Adapter_TCR_
```

查看主日志：

```bash
tail -f pipeline_250129.log
```

查看进程：

```bash
ps -ef | grep tcr_process_v3
```

## 5. 输出文件

假设：

```bash
batch_id=250129_TCR_batch01
work-root=/haplox/users/chenya/TCR_TEST
```

批次结果目录为：

```bash
/haplox/users/chenya/TCR_TEST/250129_TCR_batch01
```

关键文件：

| 路径 | 含义 |
| --- | --- |
| `filtered_info.csv` | 从 `info.csv` 筛出的本批次待分析样本 |
| `fq_matched.tsv` | 追加本地 R1/R2 FASTQ 路径后的分析输入表 |
| `get_shell.sh` | 每个样本的单样本分析命令 |
| `run2.sh` | 实际批量执行脚本 |
| `log/run2.log` | 批量任务日志 |
| `total_result.tsv` | 全批次 TCR QC 汇总 |
| `fastp_qc.csv` | fastp QC 汇总 |
| `<sample_id>/` | 单样本结果目录 |

单样本目录中常用结果（`<CHAIN>` 对 TCR 为 TRA/TRB，对 BCR 为 IGH/IGK/IGL）：

| 路径 | 含义 |
| --- | --- |
| `Map_Clone_Analysis/*.report` | MiXCR align report |
| `Map_Clone_Analysis/*.log` | MiXCR align/assemble/export 日志 |
| `Map_Clone_Analysis/*.clonotypes.<CHAIN>.raw.txt` | MiXCR 原始 clonotype 表 |
| `Map_Clone_Analysis/*.clonotypes.<CHAIN>.txt` | CDR3 长度过滤后的 clonotype 表 |
| `Map_Clone_Analysis/convert.*.txt` | VDJTools 格式文件 |
| `Stat_Picture/*.cdr3.stat.txt` | CDR3 长度统计 |
| `Stat_Picture/*.VJ.stat.txt` | VJ 使用统计 |
| `Stat_Picture/*.png` | CDR3 长度图、VJ heatmap 等 |
| `Stat_Picture/vdjtools/` | VDJTools diversity、segment、spectratype 结果 |
| `log/` | 单样本各步骤日志 |

## 6. 常见问题排查

### 6.1 filtered_info.csv 为空

表现：

```text
ERROR: no samples matched sed_id=...
```

检查：

```bash
head /path/to/info.csv
```

重点确认：

- `--sed-id` 是否等于 `info.csv` 第 1 列对应下机批次
- `--adapter-pattern` 是否能匹配第 12 列（TCR 默认 `Adapter_TCR_`，BCR 需显式指定如 `Adapter_BCR_`）
- `--project-keyword` 是否能匹配第 8 列
- `info.csv` 是否为逗号分隔 CSV
- BCR 分析是否遗漏了 `--receptor BCR` 和 `--adapter-pattern`

### 6.2 下载失败

检查：

```bash
tail -100 pipeline_250129.log
ls -lh /haplox/rawfq/TCR/<batch_id>
```

重点确认：

- `download_fq.py` 是否可访问
- `filtered_info.csv` 中 `Data_ID` 是否包含有效云端路径
- 当前用户是否有下载权限

### 6.3 BCR primer 文件未配置

BCR 分析的 primer 路径当前使用占位路径（`TCR_analysis_pipeline.v7.sh` 中 `BCR_primers/`），首次运行 BCR 前需确认并将路径替换为实际 BCR V/J primer 文件位置：

```bash
grep BCR_primers scripts/TCR_analysis_pipeline.v7.sh
```

### 6.4 fq_matched.tsv 生成失败

检查：

```bash
cat /haplox/users/chenya/TCR_TEST/<batch_id>/filtered_info.csv
ls -lh /haplox/rawfq/TCR/<batch_id>
```

常见原因：

- FASTQ 没有下载到 `<rawfq-root>/<batch_id>`
- FASTQ 文件名不以 `Sample_ID` 开头
- R1/R2 文件不成对

### 6.5 get_shell.sh 为空

检查：

```bash
wc -l /haplox/users/chenya/TCR_TEST/<batch_id>/get_shell.sh
cat /haplox/users/chenya/TCR_TEST/<batch_id>/fq_matched.tsv
```

常见原因：

- `Raw_Path_R1` 或 `Raw_Path_R2` 为空
- `Chain` 与 `--receptor` 不匹配（TCR: TRA/TRB/BOTH；BCR: IGH/IGK/IGL/BCR_ALL）
- `TCR_analysis_pipeline.v7.sh` 没有和 `tcr_get_shell_fixed_primers_v3.py` 放在同一目录

### 6.6 MiXCR license 报错

在Map_Clone_Analysis文件夹中的 .log 文件报错类似：

```text
License manager thread died.
=== No License ===
```

说明 MiXCR jar 能运行，但当前用户没有可用 license。

检查：

```bash
ls -lh ~/.mixcr/license.lic
```

### 6.7 某个样本中途失败

先看批量日志：

```bash
tail -100 /haplox/users/chenya/TCR_TEST/<batch_id>/log/run2.log
```

再进入样本目录看分步骤日志：

```bash
cd /haplox/users/chenya/TCR_TEST/<batch_id>/<sample_id>
ls log
tail -100 log/<sample_id>.fastp.log
tail -100 log/<sample_id>.flash.log
tail -100 Map_Clone_Analysis/<sample_id>.log
```

定位原则：

- fastp 失败：检查 FASTQ 文件是否损坏、路径是否正确
- FLASH 失败：检查 reads 长度、双端配对和输入文件
- primer 过滤后 reads 很少：检查引物文件或样本链型
- MiXCR 失败：检查 license、preset、输入 FASTQ 是否为空
- VDJTools 失败：检查 `vdjtools.jar` 路径和 convert 文件是否存在

## 7. 推荐每批运行清单

```text
1. 准备最新 info.csv
2. 确认 batch_id、sed_id、receptor、chain、project_keyword、adapter_pattern
3. 运行 tcr_process_v3.sh --prepare-only 生成 filtered_info.csv
4. 检查 filtered_info.csv 样本数和 Chain / Receptor
5. BCR 分析：确认 TCR_analysis_pipeline.v7.sh 中 BCR primer 路径已配置
6. 测试 MiXCR license: java -jar mixcr.jar --version
7. nohup 启动 tcr_process_v3.sh，自动下载 FASTQ 并分析
8. tail -f pipeline.log 查看主日志
9. 查看 <batch_id>/filtered_info.csv、fq_matched.tsv、log/run2.log
10. 检查 total_result.tsv、fastp_qc.csv 和各样本 Map_Clone_Analysis/Stat_Picture 输出
```

