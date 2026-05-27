#!/x03_haplox/users/donglf/miniconda3/bin/python
# coding: utf-8
# RNA测序TCR分析 pipeline v3.1，删除了对V、J的过滤，适用于TCR不同链的分析，增加了对链类型的兼容性（Chain / chain），并且对链类型进行了标准化处理（大写）。使用时需要提供样本信息文件和工作目录，以及线程数。样本信息文件需要包含 Sample_ID（或 Lib_number）和 Chain（或 chain）列，以及 Raw_Path_R1 和 Raw_Path_R2 列。脚本会根据样本信息生成一个 get_shell.sh 脚本，里面包含了每个样本的分析命令，可以直接运行该脚本进行分析。
# usage: /haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/tcr_get_shell_fixed_primers_v2.py  <sample_info_file> <working_dir>
# example working dir: /x05_haplox/users/donglf/Project/TCR/220809
# example: /haplox/users/xuliu/TCR_Project/scripts/TCR_rna_pipeline/tcr_get_shell_fixed_primers_v2.py /x05_haplox/rawfq/TCR_2022/220809/fq_matched.tsv /x05_haplox/users/donglf/Project/TCR/220809 10 


import pandas as pd
import sys
import os
import utils

TCR_PIPELINE = '/haplox/users/chenya/TCR_TEST/scripts/TCR_analysis_pipeline.v7.sh'

sample_info_table = pd.read_csv(sys.argv[1], sep='\t')
cols = sample_info_table.columns.tolist()
if 'Sample_ID' not in cols and 'Lib_number' in cols:
    sample_info_table['Sample_ID'] = sample_info_table['Lib_number']

working_dir = os.path.abspath(sys.argv[2])
THREAD = int(sys.argv[3])
PRESET = sys.argv[4] if len(sys.argv) > 4 else 'generic-amplicon'

if PRESET not in ('generic-amplicon', 'rna-seq'):
    print(f"ERROR: preset must be 'generic-amplicon' or 'rna-seq', got '{PRESET}'")
    sys.exit(1)

print(f"Pipeline: {TCR_PIPELINE}")
print(f"Preset:   {PRESET}")
print(f"Threads:  {THREAD}")
all_get_shell_script = os.path.join(working_dir, 'get_shell.sh')
all_get_shell_script_fh = open(all_get_shell_script, 'w')

for _, row in sample_info_table.iterrows():
    sample_id = row['Sample_ID']
    sample_dir = os.path.join(working_dir, sample_id)

    # ==========================
    # 兼容 Chain / chain
    # ==========================
    if 'Chain' in sample_info_table.columns:
        chain = row['Chain']
    else:
        chain = row['chain']

    # ==========================
    # 标准化（关键）
    # ==========================
    chain = str(chain).strip().upper()

    if chain == 'BOTH':
        chain = 'BOTH'
    elif chain == 'TRA':
        chain = 'TRA'
    elif chain == 'TRB':
        chain = 'TRB'
    else:
        print(f"WARNING: {sample_id} invalid chain '{chain}', skip")
        continue

    os.makedirs(sample_dir, exist_ok=True)

    if utils.check_file(row['Raw_Path_R1']) and utils.check_file(row['Raw_Path_R2']):
        all_get_shell_script_fh.write(
            f'{TCR_PIPELINE} {row["Raw_Path_R1"]} {row["Raw_Path_R2"]} '
            f'{sample_dir} {sample_id} {THREAD} {chain} {PRESET}\n'
        )

all_get_shell_script_fh.close()