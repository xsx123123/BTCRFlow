# BTCRFlow 框架说明

本文档描述 BTCRFlow（Bulk BCR/TCR Repertoire Analysis Pipeline）的实现结构。
BTCRFlow 参考 SCRNAFlow 的 Flow 家族框架构建，将原有 `07.bulkBCR_v4` 和 `08.bulkTCR_v4` 两个独立流程统一为单一可配置流程。

---

## 1. 顶层职责

BTCRFlow 的边界：

- 从 `analysisyaml` 读取项目配置、原始数据目录、工作目录和样本表。
- 使用编译好的 `seq_preprocessor` 标准化原始 FASTQ 命名，生成 `raw_data_md5.json`。
- 使用编译好的 `json_md5_verifier` 验证 MD5，生成 `01.qc/md5_check.tsv`（硬性门控）。
- 对标准化后的 FASTQ 做 FastQC/MultiQC 质控（依赖 MD5 校验通过）。
- 使用 cutadapt 进行 Clontech SMARTer 协议 trimming（-U 30）。
- 按平台（BCR/TCR）选择 MiXCR preset 运行 `mixcr analyze`。
- 可选 VDJtools 下游分析（downsample、diversity、V-J pairing）。
- 可选交付到 `data_deliver/`，可选生成 `project_summary.json`。

不在当前边界内：单细胞 10x/Cell Ranger、DNB、STAR/RSEM mapping、DEG、Seurat/Scanpy。

---

## 2. 目录结构

```text
BTCRFlow/
├── snakefile
├── config/
│   ├── config.yaml
│   ├── reference.yaml
│   ├── run_parameter.yaml
│   ├── cluster_config.yaml
│   ├── samples.csv
│   └── contrasts.csv
├── schema/
│   └── config.schema.yaml
├── rules/
│   ├── 01.common.smk
│   ├── 02.file_convert_md5.smk
│   ├── 03.short_read_qc.smk
│   ├── 04.trimming.smk
│   ├── 05.mixcr.smk
│   ├── 06.vdj_analysis.smk
│   ├── 07.deliver.smk
│   ├── 08.report.smk
│   └── utils/
│       ├── __init__.py
│       ├── common.py
│       ├── datadeliver.py
│       └── tools.py
├── src/src/md5/Compile_release/
│   ├── seq_preprocessor_x86_64/release/seq_preprocessor
│   └── json_md5_verifier_x86_64/release/json_md5_verifier
├── scripts/
│   ├── prepare_btcr_links.py
│   ├── btcrflow_deliver.py
│   ├── spectratype.R
│   ├── vj_pairing_plot.R
│   ├── pieplot_top20.R
│   └── Vgene_Freq.pl
├── envs/
│   ├── qc.yaml
│   ├── trimming.yaml
│   ├── mixcr.yaml
│   ├── vdjtools.yaml
│   └── report.yaml
├── docs/
└── report/
```

---

## 3. 入口 snakefile

保持 Flow 家族配置优先级：

```bash
snakemake -s /path/to/BTCRFlow/snakefile \
  --cores 32 \
  --use-conda \
  --config analysisyaml=/project/01.workflow/config.yaml
```

流程：
1. `min_version("9.9.0")`
2. 加载四件套默认配置
3. `load_user_config(analysisyaml)` 最高优先级覆盖
4. `validate(config, schema)`
5. `workdir: config["workflow"]`
6. `load_samples()` 要求 `sample,sample_name,group,platform,path`
7. 按 platform 分为 BCR_SAMPLES / TCR_SAMPLES
8. `DataDeliver()` 动态生成 rule all 目标

---

## 4. 规则链路

```text
01.common.smk              # 导入 utils，无业务 rule
02.file_convert_md5.smk    # seq_preprocessor 标准化命名 -> json_md5_verifier MD5 校验
03.short_read_qc.smk       # FastQC + MultiQC（依赖 md5_check.tsv 门控）
04.trimming.smk            # cutadapt (Clontech -U 30)
05.mixcr.smk               # MiXCR analyze + postanalysis
06.vdj_analysis.smk        # downsample + VDJtools diversity
07.deliver.smk             # 交付 manifest + MD5
08.report.smk              # project_summary.json
```

---

## 5. 样本表约定

必填列：

```csv
sample,sample_name,group,platform,path,index
```

- `platform`: `bcr` 或 `tcr`（也接受 `bulk_bcr`、`bulk_tcr`）
- `path`: 原始 FASTQ 所在目录
- `index`: 文件前缀（用于 glob 匹配 R1/R2）

BCR 样本使用 MiXCR preset: `takara-human-rna-bcr-umi-smarter`
TCR 样本使用 MiXCR preset: `generic-amplicon`

---

## 6. 数据流

```text
analysisyaml
  ├── raw_data_path / sample_csv
  ├── workflow
  └── data_deliver
        │
        ▼
02.file_convert_md5
  ├── seq_preprocessor
  │     └── 00.raw_data/link_dir/raw_data_md5.json
  │     └── 00.raw_data/link_dir/{sample}/{sample}_R1/R2.fq.gz
  └── json_md5_verifier
        └── 01.qc/md5_check.tsv  ← 硬性门控
        │
        ├──> 03.short_read_qc (input: md5_check.tsv)
        │      └── 01.qc/multiqc/multiqc.html
        │
        ├──> 04.trimming (input: md5_check.tsv)
        │      └── 02.trimmed/{sample}_trimmed_R1/R2.fastq.gz
        │
        └──> 05.mixcr (input: md5_check.tsv + link_dir FASTQ)
               ├── 03.mixcr/{sample}.clns
               └── 03.mixcr/postanalysis.json.gz
                      │
                      ▼
               06.vdj_analysis
                 ├── 04.vdj_analysis/downsampled/{sample}.{chain}.downsampled.tsv
                 ├── 04.vdj_analysis/converted/{sample}/
                 └── 04.vdj_analysis/diversity/{sample}.diversity.done
                      │
                      ▼
               07.deliver / 08.report
                 ├── data_deliver/delivery_manifest.json
                 └── 05.report/project_summary.json
```

---

## 7. 与旧版 bulkBCR_v4 / bulkTCR_v4 的差异

| 旧版 | BTCRFlow | 说明 |
|---|---|---|
| 两个独立流程，各自有 workflow/ 入口 | 统一 snakefile，platform 列区分 | 一套配置管理 BCR+TCR |
| 硬编码工具路径（/nextcode/, /glusterfs/） | config.software 配置化 + conda env | 可移植性大幅提升 |
| essential_config.yaml + config.yaml 双文件 | Flow 家族四件套 + analysisyaml 覆盖 | 配置优先级清晰 |
| 无 schema 验证（或 draft-04） | JSON Schema draft-07 验证 | 配置错误提前暴露 |
| vdjtools 全链路（basic/overlap/clonetype/annotation） | 精简为 downsample + diversity | 按需扩展 |
| 手动 shell 脚本做下游图表 | R 脚本参数化，可集成到 rule | 可重复性提升 |
| SGE cluster.json | cluster_config.yaml（通用） | 适配多种调度器 |
| dataDelivery/ 平铺输出 | 编号目录（00-05） | 与 SCRNAFlow 一致 |

---

## 8. 模块开关

```yaml
qc_clean: true      # FastQC + MultiQC
trimming: true      # cutadapt
mixcr: true         # MiXCR analyze + postanalysis
vdj_analysis: true  # VDJtools downsample + diversity
deliver: true       # 交付 manifest
report: true        # project_summary.json
```

模块依赖：`mixcr` 依赖 `qc_clean` + `trimming`；`vdj_analysis` 依赖 `mixcr`。

---

## 9. 验证命令

```bash
# Dry-run
snakemake -s snakefile -n --cores 1 \
  --config analysisyaml=/path/to/config.yaml

# 只运行到 MiXCR
snakemake -s snakefile --cores 12 --use-conda \
  --config analysisyaml=/path/to/config.yaml \
  03.mixcr/postanalysis.json.gz

# 完整运行
snakemake -s snakefile --cores 32 --use-conda \
  --config analysisyaml=/path/to/config.yaml
```

---

## 10. 后续优化 Checklist

1. 为 `06.vdj_analysis.smk` 补充 BCR 链（IGH/IGK/IGL）的完整 downsample 分支。
2. 将 `spectratype.R`、`vj_pairing_plot.R` 集成为 Snakemake rule，而非手动运行。
3. 补充 `mixcr exportClones` 的完整 clone table 输出（当前只有 downsampled TSV）。
4. 考虑容器化 MiXCR（类似 SCRNAFlow 对 Cell Ranger 的策略）。
5. 为 `prepare_btcr_links.py` 增加单元测试。
6. 报告模块接入正式 R Markdown / Quarto 模板。
7. 支持多 lane 样本合并（当前假设单 lane）。
