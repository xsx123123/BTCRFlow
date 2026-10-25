# Flow Framework — 多组学通用流程框架

本文档描述 Flow 家族流程的通用架构规范。所有 Flow 系列流程（RNAFlow、SCRNAFlow、BTCRFlow、ATACFlow、WGSFlow 等）共享同一套入口、配置、数据完整性校验、样本表、交付和报告骨架，仅在分析层（mapping/count/variant 等）按组学场景分化。

---

## 1. 框架总览

### 1.1 设计哲学

- **共享骨架，分化分析**：入口配置、FASTQ 标准化、MD5 校验、QC、交付、报告为通用层；mapping/count/variant 为组学专用层。
- **配置优先级链**：默认四件套 → 项目级 `analysisyaml` 覆盖 → 命令行 `--config` 最终覆盖。
- **硬性数据门控**：MD5 校验不通过则流程中止，下游分析不启动。
- **模块开关**：每个分析模块可独立启停，`DataDeliver()` 动态生成目标。
- **零拷贝标准化**：原始 FASTQ 通过 symlink 标准化命名，不消耗额外磁盘。

### 1.2 Flow 家族成员

| 流程 | 组学场景 | 核心分析层 | 状态 |
|---|---|---|---|
| RNAFlow | Bulk RNA-seq | STAR mapping → RSEM count → DEG → variant → rMATS → fusion | 生产 |
| SCRNAFlow | 单细胞 RNA | Cell Ranger / dnbc4tools → count matrix | 生产 |
| BTCRFlow | Bulk BCR/TCR | cutadapt → MiXCR → VDJtools diversity | 生产 |
| ATACFlow | 染色质可及性 | BWA → MACS2 peak → annotation | 规划 |
| WGSFlow | 全基因组 | BWA → GATK → CNV → SV | 规划 |

### 1.3 通用目录结构

所有 Flow 流程遵循统一的目录布局：

```text
<X>Flow/
├── snakefile                        # 入口（小写）
├── config/
│   ├── config.yaml                  # 默认配置（software、raw_data、模块开关）
│   ├── reference.yaml               # 参考基因组路径
│   ├── run_parameter.yaml           # 运行参数（线程、内存、工具参数）
│   ├── cluster_config.yaml          # HPC 资源配置
│   ├── samples.csv                  # 样本表
│   └── contrasts.csv                # 对比表（可选）
├── schema/
│   └── config.schema.yaml           # JSON Schema 验证
├── rules/
│   ├── 01.common.smk                # 导入 utils，无业务 rule
│   ├── 02.file_convert_md5.smk      # FASTQ 标准化 + MD5 校验（通用）
│   ├── 03.short_read_qc.smk         # FastQC + MultiQC（通用）
│   ├── 04.*.smk                     # 组学专用（trimming/contamination/...）
│   ├── 05.*.smk                     # 组学专用（clean/merge_qc/...）
│   ├── 06.*.smk                     # 核心分析（mapping/count/mixcr/...）
│   ├── ...
│   ├── N-1.deliver.smk              # 交付（通用）
│   ├── N.Report.smk                 # 报告（通用）
│   └── utils/
│       ├── __init__.py
│       ├── common.py                # DataDeliver、模块依赖、辅助函数
│       ├── datadeliver.py           # 模块目标文件收集器
│       ├── id_convert.py            # 样本表/对比表解析
│       ├── validate.py              # 配置验证、参考路径检查
│       ├── reference_update.py      # resolve_reference_paths
│       ├── resource_manager.py      # HPC 资源 profile
│       └── tools.py                 # 命令行构建辅助
├── src/
│   └── JZ_Tools/                    # Git submodule（Rust 工具源码）
│       └── src/md5/
│           ├── seq_preprocessor/
│           └── json_md5_verifier/
├── scripts/                         # 组学专用脚本
├── envs/                            # Conda 环境 YAML
├── docs/
├── report/
└── container_env/                   # 容器配置（可选）
```

---

## 2. 入口 snakefile（通用）

所有 Flow 流程的 `snakefile` 遵循相同的初始化序列：

```python
from snakemake.utils import validate, min_version
min_version("9.9.0")

# 1. 加载默认四件套
configfile: "config/config.yaml"
configfile: "config/reference.yaml"
configfile: "config/run_parameter.yaml"
configfile: "config/cluster_config.yaml"

# 2. 用户 analysisyaml 最高优先级覆盖
load_user_config(config, cmd_arg_name="analysisyaml")

# 3. 解析参考基因组路径
resolve_reference_paths(config, index_keys=[...])

# 4. Schema 验证
validate(config, "schema/config.schema.yaml")
validate_genome_version(config)
validate_species(config)

# 5. 工作目录重定向
workdir: config["workflow"]

# 6. 加载样本表
SAMPLES = load_samples(config, required_columns=[...])
CONTRASTS = load_contrasts(config)

# 7. 动态目标生成
rule all:
    input: DataDeliver(config, SAMPLES, ...)

# 8. 包含规则文件
include: "rules/01.common.smk"
include: "rules/02.file_convert_md5.smk"
...
```

运行方式（所有 Flow 统一）：

```bash
snakemake -s /path/to/<X>Flow/snakefile \
  --cores 32 \
  --use-conda \
  --config analysisyaml=/project/01.workflow/config.yaml
```

---

## 3. 配置系统（通用）

### 3.1 配置优先级

```text
config/config.yaml (仓库默认)
  ← config/reference.yaml
  ← config/run_parameter.yaml
  ← config/cluster_config.yaml
  ← analysisyaml (项目级，最高优先级)
  ← --config (命令行，最终覆盖)
```

### 3.2 项目级 analysisyaml 最小结构

```yaml
project_name: my_project
Genome_Version: hg38
species: Homo_sapiens
client: demo

raw_data_path:
  - /data/project/00.raw_data
sample_csv: /data/project/01.workflow/samples.csv
paired_csv: /data/project/01.workflow/contrasts.csv
workflow: /data/project/01.workflow
data_deliver: /data/project/02.data_deliver

reference_path: /data/reference/<X>Flow_reference

# 模块开关（按需覆盖）
qc_clean: true
deliver: true
report: true
```

### 3.3 通用 config/config.yaml 关键块

```yaml
software:
  seq_preprocessor: ../src/src/md5/Compile_release/seq_preprocessor_x86_64/release/seq_preprocessor
  json_md5_verifier: ../src/src/md5/Compile_release/json_md5_verifier_x86_64/release/json_md5_verifier
  fastqc: fastqc
  multiqc: multiqc
  # ... 组学专用工具

raw_data:
  md5: md5.txt              # 原始数据目录中的 MD5 文件名
  library_type: short-read  # 传给 seq_preprocessor 的库类型
  sample_sheet_rename: false

convert_md5: link_dir       # 标准化链接输出子目录名
r1_suffix: _R1.fq.gz       # 标准化 R1 后缀
r2_suffix: _R2.fq.gz       # 标准化 R2 后缀
```

### 3.4 参考基因组

通过 `resolve_reference_paths(config, index_keys=[...])` 解析，不同流程传不同 `index_keys`：

| 流程 | index_keys |
|---|---|
| RNAFlow | `["STAR_index", "RSEM_index", "GTF", "genome_fasta"]` |
| SCRNAFlow | `["CellRanger_reference", "DNB_reference"]` |
| BTCRFlow | `["mixcr_species", "vdjtools_jar"]` |
| ATACFlow | `["BWA_index", "genome_fasta", "chrsize"]` |

---

## 4. 02 步：FASTQ 标准化与 MD5 校验（通用核心）

### 4.1 数据流

```text
raw_data_path (任意命名的原始 FASTQ)
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  seq_preprocessor (Rust)                                     │
│                                                              │
│  1. 递归扫描 raw_data_path 目录                              │
│  2. 识别多种 FASTQ 命名模式：                                │
│     - Illumina: {sample}_S1_L001_R1_001.fastq.gz            │
│     - Generic:  {sample}_R1.fq.gz / {sample}_1.fq.gz        │
│     - SRA PE:   SRRxxxxxx_1.fq.gz / _2.fq.gz               │
│     - .raw 后缀: {sample}.R1.raw.fq.gz                      │
│     - Long-read SE: {sample}.fq.gz                          │
│  3. 创建标准化软链接（零拷贝）：                             │
│     00.raw_data/link_dir/{sample}/{sample}_R1.fq.gz         │
│     00.raw_data/link_dir/{sample}/{sample}_R2.fq.gz         │
│  4. 解析原始目录中的 md5.txt，写入 per-sample md5.txt        │
│  5. 生成审计报告：                                           │
│     00.raw_data/link_dir/raw_data_md5.json                   │
└─────────────────────────────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────┐
│  json_md5_verifier (Rust)                                    │
│                                                              │
│  1. 读取 raw_data_md5.json                                   │
│  2. 多线程并行计算实际 MD5（Rayon，默认 40 线程）            │
│  3. 逐文件比对 expected vs actual MD5                        │
│  4. SRA 数据（无厂商 MD5）：计算并记录，标记 PASS            │
│  5. 输出校验报告：01.qc/md5_check.tsv                        │
│  6. 任何文件校验失败 → 非零退出码 → 流程中止                │
└─────────────────────────────────────────────────────────────┘
        │
        ▼
  01.qc/md5_check.tsv  ← 硬性门控（Hard Gate）
        │
        ├──> 03.short_read_qc (FastQC/MultiQC)
        └──> 组学专用下游规则
```

### 4.2 设计原则

- **命名标准化**：无论输入是何种命名 convention，统一为 `{sample}/{sample}_R1.fq.gz`。
- **零拷贝**：使用 Unix symlink，不消耗额外磁盘空间。
- **硬性门控**：`01.qc/md5_check.tsv` 是所有下游 rule 的显式 input 依赖，MD5 不通过则流程中止。
- **审计追踪**：`raw_data_md5.json` 保留原始路径（绝对）→ 标准化路径（相对）的完整映射。
- **SRA 兼容**：公共数据库下载的数据无厂商 MD5，工具自动计算并记录，建立基线。

### 4.3 raw_data_md5.json 格式

```json
[{
  "sample_name": "SampleA",
  "library_type": "PE",
  "new_r1_path_relative": "SampleA/SampleA_R1.fq.gz",
  "original_r1_path_absolute": "/data/raw/SampleA/SampleA_R1.clean.fq.gz",
  "md5_r1": "d41d8cd98f00b204e9800998ecf8427e",
  "new_r2_path_relative": "SampleA/SampleA_R2.fq.gz",
  "original_r2_path_absolute": "/data/raw/SampleA/SampleA_R2.clean.fq.gz",
  "md5_r2": "098f6bcd4621d373cade4e832627b4f6"
}]
```

### 4.4 md5_check.tsv 格式

```text
CheckTime	SampleName	FilePath	ExpectedMD5	ActualMD5	Status	Message
2024-01-01T00:00:00	SampleA	SampleA/SampleA_R1.fq.gz	d41d8c...	d41d8c...	PASS	
```

### 4.5 工具来源与编译

`seq_preprocessor` 和 `json_md5_verifier` 均为 Rust 编写的高性能命令行工具，源码托管在：

```
https://github.com/xsx123123/JZ_Tools.git
```

源码位于仓库的 `src/md5/` 目录下：

```text
JZ_Tools/
└── src/md5/
    ├── seq_preprocessor/
    │   ├── Cargo.toml
    │   └── src/main.rs
    └── json_md5_verifier/
        ├── Cargo.toml
        └── src/main.rs
```

### 4.6 Git 子模块集成

所有 Flow 流程将 JZ_Tools 作为 git submodule 引入，路径为 `src/JZ_Tools`：

```bash
cd /path/to/<X>Flow
git submodule add https://github.com/xsx123123/JZ_Tools.git src/JZ_Tools
git submodule update --init --recursive
```

`.gitmodules` 内容：

```ini
[submodule "src/JZ_Tools"]
    path = src/JZ_Tools
    url = https://github.com/xsx123123/JZ_Tools.git
```

### 4.7 编译方法

需要 Rust 工具链（rustc + cargo）：

```bash
# 安装 Rust（如未安装）
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source ~/.cargo/env

# 编译 seq_preprocessor
cd src/JZ_Tools/src/md5/seq_preprocessor
cargo build --release

# 编译 json_md5_verifier
cd ../json_md5_verifier
cargo build --release
```

编译后复制到标准路径：

```bash
mkdir -p src/src/md5/Compile_release/seq_preprocessor_x86_64/release/
mkdir -p src/src/md5/Compile_release/json_md5_verifier_x86_64/release/

cp src/JZ_Tools/src/md5/seq_preprocessor/target/release/seq_preprocessor \
   src/src/md5/Compile_release/seq_preprocessor_x86_64/release/

cp src/JZ_Tools/src/md5/json_md5_verifier/target/release/json_md5_verifier \
   src/src/md5/Compile_release/json_md5_verifier_x86_64/release/
```

---

## 5. 通用层规则

### 5.1 01.common.smk

导入 `rules/utils/` 中的公共函数，不定义业务 rule。所有后续 rule 文件依赖此文件提供的 helper。

### 5.2 03.short_read_qc.smk

- 对标准化后的 `00.raw_data/link_dir/{sample}/{sample}_R1/R2.fq.gz` 运行 FastQC。
- MultiQC 汇总所有样本的 FastQC 结果。
- 可选 FastQ Screen 污染检测（由 `parameter.fastq_screen.run` 控制）。
- **依赖 `01.qc/md5_check.tsv` 作为 input gate**。

### 5.3 交付规则（N-1.deliver.smk）

- 生成 `delivery_manifest.json`（交付文件清单）。
- 生成 `delivery_manifest.md5`（交付文件 MD5）。
- 可选 `delivery_details.log`（交付日志）。
- 由 `deliver: true` 开关控制。

### 5.4 报告规则（N.Report.smk）

- 最小实现：`project_summary.json`（项目元数据摘要）。
- 完整实现：Quarto / R Markdown HTML 报告。
- 由 `report: true` 开关控制。

---

## 6. 组学专用层

### 6.1 RNAFlow（Bulk RNA-seq）

```text
02.file_convert_md5 → 03.QC → 04.Contamination → 05.fastp clean
→ 06.STAR mapping → 07.RSEM count → 08.GATK variant
→ 09.StringTie → 10.DEG/Enrichment → 11.rMATS → 12.Fusion
→ 13.deliver → 14.Report
```

核心工具：STAR、RSEM、fastp、DESeq2、GATK、rMATS、Arriba。

样本表必填：`sample,sample_name,group,bam`（bam 可自动生成）。

### 6.2 SCRNAFlow（单细胞 RNA）

```text
02.file_convert_md5 → prepare_scrna_fastq_links
→ 03.QC → 04.Contamination → 05.merge_qc_report
→ 06.count (Cell Ranger / dnbc4tools) → 07.deliver → 08.Report
```

核心工具：Cell Ranger、dnbc4tools。

样本表必填：`sample,sample_name,group,platform`（platform: 10x/dnb）。

特殊点：
- 不做 read clean，原始 FASTQ 直接交给 count backend。
- 需要二次链接（bcl2fastq 命名 for Cell Ranger，cDNA/oligo for DNB）。
- 输出统一矩阵：`03.count/matrix/{sample}/filtered_feature_bc_matrix/`。

### 6.3 BTCRFlow（Bulk BCR/TCR）

```text
02.file_convert_md5 → 03.QC → 04.trimming (cutadapt -U 30)
→ 05.MiXCR analyze (BCR/TCR preset) → 06.VDJtools diversity
→ 07.deliver → 08.Report
```

核心工具：MiXCR v4、VDJtools 1.2.1、cutadapt。

样本表必填：`sample,sample_name,group,platform,path`（platform: bcr/tcr）。

特殊点：
- BCR 使用 preset `takara-human-rna-bcr-umi-smarter`。
- TCR 使用 preset `generic-amplicon`。
- Clontech SMARTer 协议需要 `-U 30` 硬切 5' 端。
- VDJtools 做 downsample + diversity 统计。

### 6.4 ATACFlow（染色质可及性，规划中）

```text
02.file_convert_md5 → 03.QC → 04.trimming → 05.BWA mapping
→ 06.dedup/filter → 07.MACS2 peak → 08.annotation → 09.deliver → 10.Report
```

### 6.5 WGSFlow（全基因组，规划中）

```text
02.file_convert_md5 → 03.QC → 04.trimming → 05.BWA mapping
→ 06.dedup/BQSR → 07.GATK variant → 08.CNV/SV → 09.deliver → 10.Report
```

---

## 7. 样本表约定（通用）

### 7.1 通用必填列

```csv
sample,sample_name,group
```

### 7.2 组学扩展列

| 流程 | 额外必填列 | 可选列 |
|---|---|---|
| RNAFlow | `bam`（可自动生成） | `condition, batch` |
| SCRNAFlow | `platform` | `fastq_id, chemistry, lanes, expect_cells` |
| BTCRFlow | `platform, path` | `index` |
| ATACFlow | `path` | `condition, batch` |

### 7.3 路径列多 lane 支持

路径列支持多个 lane 用逗号、分号或竖线分隔：

```csv
sample,path
SampleA,/data/lane1/SampleA,/data/lane2/SampleA
```

---

## 8. 模块开关与目标生成

### 8.1 通用模块开关

```yaml
qc_clean: true      # FastQC + MultiQC（所有流程）
deliver: true       # 交付 manifest（所有流程）
report: true        # 项目报告（所有流程）
```

### 8.2 组学专用开关

| 流程 | 专用开关 |
|---|---|
| RNAFlow | `trimming, mapping, count, deg, variant, rmats, fusion` |
| SCRNAFlow | `count` |
| BTCRFlow | `trimming, mixcr, vdj_analysis` |

### 8.3 模块依赖

```python
MODULE_DEPENDENCIES = {
    # RNAFlow
    "mapping": ["qc_clean", "trimming"],
    "count": ["mapping"],
    "deg": ["count"],
    # SCRNAFlow
    "count": ["qc_clean"],
    # BTCRFlow
    "mixcr": ["qc_clean", "trimming"],
    "vdj_analysis": ["mixcr"],
}
```

### 8.4 DataDeliver() 固定目标

所有流程的 `DataDeliver()` 固定加入：

```text
00.raw_data/link_dir/raw_data_md5.json
01.qc/md5_check.tsv
```

其余目标由各模块收集器按开关动态添加。

---

## 9. 容器化策略（通用）

### 9.1 分层原则

| 层级 | 工具 | 管理方式 |
|---|---|---|
| 通用层 | seq_preprocessor, json_md5_verifier | 仓库内编译二进制 |
| QC 层 | FastQC, MultiQC, fastq_screen, cutadapt | Conda env |
| 分析层（可 conda） | STAR, RSEM, MiXCR, BWA, GATK | Conda env |
| 分析层（不可 conda） | Cell Ranger, dnbc4tools | 容器 / bind mount |
| 参考基因组 | 所有 index | Bind mount，不打入镜像 |

### 9.2 推荐配置

```yaml
container:
  enabled: false
  engine: docker          # docker / apptainer / singularity
  images:
    cellranger: cellranger:10.1.0
    dnbc4tools: dnbc4tools:3.1
  bind:
    - /data:/data
    - /tmp:/tmp
  docker_args: ""
```

### 9.3 注意事项

- 参考基因组体积大、更新频繁，永远通过 bind mount 提供。
- Cell Ranger 有许可限制，不放入公共镜像仓库。
- 容器内路径必须和配置路径一致（通过 bind mount 保证）。
- `logs/`、`benchmarks/`、输出目录必须在容器内可写。

---

## 10. 验证命令（通用）

```bash
# Dry-run（验证 DAG）
snakemake -s snakefile -n --cores 1 \
  --config analysisyaml=/path/to/config.yaml

# 只运行到 MD5 校验
snakemake -s snakefile --cores 4 \
  --config analysisyaml=/path/to/config.yaml \
  01.qc/md5_check.tsv

# 完整运行
snakemake -s snakefile --cores 32 --use-conda \
  --config analysisyaml=/path/to/config.yaml
```

---

## 11. 维护原则

### 11.1 通用层维护

- `raw_data_md5.json` 是 MD5 校验层，不要让下游业务直接依赖它作为数据源。
- 改 rule `output` 时必须同步改 `rules/utils/datadeliver.py`。
- 改样本表字段时必须同步检查 `id_convert.py`、相关 scripts、`tools.py`。
- `seq_preprocessor` / `json_md5_verifier` 升级时，所有 Flow 流程同步更新 submodule。
- 文档以当前实现为准；通用规范不应覆盖各流程的实际边界。

### 11.2 组学层维护

- 新增组学流程时，复用 01-03 通用规则，只编写组学专用 rule 文件。
- 组学专用 manifest（如 `scrna_fastq_manifest.json`）不要和 `raw_data_md5.json` 耦合。
- 参考基因组 `index_keys` 按流程独立维护，不要交叉引入。

### 11.3 工具链维护

- JZ_Tools submodule 统一追踪 `main` 分支。
- 编译产物放入 `src/src/md5/Compile_release/` 并加入 `.gitignore`（或 LFS）。
- Rust 工具升级后需重新编译并验证所有 Flow 流程的 02 步 dry-run。

---

## 12. 各流程 02 步后续对比

| 流程 | 02 步输出 | 后续业务 manifest | 核心分析入口 |
|---|---|---|---|
| RNAFlow | `link_dir/{sample}/{sample}_R1/R2.fq.gz` | 无（直接用 link_dir） | STAR mapping |
| SCRNAFlow | `link_dir/` + `scrna_fastq_manifest.json` | `scrna_fastq_manifest.json` | Cell Ranger / dnbc4tools |
| BTCRFlow | `link_dir/{sample}/{sample}_R1/R2.fq.gz` | 无（直接用 link_dir） | MiXCR analyze |
| ATACFlow | `link_dir/{sample}/{sample}_R1/R2.fq.gz` | 无（直接用 link_dir） | BWA mapping |

三者共享同一套 `seq_preprocessor` + `json_md5_verifier` 二进制和 `raw_data_md5.json` 格式，仅下游业务 manifest 和分析入口不同。
