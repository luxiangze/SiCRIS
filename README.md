# SiCRIS: Silkworm integrated CRISPR analysis

SiCRIS is a specialized toolkit designed for the analysis of CRISPR/Cas9 screening libraries in Bombyx mori (silkworm). It provides a streamlined pipeline for processing sequencing reads, quantifying gRNA abundances, performing differential analysis, and identifying candidate target genes with functional relevance.

## Features
- Optimized for Bombyx mori genome and gRNA libraries
- High-throughput read count extraction and normalization
- Statistical analysis for screening experiments (e.g., positive/negative selection)
- Visualization tools for gRNA and gene-level enrichment
- Customizable parameters and modular workflow

## Use Cases
- Genome-wide CRISPR knockout screening in silkworm
- Functional genomics studies in lepidopteran insects
- Identification of essential genes under specific conditions

## Installation

### 环境配置

项目依赖可通过conda/mamba环境管理器安装。我们提供了完整的环境配置文件`environment.yml`。

```bash
# 使用conda创建环境
conda env create -f environment.yml

# 或使用mamba创建环境（推荐，速度更快）
mamba env create -f environment.yml

# 激活环境
conda activate sicris
```

### 断点续传功能

所有脚本都支持断点续传功能，当脚本执行中断后，可以从上次停止的地方继续执行，无需重新开始，节省时间。

## Usage

```bash
mamba activate sicris
# 仅处理样本并准备mageck输入文件
/home/gyk/project/SiCRIS/scripts/batch_process_for_mageck.sh /home/gyk/project/SiCRIS/data/samplesheet.csv 100

# 处理样本，合并gene和promoter数据，并自动运行mageck mle分析
/home/gyk/project/SiCRIS/scripts/batch_process_for_mageck.sh /home/gyk/project/SiCRIS/data/samplesheet.csv 100 --merge --run-mageck --mode=mle
```