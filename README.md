# circRNA_hepatic

A bioinformatics pipeline for the identification, quantification, and differential expression analysis of circular RNAs (circRNAs) in hepatic tissues.

## Overview

This project provides a comprehensive workflow for analyzing circRNA expression profiles from RNA-seq data. It integrates several specialized tools to handle the unique challenges of circRNA discovery and quantification, specifically tailored for studies involving drug treatments on hepatic models.

## Project Structure

This repository follows standard bioinformatics organizational practices. **Note: Due to privacy and security guidelines, no experimental data is hosted in this repository.**

```text
circRNA_hepatic/
├── data/                  # Ignored by git: Directory for datasets
│   ├── raw/               # Raw FASTQ files (do not modify)
│   └── processed/         # Intermediate and processed files
├── results/               # Ignored by git: Final analysis results and plots
├── envs/                  # Conda environment specifications
├── docs/                  # Project documentation, protocols, and metadata templates
├── logs/                  # Ignored by git: Log files from pipeline executions
└── scripts/               # Core pipeline scripts
    ├── trimmomatic/       # Read preprocessing scripts
    ├── ciri2/             # CircRNA identification scripts
    ├── salmon_circbase/   # Quantification scripts
    ├── deSeq2/            # Differential expression analysis
    ├── pre-DESeq2/        # Merging Salmon outputs for DESeq2
    ├── ranking_circrna/   # Prioritization and plotting scripts
    └── functions/         # Common R utility functions
```

## Pipeline Stages

The analysis is divided into several main stages. Execution relies on master scripts found in each stage's directory:

1.  **Preprocessing**: Trimming and quality control of raw reads using Trimmomatic.
    ```bash
    ./scripts/trimmomatic/trimming_v7_JOA.sh
    ```
2.  **Discovery**: De novo identification of circRNAs using CIRI2.
    ```bash
    ./scripts/ciri2/run_all_compounds.sh
    ```
3.  **Quantification**: Expression level estimation using Salmon.
    ```bash
    ./scripts/salmon_circbase/run_all_compounds_salmon.sh
    ```
4.  **Differential Expression**: Statistical analysis of circRNAs, genes, and miRNAs using DESeq2.
    ```bash
    Rscript scripts/deSeq2/run_all_compounds.R circ   # For circRNAs
    Rscript scripts/deSeq2/run_all_compounds.R gene   # For genes
    Rscript scripts/deSeq2/run_all_compounds.R mirna  # For miRNAs
    ```
5.  **Downstream Analysis**: Ranking and visualization of top candidates in `scripts/ranking_circrna/`.

## Dependencies and Installation

A generic Conda environment configuration is provided in `envs/environment.yml`. To recreate the analysis environment (excluding specific manual installations like CIRI2):

```bash
conda env create -f envs/environment.yml
conda activate circrna-pipeline
```

*Note: CIRI2 requires manual installation. Ensure all binaries (CIRI2, BWA, Salmon, Trimmomatic) are available in your system `$PATH` or appropriately referenced in the scripts.*

## Configuration

A template for managing pipeline parameters and paths is provided in `config.yaml.template`. To configure your pipeline:

1.  Copy `config.yaml.template` to `config.yaml`:
    ```bash
    cp config.yaml.template config.yaml
    ```
2.  Edit `config.yaml` to specify the correct paths for your data, results, logs, and reference files, as well as any tool-specific parameters.

## Configuration and Adaptability

**Important Note for Portability:** Many scripts in this pipeline may contain hardcoded paths specific to the original computing cluster environment (e.g., `/share/analysis/hecatos/...`). The `config.yaml` file is designed to help centralize these settings.

Before executing the pipeline in a new environment, you must:
1.  **Set up `config.yaml`** as described above, customizing all relevant paths and parameters.
2.  Review the `.sh` and `.R` master scripts to ensure they reference variables from `config.yaml` (this may require manual updates to the scripts).
3.  Ensure experimental metadata is correctly defined (often located in `scripts/deSeq2/metadata_deseq2.R`, though ideally, this would also be configurable).

## License

This project is intended for research purposes. (Include appropriate License details here).