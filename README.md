# circRNA_hepatic

A bioinformatics pipeline for the identification, quantification, and differential expression analysis of circular RNAs (circRNAs) in hepatic tissues.

## Overview

This project provides a comprehensive workflow for analyzing circRNA expression profiles from RNA-seq data. It integrates several specialized tools to handle the unique challenges of circRNA discovery and quantification.

## Pipeline Stages

The analysis is divided into several main stages, each located in the `scripts/` directory:

1.  **Preprocessing**: Trimming and quality control of raw reads using Trimmomatic.
2.  **Discovery**: De novo identification of circRNAs using CIRI2.
3.  **Quantification**: Expression level estimation using Salmon.
4.  **Differential Expression**: Statistical analysis of circRNAs, genes, and miRNAs using DESeq2.
5.  **Downstream Analysis**: Ranking and visualization of top candidates.

## Directory Structure

*   `scripts/trimmomatic/`: Read preprocessing scripts.
*   `scripts/ciri2/`: CircRNA identification scripts.
*   `scripts/salmon_circbase/`: Quantification scripts.
*   `scripts/deSeq2/`: Differential expression analysis.
*   `scripts/functions/`: Common R utility functions.
*   `scripts/ranking_circrna/`: Prioritization and plotting scripts.

## Usage

Most analysis stages can be run using the "run_all" scripts provided in each directory. For example, to run the DESeq2 analysis for circRNAs:

```bash
Rscript scripts/deSeq2/run_all_compounds.R circ
```

*Note: Many scripts contain hardcoded paths specific to the original computing cluster environment. Ensure paths are updated before running in a new environment.*

## Dependencies

*   Bash / R
*   CIRI2
*   Salmon
*   Trimmomatic
*   R Packages: `DESeq2`, `dplyr`, `biomaRt`, etc.
