# GEMINI.md - circRNA_hepatic Project Overview

## Project Overview
The `circRNA_hepatic` project is a bioinformatics pipeline designed for the identification, quantification, and differential expression analysis of circular RNAs (circRNAs) in hepatic samples. This project is part of a larger study (likely the "hecatos" project) involving drug treatments (compounds) on hepatic models.

The pipeline integrates several standard bioinformatics tools and custom R/Bash scripts to transform raw sequencing data into biological insights.

### Main Technologies
- **Bash**: Used for pipeline orchestration and automation.
- **R**: Primary language for statistical analysis and data visualization.
- **CIRI2**: De novo identification of circRNAs from RNA-seq data.
- **Salmon**: Quantification of transcript/circRNA expression.
- **Trimmomatic**: Pre-processing of raw FASTQ reads.
- **DESeq2**: Statistical framework for differential expression analysis.

---

## Project Structure

### `scripts/`
The `scripts/` directory is the core of the project, organized by the stage of analysis:

- **`trimmomatic/`**: Contains scripts for quality control and adapter trimming of raw sequencing reads.
- **`ciri2/`**: Includes scripts to run the CIRI2 pipeline for circRNA discovery.
- **`salmon_circbase/`**: Contains scripts for quantifying circRNAs (and potentially linear transcripts) using Salmon.
- **`deSeq2/`**: Scripts for performing differential expression analysis. It handles different types of features: circRNAs, genes, and miRNAs.
- **`functions/`**: Reusable R functions (e.g., `functions_JOA.R`) used across multiple scripts for data cleaning, ID mapping, and package management.
- **`pre-DESeq2/`**: Scripts for preparing and merging quantification output (e.g., from Salmon) into a format suitable for DESeq2.
- **`ranking_circrna/`**: Scripts for prioritizing circRNA candidates and generating visualizations like barplots.

---

## Building and Running

The pipeline is typically executed in stages using the "run_all" or specific master scripts. Note that many scripts contain hardcoded paths (e.g., `/share/analysis/hecatos/...`) specific to the original execution environment.

### 1. Trimming
```bash
./scripts/trimmomatic/trimming_v7_JOA.sh
```

### 2. CircRNA Identification (CIRI2)
```bash
./scripts/ciri2/run_all_compounds.sh
```

### 3. Quantification (Salmon)
```bash
./scripts/salmon_circbase/run_all_compounds_salmon.sh
```

### 4. Differential Expression (DESeq2)
The master R script accepts arguments to specify the analysis type:
```bash
Rscript scripts/deSeq2/run_all_compounds.R circ   # For circRNAs
Rscript scripts/deSeq2/run_all_compounds.R gene   # For genes
Rscript scripts/deSeq2/run_all_compounds.R mirna  # For miRNAs
```

---

## Development Conventions

- **Modularity**: Bash scripts frequently use `source` or `.` to include specific logic from other scripts in the same folder.
- **Shared Functions**: R scripts rely on `scripts/functions/functions_JOA.R`. Always check this file for existing helpers (e.g., `forceLibrary`, `naToZero`, `transcrToGene`) before implementing new utility functions.
- **Metadata Driven**: DESeq2 scripts use metadata files (often defined in `scripts/deSeq2/metadata_deseq2.R`) to define experimental groups and conditions.
- **Path Management**: Absolute paths are common in scripts. When porting to a new environment, these will likely need to be updated or parameterized.
