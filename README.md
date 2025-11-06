# RNA-seq Differential Expression Analysis Pipeline

A complete, automated NextFlow pipeline for RNA-seq data analysis including quality control, alignment, quantification, and differential expression analysis.

## Overview

This pipeline performs end-to-end RNA-seq analysis from raw FASTQ files to differential expression results. It processes paired-end sequencing data through quality control, genome alignment, gene quantification, and statistical analysis to identify differentially expressed genes between experimental conditions.

**Key Features:**
- Automated workflow management with NextFlow
- Containerized execution for reproducibility
- Multi-sample parallel processing
- Comprehensive quality control reporting
- Statistical analysis with DESeq2
- Pathway enrichment analysis

## Pipeline Workflow

```
FASTQ Files (Paired-end)
    ├── FastQC → Quality Control Metrics
    ├── STAR Alignment → BAM Files
    │       └── MultiQC → Aggregated QC Report
    ├── VERSE Quantification → Gene Counts
    └── Count Matrix Generation
            └── DESeq2 Analysis → Differential Expression Results
                    ├── PCA & Sample Distance Plots
                    ├── Volcano Plots
                    └── GSEA/Pathway Enrichment
```

## Requirements

### Software Dependencies
- NextFlow (≥ 21.0)
- Singularity/Docker
- R (≥ 4.0) with Bioconductor packages

### Containerized Tools
All bioinformatics tools are provided as Docker containers:
- FastQC (`ghcr.io/bf528/fastqc:latest`)
- MultiQC (`ghcr.io/bf528/multiqc:latest`)
- STAR (`ghcr.io/bf528/star:latest`)
- VERSE (`ghcr.io/bf528/verse:latest`)
- Pandas/Python (`ghcr.io/bf528/pandas:latest`)

### Input Data
- Paired-end FASTQ files (`.fastq.gz`)
- Reference genome (FASTA format)
- Gene annotation (GTF format)

## Installation

1. **Clone the repository:**
```bash
git clone <repository-url>
cd rnaseq-pipeline
```

2. **Configure paths:**
Edit `nextflow.config` to specify:
```groovy
params {
    reads = "/path/to/fastq/*_R{1,2}.fastq.gz"
    genome = "/path/to/reference/genome.fa"
    gtf = "/path/to/annotation/genes.gtf"
    results = "./results"
}
```

3. **Ensure Singularity is available** (for container execution)

## Usage

### Basic Execution

Run the complete pipeline:
```bash
nextflow run main.nf -profile singularity,cluster
```

### Local Testing
For testing with subset data:
```bash
nextflow run main.nf -profile singularity,local
```

### Resume Failed Runs
```bash
nextflow run main.nf -profile singularity,cluster -resume
```

## Pipeline Components

### 1. Quality Control (FastQC)
Analyzes raw sequencing data quality including:
- Per-base sequence quality
- Sequence duplication levels
- Adapter content
- GC content distribution

**Output:** `*.html` and `*.zip` reports per sample

### 2. Genome Index Generation (STAR)
Creates an indexed reference genome for efficient alignment.

**Command:**
```bash
STAR --runMode genomeGenerate \
     --genomeDir <index_dir> \
     --genomeFastaFiles <genome.fa> \
     --sjdbGTFfile <genes.gtf> \
     --runThreadN 16
```

### 3. Read Alignment (STAR)
Aligns paired-end reads to the reference genome.

**Features:**
- Splice-aware alignment
- Multi-threaded processing
- BAM format output
- Alignment statistics logging

**Expected alignment rate:** >70% for well-prepared samples

### 4. Aggregated QC Report (MultiQC)
Combines FastQC and STAR logs into a single interactive HTML report for easy comparison across samples.

### 5. Gene Quantification (VERSE)
Counts reads mapping to gene features from GTF annotation.

**Output:** Gene-level count tables per sample

### 6. Count Matrix Assembly
Python script using Pandas to concatenate individual sample counts into a single matrix:
- Rows: Genes (Ensembl IDs)
- Columns: Samples
- Values: Raw read counts

### 7. Differential Expression Analysis (DESeq2)
Statistical analysis performed in R using DESeq2:

```r
# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = metadata,
                               design = ~ condition)

# Run analysis
dds <- DESeq(dds)
results <- results(dds, alpha = 0.05)
```

**Outputs:**
- Differential expression results table
- Normalized counts
- PCA plots
- Sample distance heatmaps
- Volcano plots

### 8. Pathway Enrichment Analysis
- **DAVID/ENRICHR:** GO terms and pathway enrichment
- **GSEA (fgsea):** Gene Set Enrichment Analysis using ranked gene lists

## Project Structure

```
rnaseq-pipeline/
├── main.nf                 # Main NextFlow workflow
├── nextflow.config         # Configuration file
├── modules/                # Process definitions
│   ├── fastqc.nf
│   ├── star_index.nf
│   ├── star_align.nf
│   ├── verse.nf
│   ├── multiqc.nf
│   └── concat_counts.nf
├── bin/                    # Custom scripts
│   ├── parse_gtf.py
│   └── concat_verse_counts.py
└── analysis/               # R analysis scripts
    └── deseq2_analysis.Rmd
```

## Key Scripts

### parse_gtf.py
Extracts Ensembl ID to gene symbol mappings from GTF file:
```python
# Usage
python parse_gtf.py --gtf genes.gtf --output gene_names.txt
```

### concat_verse_counts.py
Merges individual sample count files into a count matrix:
```python
# Usage
python concat_verse_counts.py --input counts/*.txt --output count_matrix.csv
```

## Outputs

### Directory Structure
```
results/
├── fastqc/                 # Individual FastQC reports
├── multiqc_report.html     # Aggregated QC report
├── bam/                    # Alignment files
├── verse/                  # Gene count files
├── count_matrix.csv        # Combined counts matrix
└── deseq2/                 # Differential expression results
    ├── results_table.csv
    ├── pca_plot.pdf
    ├── volcano_plot.pdf
    └── enrichment_results.csv
```

### Key Results Files

**count_matrix.csv** - Raw gene counts matrix  
**results_table.csv** - DESeq2 output with log2FC and adjusted p-values  
**significant_genes.csv** - Filtered list of DE genes (padj < 0.05)  
**enrichment_results.csv** - Pathway analysis results  

## Analysis Parameters

### Filtering Strategy
Genes are filtered before DESeq2 analysis:
- Minimum count threshold: 10 reads across all samples
- Rationale: Removes low-abundance genes unlikely to be biologically meaningful

### Statistical Thresholds
- **Adjusted p-value (padj):** < 0.05 (Benjamini-Hochberg correction)
- **Log2 Fold Change:** Optional threshold (|log2FC| > 1)

### Normalization
DESeq2 uses median-of-ratios normalization. For visualization:
- **rlog:** For samples with n < 30
- **vst:** For samples with n ≥ 30

## Computational Resources

### Resource Requirements
Configured via process labels:

**process_low:** 1 CPU, 4 GB RAM (FastQC, MultiQC)  
**process_medium:** 8 CPUs, 32 GB RAM (VERSE)  
**process_high:** 16 CPUs, 128 GB RAM (STAR alignment/indexing)  

### Estimated Runtime
For 6 samples (~50M reads each):
- FastQC: ~10 min per sample
- STAR indexing: ~30 min (one-time)
- STAR alignment: ~20 min per sample
- VERSE: ~5 min per sample
- **Total:** ~3-4 hours on cluster

## Quality Control Metrics

### Acceptable Ranges
- **Alignment rate:** >70% for human/mouse genomes
- **Duplication rate:** <50% (library complexity)
- **Per-base quality:** Phred score >28
- **GC content:** ~50% for human (expect normal distribution)

### Red Flags
- Low alignment rates (<60%) → contamination or wrong reference
- High duplication (>70%) → PCR over-amplification
- Adapter content detected → incomplete trimming

## Reproducibility

### Version Control
All software versions are locked in containers for reproducibility.

### Session Information
R analysis includes `sessionInfo()` output documenting all package versions.

### Container Specifications
Containers are built from specifications at: `github.com/BF528/pipeline_containers`

## Troubleshooting

### Common Issues

**Issue:** STAR alignment fails with memory error  
**Solution:** Increase memory allocation in `process_high` label or reduce `--genomeSAindexNbases`

**Issue:** DESeq2 reports "every gene contains at least one zero"  
**Solution:** Filtering threshold too stringent; relax minimum count requirement

**Issue:** MultiQC report is empty  
**Solution:** Check that FastQC `.zip` files (not `.html`) are being passed to MultiQC

## Example Results

### Sample Dataset
Analysis of pancreatic beta cell-derived iPSCs (n=6):
- 3 control samples
- 3 treated samples
- ~23,000 genes quantified

### Key Findings
- **1,227 significantly differentially expressed genes** (padj < 0.05)
- Enriched pathways: PI3K-Akt signaling, p53 pathway, neuronal system
- Strong sample clustering by condition in PCA (PC1: 68% variance)

## Citation

If using this pipeline, please cite:

- **STAR aligner:** Dobin et al. (2013) STAR: ultrafast universal RNA-seq aligner. Bioinformatics.
- **DESeq2:** Love et al. (2014) Moderated estimation of fold change and dispersion for RNA-seq data. Genome Biology.
- **FastQC:** Andrews, S. (2010) FastQC: a quality control tool for high throughput sequence data.
- **MultiQC:** Ewels et al. (2016) MultiQC: summarize analysis results for multiple tools and samples. Bioinformatics.

## License

This pipeline is provided for educational and research purposes.

## Contact

For questions or issues, please open a GitHub issue or contact:
- **Author:** Shreya Das
- **Email:** dshreya@bu.edu
- **Institution:** Boston University

## Acknowledgments

Pipeline developed as part of BF528: Applications in Translational Bioinformatics at Boston University.
