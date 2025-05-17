# Some Common Pipelines for Sequencing Data Processing (hg38)

## Somatic Mutation Calling in WGS data

This Snakemake workflow automates whole-genome sequencing (WGS) data processing for human samples aligned to the hg38 reference genome. It covers the following major steps:

- Quality control with FastQC and fastp
- Read mapping with BWA-MEM
- BAM processing (adding read groups, sorting, marking duplicates)
- Base Quality Score Recalibration (BQSR) with GATK
- Somatic variant calling using GATK Mutect2
- Summary statistics with samtools and MultiQC
- Functional annotation with snpEff

## Requirements

- Snakemake (>=5.0)
- Java 1.6 or later
- GATK 4.6.1.0
- Picard 2.26.10
- FastQC 0.11.9
- MultiQC 1.14
- Sambamba 0.7.1
- Fastp 0.12.4
- BWA 0.7.18
- Samtools 1.21
- snpEff (configured for hg38)

## Configuration

Edit `config.yaml` to specify:

- Working directory (`wkdir`)
- Reference genome path (`reference`)
- Samples list and FASTQ file paths

Example `config.yaml` is provided.

## Usage

1. Prepare `WGS-SomaticMutations.yaml` with your samples and paths.
2. Run the pipeline with:

```bash
snakemake --snakefile WGS-SomaticMutations.snakefile --configile WGS-SomaticMutations.yaml --cores 8
```

Adjust cores and resources as needed.

## Directory Structure

- FastqDir: raw FASTQ files
- FastQCDir: FastQC outputs
- FastpDir: trimmed FASTQ files by fastp
- BamDir: BAM files (mapped and processed)
- VCFDir: Variant calling outputs
- MultiQCDir: MultiQC reports and summaries
- AdapterF: Adapter sequences (optional)
- PreprocessfastqDir, GATKBundleDir, GVCFDir, TempDir: intermediate files

## Contact

For issues, please open an issue on the GitHub repository or contact the maintainer.

---

**Author:** Mengyue Zheng  
**Date:** 2025-05-17
