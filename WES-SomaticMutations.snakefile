### gatk/4.6.1.0
### picard/2.26.10
### fastqc/0.11.9
### multiqc/1.14
### sambamba/0.7.1
### fastp/0.12.4
### bwa/0.7.18
### samtools/1.21
### java/1.6.0_27

configfile: "config.yaml"

import os

wkdir = config["wkdir"]
REFERENCE = config.get("reference", "/home/reference/hg38/Homo_sapiens_assembly38.fasta")

# Define working subdirectories
wkFastqDir = os.path.join(wkdir, "FastqDir")
wkFastQCDir = os.path.join(wkdir, "FastQCDir")
wkFastpDir = os.path.join(wkdir, "FastpDir")
wkVCFDir = os.path.join(wkdir, "VCFDir")
wkBamDir = os.path.join(wkdir, "BamDir")
wkMultiQCDir = os.path.join(wkdir, "MultiQCDir")
wkAdapterF = os.path.join(wkdir, "AdapterF")
wkPreprocessfastqDir = os.path.join(wkdir, "PreprocessfastqDir")
wkGATKBundleDir = os.path.join(wkdir, "GATKBundleDir")
wkGVCFDir = os.path.join(wkdir, "GVCFDir")
wkTempDir = os.path.join(wkdir, "TempDir")

chrs = [str(i) for i in range(1, 23)] + ["X", "Y"]

rule all:
    input:
        expand(os.path.join(wkVCFDir, "{sample}.Mutect2Soma.snpEff"), sample=config["samples"])

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule FastQC:
    """
    Run FastQC for quality control of raw FASTQ files
    """
    input:
        fastq1=os.path.join(wkFastqDir, "{sample}.R1.clean.fastq.gz"),
        fastq2=os.path.join(wkFastqDir, "{sample}.R2.clean.fastq.gz")
    output:
        FastQCzip1=os.path.join(wkFastQCDir, "{sample}.1_fastqc.zip"),
        FastQCzip2=os.path.join(wkFastQCDir, "{sample}.2_fastqc.zip"),
        FastQChtml1=os.path.join(wkFastQCDir, "{sample}.1_fastqc.html"),
        FastQChtml2=os.path.join(wkFastQCDir, "{sample}.2_fastqc.html")
    params:
        OutputDir=wkFastQCDir
    threads: 8
    resources:
        mem_mb=20000
    shell:
        """
        mkdir -p {params.OutputDir}
        fastqc -o {params.OutputDir} -t {threads} {input.fastq1} {input.fastq2}
        """

rule fastp:
    """
    Perform adapter trimming and quality filtering with fastp
    """
    wildcard_constraints:
        sample = "[^\.]+"
    params:
        OutputDir=wkFastpDir,
        sample="{sample}"
    input:
        fastq1=os.path.join(wkFastqDir, "{sample}.R1.clean.fastq.gz"),
        fastq2=os.path.join(wkFastqDir, "{sample}.R2.clean.fastq.gz")
    output:
        fastqpreprocess1=os.path.join(wkFastpDir, "{sample}.fastp_1.fastq.gz"),
        fastqpreprocess2=os.path.join(wkFastpDir, "{sample}.fastp_2.fastq.gz"),
        fastpjson=os.path.join(wkFastpDir, "{sample}.fastp.json"),
        fastphtml=os.path.join(wkFastpDir, "{sample}.fastp.html")
    threads: 8
    resources:
        mem_mb=20000
    shell:
        """
        mkdir -p {params.OutputDir}
        fastp --detect_adapter_for_pe -w {threads} \
            -i {input.fastq1} -o {output.fastqpreprocess1} \
            -I {input.fastq2} -O {output.fastqpreprocess2} \
            -j {output.fastpjson} -h {output.fastphtml} \
            -R {params.sample}
        """

rule BWAmemAnd2BAM:
    """
    Map reads with bwa mem and convert to BAM
    """
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        fq1=os.path.join(wkFastpDir, "{sample}.fastp_1.fastq.gz"),
        fq2=os.path.join(wkFastpDir, "{sample}.fastp_2.fastq.gz")
    output:
        bam=os.path.join(wkBamDir, "{sample}.bam")
    threads: 8
    resources:
        mem_mb=20000
    shell:
        """
        bwa mem -t {threads} {REFERENCE} {input.fq1} {input.fq2} | samtools view -Sb - > {output.bam}
        """

rule AddOrReplaceReadGroups:
    """
    Add or replace read groups using GATK
    """
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        bam=os.path.join(wkBamDir, "{sample}.bam")
    output:
        addRGbam=os.path.join(wkBamDir, "{sample}.RG.bam")
    params:
        ID="{sample}",
        LB="{sample}",
        PL="illumina",
        PU="{sample}",
        SM="{sample}"
    threads: 8
    resources:
        mem_mb=20000
    shell:
        """
        gatk --java-options '-Xms15G -Xmx15G -XX:ParallelGCThreads=8' AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.addRGbam} \
            -ID {params.ID} -LB {params.LB} -PL {params.PL} -PU {params.PU} -SM {params.SM}
        """

rule sambambaSortMarkdupIndex:
    """
    Sort, mark duplicates and index BAM with sambamba
    """
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        addRGbam=os.path.join(wkBamDir, "{sample}.RG.bam")
    output:
        sortbam=os.path.join(wkBamDir, "{sample}.sorted.bam"),
        markdupbam=os.path.join(wkBamDir, "{sample}.sorted.markdup.bam"),
        bamindex=os.path.join(wkBamDir, "{sample}.sorted.markdup.bam.bai")
    params:
        tmpDir=os.path.join(wkBamDir, "{sample}_temp")
    threads: 8
    shell:
        """
        mkdir -p {params.tmpDir}
        sambamba sort -t {threads} -m 20GB -p --tmpdir {params.tmpDir} -o {output.sortbam} {input.addRGbam}
        sambamba markdup -t {threads} -p --tmpdir {params.tmpDir} {output.sortbam} {output.markdupbam}
        sambamba index -t {threads} -p {output.markdupbam} {output.bamindex}
        """

rule BQSR:
    """
    Base quality score recalibration (BQSR) with GATK BaseRecalibrator
    """
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam=os.path.join(wkBamDir, "{sample}.sorted.markdup.bam"),
        dbSNP_vcf="/home/resource/hg38/Homo_sapiens_assembly38.dbsnp138.vcf",
        known_indels_sites_vcf="/home/reference/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz",
        Mills_and_1000G_gold_standard_indels_vcf="/home/resource/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    output:
        recalibration_report_filename=os.path.join(wkBamDir,"{sample}.sorted.markdup.BQSRrecalibration.table")
    params:
        tmpDir=os.path.join(wkBamDir, "{sample}_temp")
    threads: 8
    shell:
        """
        mkdir -p {params.tmpDir}
        /home/bin/java -Xms15G -Xmx15G -XX:ParallelGCThreads=8 \
            -jar /home/software/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar BaseRecalibrator \
            -I {input.markdupbam} \
            -R {REFERENCE} \
            --use-original-qualities \
            -O {output.recalibration_report_filename} \
            --known-sites {input.dbSNP_vcf} \
            --known-sites {input.known_indels_sites_vcf} \
            --known-sites {input.Mills_and_1000G_gold_standard_indels_vcf} \
            --tmp-dir {params.tmpDir}
        """

rule ApplyBQSR:
    """
    Apply BQSR to BAM file
    """
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam=os.path.join(wkBamDir, "{sample}.sorted.markdup.bam"),
        recalibration_report_filename=os.path.join(wkBamDir,"{sample}.sorted.markdup.BQSRrecalibration.table")
    output:
        BQSRbam=os.path.join(wkBamDir,"{sample}.sorted.markdup.BQSR.bam"),
        BQSRbambai=os.path.join(wkBamDir,"{sample}.sorted.markdup.BQSR.bai")
    params:
        tmpDir=os.path.join(wkBamDir, "{sample}_temp")
    threads: 8
    shell:
        """
        /home/bin/java -Xms15G -Xmx15G -XX:ParallelGCThreads=8 \
            -jar /home/software/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar ApplyBQSR \
            -I {input.markdupbam} \
            -R {REFERENCE} \
            --use-original-qualities \
            --bqsr-recal-file {input.recalibration_report_filename} \
            -O {output.BQSRbam} \
            --tmp-dir {params.tmpDir}
        """

rule GetSummary1:
    """
    Generate flagstat summary with samtools
    """
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam=os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam")
    output:
        flagstat=os.path.join(wkMultiQCDir,"{sample}.flagstat")
    params:
        tmpDir=os.path.join(wkMultiQCDir, "{sample}_temp")
    threads: 8
    shell:
        """
        mkdir -p {params.tmpDir}
        /home/samtools flagstat {input.markdupbam} > {output.flagstat}
        """

rule GetSummary2:
    """Run samtools stats to summarize BAM file statistics."""
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam")
    output:
        stats = os.path.join(wkMultiQCDir, "{sample}.stats")
    threads: 8
    shell:
        """
        /home/samtools stats {input.markdupbam} > {output.stats}
        """


rule GetSummary3:
    """Run pandepth to generate per-chromosome depth statistics."""
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam")
        bed = "/home/resource/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"
    output:
        chrstats = os.path.join(wkMultiQCDir, "{sample}.chrstats")
    threads: 8
    shell:
        """
        /home/pandepth -b {input.bed} -i {input.markdupbam} -o {output.chrstats}
        """

rule Mutect2RMGerm:
    """Run Mutect2 for somatic variant calling using only tumor sample and remove germline."""
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam")
    output:
        mutect2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2RMGerm.vcf.gz")
    threads: 8
    shell:
        """
        /home/bin/java -Xms15G -Xmx15G -XX:ParallelGCThreads=8 \
            -jar /home/software/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar Mutect2 \
            -R /home/resource/reference/hg38/Homo_sapiens_assembly38.fasta \
            -I {input.markdupbam} \
            --germline-resource /home/resource/hg38/af-only-gnomad.hg38.vcf.gz \
            --panel-of-normals /home/resource/hg38/1000g_pon.hg38.vcf.gz \
            -O {output.mutect2vcf}
        """


rule Mutect2Soma:
    """Run Mutect2 with matched normal and tumor BAM files."""
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam")
    output:
        mutect2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.vcf.gz")
    threads: 8
    shell:
        """
        /home/bin/java -Xms15G -Xmx15G -XX:ParallelGCThreads=8 \
            -jar /home/software/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar Mutect2 \
            -R /home/resource/reference/hg38/Homo_sapiens_assembly38.fasta \
            -I {input.markdupbam} \
            -I Normal.sorted.markdup.BQSR.bam \
            --germline-resource /home/resource/hg38/af-only-gnomad.hg38.vcf.gz \
            --panel-of-normals /home/resource/hg38/1000g_pon.hg38.vcf.gz \
            -O {output.mutect2vcf}
        """


rule Filtered2Soma:
    """Filter Mutect2 somatic calls using GATK FilterMutectCalls."""
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        mutect2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.vcf.gz")
    output:
        filtered2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.filtered.vcf")
    threads: 8
    shell:
        """
        /home/bin/java -Xms15G -Xmx15G -XX:ParallelGCThreads=8 \
            -jar /home/software/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar FilterMutectCalls \
            -R /home/reference/hg38/Homo_sapiens_assembly38.fasta \
            --min-median-base-quality 20 \
            --min-median-mapping-quality 30 \
            -V {input.mutect2vcf} \
            -O {output.filtered2vcf}
        """


rule snpEff2Soma:
    """Annotate filtered VCF with SnpEff."""
    input:
        filtered2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.filtered.vcf")
    output:
        snpEff2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.snpEff")
    threads: 8
    shell:
        """
        genome="GRCh38.86"
        /home/bin/java -Xmx4g -jar /home/software/snpEff/snpEff.jar \
            -v $genome {input.filtered2vcf} > {output.snpEff2vcf}
        """
