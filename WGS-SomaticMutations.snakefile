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
import glob  

wkdir = config["wkdir"]
REFERNCE = "/home/reference/hg38/Homo_sapiens_assembly38.fasta"
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
            expand(os.path.join(wkVCFDir, "{sample}.Mutect2Soma.snpEff"),sample=config["samples"]),

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule FastQC:
        ## Use FastQC for Quality Control of FASTQ Files 
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
        threads:8
        resources:
                mem_mb = 20000
        shell:
                """
                mkdir -p {params.OutputDir} # -p create the directory and, if required, all parent directories.                
                fastqc -o {params.OutputDir} -t 2 {input.fastq1} {input.fastq2} # -o --outdir # -t --threads
                """


rule fastp:
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
        threads:8
        resources:
                mem_mb = 20000
        shell:
                "mkdir -p {params.OutputDir} ;"
                "fastp" 
                " --detect_adapter_for_pe" # --detect_adapter_for_pe \
                " -w 8" # -w, --thread
                " -i {input.fastq1}" # -i, --in1
                " -o {output.fastqpreprocess1}" # -o, --out1
                " -I {input.fastq2}" # -I, --in2
                " -O {output.fastqpreprocess2}" # -O, --out2
                " -j {output.fastpjson}" # -j, --json the json format report file name
                " -h {output.fastphtml}" # -h, --html the html format report file name
                " -R {params.sample}" # -R, --report_title should be quoted with ' or ", default is "fastp report"


rule BWAmemAnd2BAM:
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        fq1 = os.path.join(wkFastpDir, "{sample}.fastp_1.fastq.gz"),
        fq2 = os.path.join(wkFastpDir, "{sample}.fastp_2.fastq.gz")
    output:
        bam = os.path.join(wkBamDir, "{sample}.bam")
    threads: 8
    resources:
        mem_mb = 20000
    shell:
        """
        bwa mem -t {threads} {REFERNCE} {input.fq1} {input.fq2} | samtools view -Sb - > {output.bam}
        """

rule AddOrReplaceReadGroups:
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
    threads:8
    resources:
        mem_mb = 20000
    shell:
        "gatk --java-options '-Xms15G -Xmx15G -XX:ParallelGCThreads=8' AddOrReplaceReadGroups" 
        " -I {input.bam}" 
        " -O {output.addRGbam}" 
        " -ID {params.ID}"
        " -LB {params.LB}"
        " -PL {params.PL}" 
        " -PU {params.PU}" 
        " -SM {params.SM}"


rule sambambaSortMarkdupIndex:
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        addRGbam = os.path.join(wkBamDir, "{sample}.RG.bam")
    output:
        sortbam = os.path.join(wkBamDir, "{sample}.sorted.bam"),
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.bam"),
        bamindex = os.path.join(wkBamDir, "{sample}.sorted.markdup.bam.bai")
    params: 
        tmpDir = os.path.join(wkBamDir, "{sample}_temp")
    threads:8
    shell:
        "mkdir -p {params.tmpDir} ;"
        "sambamba sort -t {threads} -m 20GB -p --tmpdir {params.tmpDir} -o {output.sortbam} {input.addRGbam} ;"
        "sambamba markdup -t {threads} -p --tmpdir {params.tmpDir} {output.sortbam} {output.markdupbam} ;"
        "sambamba index -t {threads} -p {output.markdupbam} {output.bamindex}"

rule BQSR:
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.bam"),
        dbSNP_vcf = "/home/resource/hg38/Homo_sapiens_assembly38.dbsnp138.vcf",
        known_indels_sites_vcf = "/home/reference/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz",
        Mills_and_1000G_gold_standard_indels_vcf = "/home/resource/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        # Maybe need a index:tbi for vcf.gz,idx for vcf  
    output:
        recalibration_report_filename = os.path.join(wkBamDir,"{sample}.sorted.markdup.BQSRrecalibration.table"),
    params:
        tmpDir = os.path.join(wkBamDir, "{sample}_temp")
    threads:8
    shell:
        ## BaseRecalibrator
        """
        mkdir -p {params.tmpDir} ;
        java -version ;
        java -Xms15G -Xmx15G -version ;
        /home/bin/java -Xms15G -Xmx15G -XX:ParallelGCThreads=8 \
            -jar /home/software/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar BaseRecalibrator \
            -I {input.markdupbam} \
            -R /home/reference/hg38/Homo_sapiens_assembly38.fasta \
            --use-original-qualities \
            -O {output.recalibration_report_filename} \
            --known-sites {input.dbSNP_vcf} \
            --known-sites {input.known_indels_sites_vcf} \
            --known-sites {input.Mills_and_1000G_gold_standard_indels_vcf} \
            --tmp-dir {params.tmpDir}
        """

rule ApplyBQSR:
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.bam"),
        recalibration_report_filename = os.path.join(wkBamDir,"{sample}.sorted.markdup.BQSRrecalibration.table"),
    output:
        BQSRbam = os.path.join(wkBamDir,"{sample}.sorted.markdup.BQSR.bam"),
        BQSRbambai = os.path.join(wkBamDir,"{sample}.sorted.markdup.BQSR.bai")
    params:
        tmpDir = os.path.join(wkBamDir, "{sample}_temp")
    threads:8
    shell:     
        ## ApplyBQSR
        """
        /home/bin/java -Xms15G -Xmx15G -XX:ParallelGCThreads=8 \
            -jar /home/software/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar ApplyBQSR \
            -I {input.markdupbam} \
            -R /home/reference/hg38/Homo_sapiens_assembly38.fasta \
            --use-original-qualities \
            --bqsr-recal-file {input.recalibration_report_filename} \
            -O {output.BQSRbam} \
            --tmp-dir {params.tmpDir} \
        """

rule GetSummary1:
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam"),
    output:
        flagstat = os.path.join(wkMultiQCDir,"{sample}.flagstat"),
    params:
        tmpDir = os.path.join(wkMultiQCDir, "{sample}_temp")
    threads:8
    shell:
       """ 
        mkdir -p {params.tmpDir} ;
        /home/samtools flagstat {input.markdupbam} > {output.flagstat}
       """

rule GetSummary2:
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam"),
    output:
        stats = os.path.join(wkMultiQCDir,"{sample}.stats"),
    params:
        tmpDir = os.path.join(wkMultiQCDir, "{sample}_temp")
    threads:8
    shell:
       """ 
        /home/samtools stats {input.markdupbam} > {output.stats}
       """

rule GetSummary3:
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam"),
    output:
        chrstats = os.path.join(wkMultiQCDir,"{sample}.chrstats")
    params:
        tmpDir = os.path.join(wkMultiQCDir, "{sample}_temp")
    threads:8
    shell:
       """ 
        /home/pandepth -i {input.markdupbam} -o {output.chrstats}
       """


rule Mutect2RMGerm:
    wildcard_constraints:
        sample = "[^\.]+"
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam"),
    output:
        mutect2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2RMGerm.vcf.gz")
    params:
        tmpDir = os.path.join(wkVCFDir, "{sample}_temp")
    threads:8
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
    input:
        markdupbam = os.path.join(wkBamDir, "{sample}.sorted.markdup.BQSR.bam"),
    output:
        mutect2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.vcf.gz"),
    params:
        tmpDir = os.path.join(wkVCFDir, "{sample}_temp")
    threads:8
    shell:
       """
          /home/bin/java -Xms15G -Xmx15G -XX:ParallelGCThreads=8 \
            -jar /home/software/gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar Mutect2 \
            -R /home/resource/reference/hg38/Homo_sapiens_assembly38.fasta \
            -I {input.markdupbam} \
            -I /Normal.sorted.markdup.BQSR.bam \
            --germline-resource /home/resource/hg38/af-only-gnomad.hg38.vcf.gz \
            --panel-of-normals /home/resource/hg38/1000g_pon.hg38.vcf.gz \
            -O {output.mutect2vcf}
       """

rule Filtered2Soma:
    wildcard_constraints:
           sample = "[^\.]+"
    input:
        mutect2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.vcf.gz")
    output:
        filtered2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.filtered.vcf"),
    params:
        tmpDir = os.path.join(wkVCFDir, "{sample}_temp")
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
    input:
        filtered2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.filtered.vcf")
    output:
        snpEff2vcf = os.path.join(wkVCFDir, "{sample}.Mutect2Soma.snpEff"),
    params:
        tmpDir = os.path.join(wkVCFDir, "{sample}_temp")
    threads:8
    shell:
       """
          genome="GRCh38.86"
          /home/bin/java -Xmx4g -jar /home/software/snpEff/snpEff.jar \
          -v $genome {input.filtered2vcf} > {output.snpEff2vcf}
        """
