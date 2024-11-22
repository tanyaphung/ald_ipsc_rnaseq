import os

configfile: "ipsc_config.json"

# SETTING UP PARAMETERS
fq_prefix = ["R1", "R2"]
trimmed_fq_prefix = ["read1", "read2"]
# --- DONE SETTING UP PARAMETERS --- #

samples = ["4388"]

rule all:
    input:
        "results/multiqc_rna/multiqc_report.html", #for rule multiqc_analysis_rna
        expand(os.path.join(config["out_directory"], "bbduk_trimmed_fastqs", "{sample}_trimmed_read1.fastq.gz"), sample=config["rna_samples"]), #for rule trim_adaptors_paired_bbduk
        expand("results/post_trimming_fastqc/{sample}_trimmed_read1_fastqc.html", sample=config["rna_samples"]), #for rule post_trimming_fastqc
        "results/post_trimming_multiqc/multiqc_report.html", #for rule multiqc_analysis_rna_posttrim
        expand("results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup.bam", sample=config["rna_samples"]),
        expand("results/default_reference/{sample}_default_featurecounts.tsv", sample=config["rna_samples"])

rule multiqc_analysis_rna:
	input:
		expand(
			"{fastq_directory}/QC/{read_group_identifier}_{prefix}_fastqc.html", fastq_directory=config["fastq_directory"], read_group_identifier=config["read_group_identifier"], prefix=fq_prefix)
	output:
		"results/multiqc_rna/multiqc_report.html"
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"multiqc --interactive -f -o results/multiqc_rna fastq_rna" #here, one needs to put in the correct path for fastq_rna

# adaptor and polyG trimming using bbduk
rule trim_bbduk:
	input:
		fq1 = os.path.join(config["fastq_directory"], "{sample}_S1_L001_R1.fastq.gz"),
		fq2 = os.path.join(config["fastq_directory"], "{sample}_S1_L001_R2.fastq.gz")
	output:
		out_fq1 = os.path.join(config["out_directory"], "bbduk_trimmed_fastqs", "{sample}_trimmed_read1.fastq.gz"),
		out_fq2 = os.path.join(config["out_directory"], "bbduk_trimmed_fastqs", "{sample}_trimmed_read2.fastq.gz")
	params:
		adapter = "bbmap-39.09-0/resources/adapters.fa" #here, put in the correct path to the adapters file
	shell:
		"bbduk.sh -Xmx3g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref={params.adapter} ktrim=r k=21 mink=11 trimpolygleft=30 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=55 maq=20"

rule post_trimming_fastqc: 
    input:
        fq1 = os.path.join(config["out_directory"], "bbduk_trimmed_fastqs", "{sample}_trimmed_read1.fastq.gz"),
        fq2 = os.path.join(config["out_directory"], "bbduk_trimmed_fastqs", "{sample}_trimmed_read2.fastq.gz")
    output:
        fq1_fastqc = "results/post_trimming_fastqc/{sample}_trimmed_read1_fastqc.html",
        fq2_fastqc = "results/post_trimming_fastqc/{sample}_trimmed_read2_fastqc.html"
    shell:
        """
        fastqc -o results/post_trimming_fastqc {input.fq1};
        fastqc -o results/post_trimming_fastqc {input.fq2};
        """

rule post_trimming_multiqc:
	output:
		"results/post_trimming_multiqc/multiqc_report.html"
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"multiqc --interactive -f -o results/post_trimming_multiqc results/post_trimming_fastqc/"

# MAPPING WITH HISAT2
# using the default reference
rule hisat2_default_reference_index:
	input:
		"Reference_Data/gencode_release46/GRCh38.p14.genome.fa"
	output:
		expand(
			"Reference_Data/gencode_release46/hisat2/GRCh38.p14.genome.{suffix}.ht2",
			suffix=[
				"1", "2", "3", "4", "5", "6", "7", "8"])
	shell:
		"hisat2-build {input} Reference_Data/gencode_release46/hisat2/GRCh38.p14.genome"

rule hisat2_map_reads_default_reference:
    input:
        fq1 = os.path.join(config["out_directory"], "bbduk_trimmed_fastqs", "{sample}_trimmed_read1.fastq.gz"),
        fq2 = os.path.join(config["out_directory"], "bbduk_trimmed_fastqs", "{sample}_trimmed_read2.fastq.gz")
    output:
        out = "results/default_reference/{sample}_GRCh38.p14.genome_default.sam"
    params:
        HISAT_Index_default = "Reference_Data/gencode_release46/hisat2/GRCh38.p14.genome"
    shell:
        "hisat2 -q --rna-strandness RF -p 8 -x {params.HISAT_Index_default} -1 {input.fq1} -2 {input.fq2} -S {output.out}"

# convert to bam and view stats
rule samtools_view:
    input:
        SAM = "results/default_reference/{sample}_GRCh38.p14.genome_default.sam"
    output:
        BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default.bam"
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule bam_sort:    
    input:
        IN_BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default.bam"
    output:
        sort_BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule MarkDups:
    input:
        sort_BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort.bam"
    output:
        BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup.bam",
        metrics = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup_metrics.txt"
    shell:
        "picard -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule AddReadGrps:
    input:
        Read_BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup.bam"
    output:
        BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup_rdgrp.bam"
    params:
        id = lambda wildcards: config[wildcards.sample]["ID"],
        sm = lambda wildcards: config[wildcards.sample]["SM"],
        lb = lambda wildcards: config[wildcards.sample]["LB"],
        pu = lambda wildcards: config[wildcards.sample]["PU"],
        pl = lambda wildcards: config[wildcards.sample]["PL"]
    shell:
        "picard -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule index_bam:
    input:
        BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup_rdgrp.bam"
    output:
        BAI = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup_rdgrp.bam.bai"
    shell:
        "bamtools index -in {input.BAM}"

rule stats_bam:
    input:
        BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup_rdgrp.bam"
    output:
        stats = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup_rdgrp_stats.txt"
    shell:
        "bamtools stats -in {input.BAM} > {output.stats}"

rule featurecounts_default:
    input:
        BAM = "results/default_reference/{sample}_GRCh38.p14.genome_default_sort_mkdup_rdgrp.bam",
        GTF = "Reference_Data/gencode_release46/gencode.v46.annotation.gtf"
    output:
        COUNTS = "results/default_reference/{sample}_default_featurecounts.tsv"
    params:
        THREADS = 5
    shell:
        """
        featureCounts -t exon -g gene_id -a {input.GTF} -o {output.COUNTS} {input.BAM}
        """