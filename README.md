## Code documentation
1. conda environment
- the packages used in this project is stored in the file `packages.txt`
- to use this spec file: 
```
conda create --name myenv --file spec-file.txt
```
2. RNAseq processing
- a snakemake pipeline is provided in the file `process_rna_shared.smk` (relative paths are included)
- config file for the snakemake pipeline: `ipsc_config_shared.json` (relative paths are included)
**Important notes:** all absolute paths in these two files have been removed. Therefore, in order to run the pipeline, one would need to supply the correct path
- The snakemake pipeline includes the following rules for RNAseq processing: 
    - aggregate fastqc using multiqc
        - rule: `multiqc_analysis_rna`
        - here, since fastqc was done by the sequencing vendor, multiqc was run to aggregate the results
    - trim adaptors and polyG
        - rule: `trim_bbduk`
    - after trimming, examine quality
        - rule: `post_trimming_fastqc`
        - rule: `post_trimming_multiqc`
    - mapping to a reference genome
        - Download reference data
            - Download the file: Nucleotide sequence of the GRCh38.p14 genome assembly version on all regions, including reference chromosomes, scaffolds, assembly patches and haplotypes
            ```
            # working directory: /gpfs/work5/0/vusr0480/Reference_Data/gencode_release46/
            wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz
            ```
            - Download the file: Comprehensive gene annotation
            ```
            # working directory: /gpfs/work5/0/vusr0480/Reference_Data/gencode_release46/
            wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
            ```
        - rule: `hisat2_default_reference_index`
        - rule: `hisat2_map_reads_default_reference`
    - after mapping, sort bam, mark duplications, add readgroups, and index
        - rule: `samtools_view` to convert sam to bam
        - rule: `bam_sort`: sort bam file
        - rule: `MarkDups`: to mark diplicates using picard
        - rule: `AddReadGrps`: to add read groups
            - to add readgroup to the config file, make sure to update the script to have the correct path, then run
            ```
            python generate_readgroup_config.py
            ```
        - rule: `index_bam`
        - rule: `stats_bam`
    - read quantification using featureCounts
        - rule: `featurecounts_default`
- After read quantification: 
    - Combine counts from 36 samples
    ```
    python combine_count_shared.py
    ```

    - Get the genes only
    ```
    grep -v "#" 4387_default_featurecounts.tsv | grep -v "Geneid" | awk '{print$1}' > transcripts.txt
    ```

- Then, generate a DGEobject using the counts. Subsequently, MDS plots and DGE analyses were performed. These were stored in the folder `DGE_analyses`, which will be made available before publication. 

3. Additional analyses: genotype calling using the mapped bam files for the gene ABCD1. 
- Commands:
```
gatk_path="bin/gatk3"
ref_fasta="Reference_Data/gencode_release46/GRCh38.p14.genome.fa"

i=$sample_id #make a list of sample ids here

input_bam="results/default_reference/${i}_GRCh38.p14.genome_default_sort_mkdup_rdgrp.bam"

# samtools index
samtools index ${input_bam}

# call variants with GATK Haplotype caller
${gatk_path} \
    -T HaplotypeCaller \
    -R ${ref_fasta} \
    -I ${input_bam} \
    -L chrX \
    -o 
    results/default_reference/${i}_GRCh38.p14.genome_default.vcf.gz \
    -U ALLOW_N_CIGAR_READS \

# vcftools subset
vcftools --gzvcf 
results/default_reference/${i}_GRCh38.p14.genome_default.vcf.gz --chr chrX --from-bp 153724856 --to-bp 153744755 --out 
results/default_reference/${i}_GRCh38.p14.genome_default_ABCD1 --recode
```

4. Additional analyses: compute gene length
- Commands: 
```
# compute gene length
for i in {1..22} X Y M; do

python convert_gtf_to_bed.py --input_gtf Reference_Data/gencode_release46/gencode.v46.annotation.gtf --output_bed chr${i}_info.txt --chromosome chr${i}; done;

python compute_genelength.py --input chr${i}_info.txt --uniqu
e_genes chr${i}_info_uniqgenes.txt --output chr${i}_genelength.txt
done;

# concat for all chr
for i in {1..22} X Y M; do
cat chr${i}_genelength.txt >> allChr_genelength.txt
done;
```