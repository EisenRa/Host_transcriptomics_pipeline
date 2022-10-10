############################################################################
#### Transcriptomics snakefile
#### Authors: Antton Alberdi/Raphael Eisenhofer
############################################################################

configfile: "0_Code/Transcriptomics_config.yml"

### Setup sample inputs
import os
from glob import glob

SAMPLE = [os.path.basename(fn).replace("_1.fastq.gz", "")
            for fn in glob(f"2_Reads/1_Untrimmed/*_1.fastq.gz")]

print("Detected the following samples:")
print(SAMPLE)

rule all:
    input:
        expand("3_Outputs/1_Mapping/{sample}_host.bam", sample=SAMPLE)

rule qualityfiltering:
    input:
        read1="2_Reads/1_Untrimmed/{sample}_1.fastq.gz",
        read2="2_Reads/1_Untrimmed/{sample}_2.fastq.gz"
    threads: 8
    resources:
        mem_gb=24
    conda: 
        "Transcriptomics_conda_env.yml"
    output:
        read1="2_Reads/2-Qualfilt/{sample}_1.fastq.gz",
        read2="2_Reads/2-Qualfilt/{sample}_2.fastq.gz",
        fastp_html="2_Reads/2-Qualfilt/{sample}.html",
        fastp_json="2_Reads/2-Qualfilt/{sample}.json"
    message:
        "Quality filtering {wildcards.sample} with fastp"
    shell:
        """
	    fastp \
            --in1 {input.read1} --in2 {input.read2} \
            --out1 {output.read1} --out2 {output.read2} \
            --trim_poly_g \
            --trim_poly_x \
            --n_base_limit 5 \
            --qualified_quality_phred 20 \
            --length_required 60 \
            --thread {threads} \
            --html {output.fastp_html} \
            --json {output.fastp_json} \
            --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            --adapter_sequence_r2  AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
        """

### Map to host reference genome using STAR
rule STAR_host_mapping:
    input:
        read1="2_Reads/2-Qualfilt/{sample}_1.fastq.gz",
        read2="2_Reads/2-Qualfilt/{sample}_2.fastq.gz",
    output:
        non_host_r1 = "3_Outputs/1_Mapping/{sample}_non_host_1.fastq.gz",
        non_host_r2 = "3_Outputs/1_Mapping/{sample}_non_host_2.fastq.gz",
        host_bam = "3_Outputs/1_Mapping/{sample}_host.bam"
    params:
        r1rn = "3_Outputs/1_Host_Mapping/{sample}_non_host_1.fastq",
        r2rn = "3_Outputs/1_Host_Mapping/{sample}_non_host_2.fastq",
        gene_counts = "3_Outputs/1_Mapping/{sample}_read_counts.tsv",
        sj = "3_Outputs/1_Mapping/{sample}_SJ.tsv",
        host_genome = expand("{host_genome}", host_genome=config['host_genome'])
    conda:
        "Transcriptomics_conda_env.yml"
    threads:
        24
    resources:
        mem_gb=150
    message:
        "Mapping {wildcards.sample} to the host genome using STAR"
    shell:
        """
        # Map reads to host genome using STAR
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.host_genome} \
            --readFilesIn {input.read1} {input.read2} \
            --outFileNamePrefix {wildcards.sample} \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --readFilesCommand zcat \
            --quantMode GeneCounts 

        # Rename files
        mv {wildcards.sample}Aligned.out.bam {output.host_bam}
        mv {wildcards.sample}ReadsPerGene.out.tab {params.gene_counts}
        mv {wildcards.sample}SJ.out.tab {params.sj}
        mv {wildcards.sample}Unmapped.out.mate1 {params.r1rn}
        mv {wildcards.sample}Unmapped.out.mate2 {params.r2rn}

        # Compress non-host reads
        pigz \
            -p {threads} \
            {params.r1rn}

        pigz \
            -p {threads} \
            {params.r2rn}

        # Clean up unwanted outputs
        rm -r {wildcards.sample}_STARtmp
        rm {wildcards.sample}*out*
        """
