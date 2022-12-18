# Host transcriptomics pipeline
Pipeline for analysing host transcriptomics data
_Last updated: Raphael Eisenhofer - 17/10/2022_

### General information:
This pipeline uses snakemake, and manages dependencies using conda (or mamba) for reproducibility and deployability. The 0_Code directory contains the snakefiles, scripts, and conda environment yamls. 

## Getting started:
Firstly, you'll want to clone this directory to the system where you want to run the analyses:
```
git clone https://github.com/EisenRa/Host_transcriptomics_pipeline.git
```

### What this pipeline does:
RAW reads have adapters trimmed and are quality filtered using Fastp, before having ribosomal reads removed using RiboDetector. These ribodepleted reads are then mapped to the host genome using STAR, which outputs the unmapped reads, sorted BAM file, and count table.

### Usage:
Currently, the snakefile searches for .fq.gz files located in this path (assuming you are launching the snakefile from the current directory):
```
2_Reads/1_Untrimmed/*_1.fq.gz
```
There are a couple of options for getting your data here:
- 1) Create symbolic links. This means you don't have to move or copy the files:
`ln -s reads_folder/*.fq.gz 2_Reads/1_Untrimmed/`
- 2) You can just put the reads in this directory.

(note that the fastq file suffixes should be **'_1.fq.gz'** and **'_2.fq.gz'**).

Next, you'll need to populate the absolute paths to the folder where the reference genome is stored, and the absolute path to the .gtf file. You must insert these into the following config file:
```
0_Code/Transcriptomics_config.yml
```

The default is currently the Gallus gallus 7 genome, but if you'd like to use a different target, point it to a new directory/reference genome and the snakefile will create the reference genome index for you.

That's all the setup required to get the pipeline running. Now you just need to launch the snakefile using snakemake. How you do this depends on your HPC server job queueing system. For Mjolnir, I use the following:
```
snakemake \
-s 0_Code/Transcriptomics.snakefile \
-j 10 \
--cluster "sbatch --mem {resources.mem_gb}G --time {resources.time} --cores {threads}" \
--use-conda \
--conda-frontend mamba \
--conda-prefix /projects/mjolnir1/people/ncl550/0_software \
--latency-wait 600
```

I recommend adding the `--dry-run` or `-n` command to the above code initially, as this will let you figure out if everything is working as expected.

I've written the pipeline such that it handles the requesting of optimised resources (RAM/CPUs) for each job based on the specific snakemake rule. The `--conda-prefix` snakemake option tells snakemake to look for conda environment in a particular directory, which saves having to reinstall them each time you run from a new directory.
