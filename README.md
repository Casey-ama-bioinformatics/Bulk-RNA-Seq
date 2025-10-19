# Bulk RNA Seq Workflow  
Casey Ho  
## Introduction
This page documents the basic workflow for bulk RNA sequencing data analysis at Elicieri Lab, UC San Diego, School of Medicine. This is written in parallel with data processing for Book 50-10 Bulk RNA-Seq Data using Windows 11 Home with Hyper-V, `Docker Desktop`, and `RStudio`.
## Directory

## 1. Fastq files 
The fastq files are uncompressed using `7zip` and organized into `Fastqfiles` folder. 

## 2. Running FastQC on fastq files
Running FastQC analysis prior to trimming allows us to visualise sequencing metrics. This sets a baseline for post-processing metrics to ensure adaquate adapter trimming. Additional trimming might be required depending on sequencing quality. I downloaded the `FastQC` package from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ and ran the command in R. Alternatively, you can run 'FastQC' in Command Prompt directly. 

### R script
```R
# Set the working directory to the folder containing FastQC.
setwd(<path to FastQC>)

# Define the FastQC batch file name or path.
fastqc_exe <- "run_fastqc.bat"

# Specify input folder (.fastq files) directory. 
input_folder <- "<path to fastq files>"

# Specify FastQC output (HTML and ZIP files) directory. 
output_folder <- "<path to output folder>"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

file_list <- list.files(
  path = input_folder,
  pattern = "\\.(fastq|fq)$",
  full.names = TRUE
)

for (f in file_list) {
  
  cmd <- paste(
    shQuote(fastqc_exe),
    "--outdir", shQuote(output_folder),
    shQuote(f)
  )
  
  message("Executing command: ", cmd)
  
  system(cmd)
}

## Ref(Share): https://chatgpt.com/share/680499ed-1b60-8002-abf0-0dbf20ae2f7b
## Ref: https://chatgpt.com/g/g-p-67b46db7124c8191900f83011a09eeaa-charlene-fastqc/c/680488c2-2944-8002-bd30-0a3a75c74d3c

```
The FastQC documentation provides a detailed explanation for each sequencing metric. (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

## 3. Trimming with Trim Galore
Samples were sequenced on the NovaSeq x Plus with Illumina Stranded mRNA Library Kit. I used Trim Galore to trim the Nextera transposase adapters (Read 1 5′ TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG ; Read 2 5′ GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG). There are other trimming tools such as `Trimmomatic` and `Skewer` but I found that `Trim Galore` worked best in trimming paired end adapter sequences. I pulled `Trim Galore v0.6.7` in the `Docker` Container and ran the following code in Command Prompt.
### Command
```Windows
# Set the working directory to your fastq files
setwd(<path to dir>)

#Pulling Trim Galore in Docker Container (use most updated tag version)
docker pull quay.io/biocontainers/trim-galore:<tag>

#Running Trim Galore
docker run --rm -v
C:<working directory>:/data quay.io/biocontainers/trim-galore:<tag>
--paired
--nextera
--output_dir /data/<output dir for fastq files>
--fastqc_args "--outdir /data/<output dir for fastqc files>" /data/<sample>_R1_001.fastq /data/<sample>_R2_001.fastq

```
Refer to (https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) for various options. 
### Output
```Windows

```
I made sure that trimming of the adapters were complete by checking the 'FastQC' html reports before proceeding to sequence alignment. 

## 4a. STAR alignment + FeatureCounts
STAR (Spliced Transcripts Alignment to a Reference) is tool which uses genome as a reference. This tool is used for genome level and splicing analysis. Before aligning the `fastq` files to the genome reference, a `star_index` is generated using your desired genome (.fa) and annotation (.gtf) references. For this analysis, I used the mouse primary assembly genome (GRCm38.primary_assembly.genome.fa) and mouse basic gene annotation (gencode.vM38.basic.annotation.gtf). These files can be downloaded from https://www.gencodegenes.org/mouse/. 

### Command
```Windows
#Pulling STAR in Docker Container (use most updated tag version)
docker pull quay.io/biocontainers/trim-galore:<tag>

```
## 4b. Salmon




