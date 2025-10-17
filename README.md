# Bulk RNA Seq Workflow  
Casey Ho  
## Introduction
This page documents the basic workflow for bulk RNA sequencing data analysis at Elicieri Lab, UC San Diego, School of Medicine. This is written in parallel with data processing for Book 50-10 Bulk RNA-Seq Data using Windows Command Prompt (Windows 11 Home License) with Hyper-v configuration and RStudio. 
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

## 3. Trimming of adapters
The samples were sequenced using the Illumina 
