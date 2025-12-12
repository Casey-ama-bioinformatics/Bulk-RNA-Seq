# Bulk RNA Seq Workflow  
Casey K. Ho  
## Introduction
This page documents the basic workflow for bulk RNA sequencing data analysis at Eliceiri Lab, UC San Diego, School of Medicine. This is written in parallel with data processing for Book 50-10 Bulk RNA-Seq Data using Windows 11 Home with Hyper-V, `Docker Desktop`, and `RStudio`. 

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
cd C:<path to dir>

#Pulling Trim Galore in Docker Container (use most updated tag version)
docker pull quay.io/biocontainers/trim-galore:<version--tag>

#Running Trim Galore
docker run --rm -v
C:<working directory>:/data quay.io/biocontainers/trim-galore:<version--tag>
--paired
--nextera
--output_dir /data/<output dir for fastq files>
--fastqc_args "--outdir /data/<output dir for fastqc files>" /data/<sample>_R1_001.fastq /data/<sample>_R2_001.fastq

```
Refer to https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md for more mode options. 

I made sure that trimming of the adapters were complete by checking the 'FastQC' html reports before proceeding to sequence alignment. 

## 4a. STAR alignment + FeatureCounts
### STAR
STAR (Spliced Transcripts Alignment to a Reference) is tool which uses genome as a reference. This tool is used for genome level and splicing analysis. Before aligning the `fastq` files to the genome reference, a `star_index` is generated using your desired genome (.fa) and annotation (.gtf) references. For this analysis, I used the mouse primary assembly genome (GRCm38.primary_assembly.genome.fa) and mouse basic gene annotation (gencode.vM38.basic.annotation.gtf). These files can be downloaded from https://www.gencodegenes.org/mouse/. I used `star 2.7.11b` for index generation and running alignment. 

### Command
```Windows
#Set working directory to your reference files
cd C:<path to dir>

#Pulling STAR in Docker Container (use most updated tag version)
docker pull quay.io/biocontainers/star:<version--tag>

#Generating star_index
docker run --rm ^
  -v C:/<working directory>:/data ^
  quay.io/biocontainers/star:<version--tag> ^
  STAR --runMode genomeGenerate ^
       --genomeDir /data/star_index ^
       --genomeFastaFiles /data/GRCm38.primary_assembly.genome.fa ^ 
       --sjdbGTFfile /data/gencode.vM38.basic.annotation.gtf ^
       --genomeSAindexNbases 12 ^
       --runThreadN 8

#Running star alignment for all paired samples and moving all BAM files to a different folder
@echo off
setlocal enabledelayedexpansion

:: ================================
:: Define Directories
:: ================================
set fastqDir=C:<path to fastq files>
set outputDir=C:<path to output dir>
set indexDir=C:<path to star_index>
set bamDir=%outputDir%\BAMfiles
set logFile=%outputDir%\STAR_alignment_log.txt


:: ================================
:: Create output folders if missing
:: ================================
if not exist "%outputDir%" mkdir "%outputDir%"
if not exist "%bamDir%" mkdir "%bamDir%"

:: ================================
:: Initialize log file
:: ================================
echo ======================================== > "%logFile%"
echo STAR Alignment Log - %date% %time% >> "%logFile%"
echo Output directory: %outputDir% >> "%logFile%"
echo ======================================== >> "%logFile%"
echo. >> "%logFile%"

:: ================================
:: Loop through all R1 FASTQ files
:: ================================
for %%F in ("%fastqDir%\*_R1_001_val_1.fq") do (
    set "r1=%%~nxF"
    set "base=%%~nF"
    set "base=!base:_R1_001_val_1=!"
    set "r2=!base!_R2_001_val_2.fq"
    set "r2path=%fastqDir%\!r2!"


    echo.
    echo ========================================
    echo Checking sample: !base!
    echo R1: !r1!
    echo R2: !r2!
    echo ========================================


    echo [%date% %time%] Checking sample: !base! >> "%logFile%"


    :: Check if R2 exists
    if exist "!r2path!" (
        echo [%date% %time%] Starting STAR alignment for !base!... >> "%logFile%"


        docker run --rm ^
          -v "%fastqDir%:/fastq" ^
          -v "%outputDir%:/output" ^
          -v "%indexDir%:/index" ^
          quay.io/biocontainers/star:<version--tag> ^
          STAR --runThreadN 4 ^
               --outFilterMismatchNmax 2 ^
               --quantMode GeneCounts ^
               --genomeDir /index ^
               --readFilesIn /fastq/!r1! /fastq/!r2! ^
               --outFileNamePrefix /output/!base!_ ^
               --outSAMtype BAM SortedByCoordinate


        if %errorlevel%==0 (
            echo [%date% %time%] SUCCESS: STAR finished for !base! >> "%logFile%"
        ) else (
            echo [%date% %time%] ERROR: STAR failed for !base! >> "%logFile%"
        )
    ) else (
        echo WARNING: Missing R2 file for !base! — skipping this sample.
        echo [%date% %time%] WARNING: Missing R2 for !base! — skipped. >> "%logFile%"
    )
)


:: ================================
:: Move all BAM files to BAMfiles folder
:: ================================
echo.
echo ========================================
echo Moving BAM files to BAMfiles folder...
echo ========================================
move "%outputDir%\*.bam" "%bamDir%\" >nul 2>&1


echo [%date% %time%] BAM files moved to %bamDir% >> "%logFile%"
echo. >> "%logFile%"
echo ======================================== >> "%logFile%"
echo STAR alignment completed for all samples. >> "%logFile%"
echo Log saved at: %logFile% >> "%logFile%"
echo ======================================== >> "%logFile%"


echo.
echo ========================================
echo All STAR alignments completed!
echo BAM files consolidated into:
echo   %bamDir%
echo Log file created at:
echo   %logFile%
echo ========================================
pause
```
Refer to https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf for more details. 

### FeatureCounts
FeatureCounts is a tool to quantify the number of reads aligned to genes using a genomic annotation. I pulled `subread 2.0.3` in the `Docker` container for quantification.  

```Windows
@echo off
setlocal enabledelayedexpansion

:: === USER PATHS ===
set BAM_DIR=C:<path to BAM files>
set OUT_DIR=C:<path to output folder> 
set ANNOTATION=C:<path to genomic annotation> #.gtf file


:: === CREATE OUTPUT DIR IF NOT EXIST ===
if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"


:: === PULL FEATURECOUNTS DOCKER IMAGE ===
docker pull quay.io/biocontainers/subread:<version-tag>


:: === BUILD LIST OF BAM FILES ===
set BAM_LIST=
for %%F in ("%BAM_DIR%\*.bam") do (
    set "BAM_LIST=!BAM_LIST! /data/BAMfiles/%%~nxF"
)


:: === RUN FEATURECOUNTS ===
docker run --rm ^
    -v "%BAM_DIR%:/data/BAMfiles" ^
    -v "%OUT_DIR%:/data/output" ^
    -v "%ANNOTATION%:/data/annotation.gtf" ^
    quay.io/biocontainers/subread:<version-tag> ^
    featureCounts -a /data/annotation.gtf ^
                  -o /data/output/final_counts.txt ^
                  -g gene_id -T 4 -M --fraction -p !BAM_LIST!


echo.
echo === FeatureCounts Completed ===
echo Results saved in: %OUT_DIR%
pause
```
Refer to https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html for more details. 

## 4b. Salmon + Trimport
`Salmon` is an alignment tool utilizing quasi-mapping to align and quantify raw sequencing reads on a transcript level. Compared to `STAR`, `Salmon` is less resource heavy. The `Salmon + Trimport` pipeline is an alternative of `STAR + FeatureCounts`. I found that this pipeline yielded better alignment results as `STAR` ignores pseudogenes and genes with paralogs (e.g. Hbb genes), which might obscure gene expression quantification. Similar to `STAR`, Salmon requires the generation of a `salmon_index` folder using a transcriptome reference (GRCm38.primary_assembly.genome.fa) and annotation (gencode.vM38.basic.annotation.gtf). I pulled Salmon version `gffread 0.12.7`  in `Docker`. While generating `salmon_index` with the provided genome and annotation, `gffread` creates a filtered.gtf file (gencode.vM38.filtered.gtf) which ensures chromosomes in the annotation file matches our fastq files. A transcripts.fa file (gencode.vM38.transcripts.fa) will also be created in the process which will be the key input for `Salmon` alignment. 

```Windows

#Pull latest Salmon version in Docker
docker pull quay.io/biocontainers/gffread:<version-tag>

#Setting path to working directory where the filtered.gtf file will be written
docker run --rm -v "C:<path to working directory>:/ref" ubuntu bash -c "awk '$1 ~ /^chr[1-9XYM][0-9]*$/' /ref/<annotation.gtf > /ref/<filtered.gtf>"

#running GFFREAD on the filtered gtf file & generating transcripts
docker run --rm -v "C:<path to working directory>:/ref" quay.io/biocontainers/gffread:<version-tag> gffread /ref/<filtered.gtf> -g /ref/<reference genome.fa> -w /ref/<genome transcripts.fa>

#the salmon_index generation
# -k: k-mers (default 31 for mammalian transcriptome)
docker run --rm -v "C:<path to working directory>:/ref" -v "C:<path to working directory>\Salmon_index:/index" combinelab/salmon:latest salmon index -t /ref/<genome transcripts.fa> -i /index -k 31 -p 8

#Example paired alignment
docker run --rm -v C:<path to fastq files>:/ref -v C:<path to Salmon_index>:/index combinelab/salmon:latest salmon quant ^
  -i /index ^
  -l A ^
  -1 /ref/<path to example sample 1_read1.fq> ^
  -2 /ref/<path to example sample 1_read2.fq> ^
  -p 8 ^
  -o /ref/<output file>
```
After generating a `salmon_index` file, I ran paired-end alignment by looping through all samples. 

```Windows
@echo off

# Set folders and directories
set FASTQ_DIR=C:<path to fastq files>
set SALMON_INDEX=C:<path to Salmon_index>
set OUTPUT_DIR=C:<path to output directory>

if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

# Essential for changing variables in for loop
setlocal enabledelayedexpansion

# Loop through R1 files and find corresponding R2
for %%F in ("%FASTQ_DIR%/*_R1_001_val_1.fq") do (
    REM Get base filename without _R1_001_val_1.fq
    set "BASENAME=%%~nF"
    set "BASE=!BASENAME:_R1_001_val_1=!"


    # set output file for sample run 
    set "SAMPLE_OUT=/out/quant_!BASE!"


    echo ==============================================
    echo Processing sample !BASE!
    echo R1: /ref/trimmed_fastqfiles/!BASE!_R1_001_val_1.fq
    echo R2: /ref/trimmed_fastqfiles/!BASE!_R2_001_val_2.fq
    echo Output: !SAMPLE_OUT!


    # Run Salmon inside Docker
    docker run --rm -v C:<path to fastq files>:/ref -v C:<path to Salmon_index>:/index -v %OUTPUT_DIR%:/out combinelab/salmon:latest salmon quant ^
        -i /index ^
        -l A ^
        -1 /ref/<path to fastq files>/!BASE!_R1_001_val_1.fq ^
        -2 /ref/<path to fastq files>/!BASE!_R2_001_val_2.fq ^
        -p 8 ^
        -o !SAMPLE_OUT!
)
```
Please refer to https://salmon.readthedocs.io/en/latest/salmon.html for more documentation details. 

After running the `Salmon` alignment, there will be `quant.sf` files generated from the `Fastq` files in your directory. These files contain the transcript level expression of the samples and are ready for import into `R` for quantification using `Trimport` and data analysis using `DESeq2` package.

## Transcript Quantification with Trimport
We will first need to load the necessary libraries into the session. If they are not installed, install them using `install.packages()`. 
```R
library(tximport)
library(rtracklayer)
library(dplyr)
library(stringr)
library(tools)
```

Setting directories for the quant and sample files makes it easier to load the `quant.sf` files for downstream processing. A `gtf` annotation file is also required for `Trimport` which is the same as your annotation file (gencode.vM38.basic.annotation.gtf) that you used for Salmon alignment. The `gtf` annotation has the same function as a query for a database such that transcripts are annotated and quantified using `Trimport`. To ensure that all files are imported, check the objects with `head()` or `length()`. 

```R
#Setting directories
quant_dir <- "<path to quant files>"
gtf_file <- "<path to gtf annotation file>" 

#Importing Quant files and sample names
files <- list.files(quant_dir, pattern = "\\.sf$", full.names = TRUE)
sample_names <- file_path_sans_ext(basename(files))
names(files) <- sample_names
length(files)

#Importing gtf annotation
gtf <- rtracklayer::import(gtf_file)
gtf_df <- as.data.frame(gtf)

#Extract `gene_name` from `gtf_df`
if ("gene_name" %in% colnames(gtf_df)) {
  message("Using 'gene_name' column from GTF.")
  gtf_df$gene_name <- gtf_df$gene_name
} else if ("Name" %in% colnames(gtf_df)) {
  message("Using 'Name' column as gene_name.")
  gtf_df$gene_name <- gtf_df$Name
} else if ("attributes" %in% colnames(gtf_df)) {
  message("Extracting gene_name from 'attributes' column...")
  gtf_df$gene_name <- stringr::str_match(gtf_df$attributes, "gene_name \"([^\"]+)\"")[,2]
} else {
  stop("No suitable gene_name column found. Please inspect your GTF file.")
}
```
The object `gene_map` will store the gene-level counts of the samples. Subsequently, it will  be used to create a metadata with the samples. 

```R
#Create clean `gene_map`to map gene ids to gene symbols for the gtf annotation
gene_map <- gtf_df %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    gene_id_clean = stringr::str_remove(gene_id, "\\.\\d+$"),
    gene_symbol = gene_name
  ) %>%
  dplyr::select(gene_id_clean, gene_symbol)
message("Gene map created with ", nrow(gene_map), " entries.")
```

The purpose of creating `tx2gene` is to map the transcript ids present in the `quant.sf` files to the `gene_map` object such that they are annotated by their gene names. In other words, this summarizes the raw counts to gene level. 
```
# Create tx2gene Mapping 
transcripts <- gtf_df %>% filter(type == "transcript")
tx2gene <- transcripts %>%
  transmute(
    transcript_id = str_remove(transcript_id, "\\.\\d+$"),
    gene_id = str_remove(gene_id, "\\.\\d+$")
  ) %>%
  distinct()
message("tx2gene table created with ", nrow(tx2gene), " transcript-to-gene mappings.")

# Import Quant Files with `tximport`
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE
)
```
The preparation of the metadata `DESeq2 ColData` is important for `DESeq2` analysis as this affects how the counts are modeled across the samples, conditions, sex, tissues, and batches. In this analysis, there are two conditions (HFD, Lean) and two tissue types (Blood, PVA).

```
#Loading libraries
library(tibble)
library(dplyr)

# Create sample table
sample_table <- tibble(
  sample = colnames(txi$counts),
  condition = ifelse(grepl("HFD", sample), "HFD", "Lean"),
  tissue = ifelse(grepl("Blood", sample), "Blood", "PVA")
) %>%
  column_to_rownames(var = "sample")  # rownames = sample names

```
When converting the conditions and tissue into factors, `DESeq2` uses the first level as the baseline/control (lean, PVA) as this package uses a pairwise comparison between the average expression of two groups. 
```R
# Convert to factors with levels
sample_table$condition <- factor(sample_table$condition, levels = c("Lean", "HFD"))
sample_table$tissue <- factor(sample_table$tissue, levels = c("PVA","Blood"))

```
The function `DESeqDataSetFromTximport` constructs the object `dds` that will be used for differential expression analysis. 
```R
library(DESeq2)
library(dplyr)

dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ tissue + condition)

```

Prior to `DESeq2` analysis, genes with low counts across samples are filtered out. In this case, the minimum count of the gene is set to 10 (min_count=10)and the gene must be present in at least 4 samples(min_samples=4). I subsetted the tissues `Blood` and `PVA` from the object `dds` for differential expression analysis. Before running `DESeq2`, I had two subsetted objects `dds_blood` and `dds_pva`. 

```R
filter_low_counts <- function(dds, min_counts = 10, min_samples = 4) {
  keep <- rowSums(counts(dds) >= min_counts) >= min_samples
  dds_filtered <- dds[keep, ]
  message("Retained ", sum(keep), " genes out of ", nrow(dds), " total genes.")
  return(dds_filtered)
}
dds_filtered <- filter_low_counts(dds)


# Subsetting for blood (Repeat with PVA)
dds_blood <- dds[, dds$tissue == "Blood"]
dds_blood$condition <- droplevels(dds_blood$condition)
dds_blood$condition <- relevel(dds_blood$condition, ref = "Lean")
dds_blood <- filter_low_counts(dds_blood, min_counts = 10, min_samples = 4)
design(dds_blood) <- ~ condition
```

Surrogate Variable Analysis (SVA) is a normalization tool to correct batch variations. This increases the accuracy of the differential analysis and ensure that only biological variations are included. 
```R
#Normalization and batch correction for PVA samples (Repeat with dds_blood)
library(sva)

dds_pva <- dds[, dds$tissue == "PVA"]
dds_pva$condition <- droplevels(dds_pva$condition)
dds_pva$condition <- relevel(dds_pva$condition, ref="Lean")
dds_pva <- filter_low_counts(dds_pva, min_counts=10, min_samples=4)
dds_pva <- estimateSizeFactors(dds_pva)
normCounts_pva <- counts(dds_pva, normalized=TRUE)
mod <- model.matrix(~ condition, data=colData(dds_pva))
mod0 <- model.matrix(~ 1, data=colData(dds_pva))
svobj <- svaseq(normCounts_pva, mod, mod0)

for(i in seq_len(ncol(svobj$sv))) {
    colData(dds_pva)[[paste0("SV", i)]] <- svobj$sv[,i]
}

design(dds_pva) <- as.formula(paste("~", paste0("SV", 1:ncol(svobj$sv), collapse=" + "), "+ condition"))
```
The function `DESeq()` runs differential expression analysis on the objects. The results are extracted using the function `results()`.
```
#Run DESeq2 with normalization (Repeat with dds_blood)
dds_pva <- DESeq(dds_pva)
res_pva <- results(dds_pva, contrast=c("condition","HFD","Lean"))
res_pva <- lfcShrink(dds_pva, coef="condition_HFD_vs_Lean", type="apeglm")
head(counts(dds_pva, normalized=TRUE))

# Annotate the DESeq2 results with gene_map
annotate_res <- function(res) {
  res %>%
    mutate(gene_id_clean = str_remove(gene_id, "\\.\\d+$")) %>%
    left_join(gene_map, by = "gene_id_clean") %>%
    select(gene_id, gene_symbol, baseMean, log2FoldChange, padj)
}

res_blood_annotated <- annotate_res(res_blood)
res_pva_annotated   <- annotate_res(res_pva)

summary(res_blood_annotated) #summary() used to inspect the object.

#Optionally, save the annotated files with function write.csv(). 

## Extracting significant genes (padj < 0.05) 
sig_blood <- res_blood_annotated %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(padj)

sig_pva <- res_pva_annotated %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(padj)

```

















