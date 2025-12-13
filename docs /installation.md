# Package installation for data processing (Book 50-10)
Ensure Docker Desktop is running prior to package installation. 

## Sequencing
Illumina Nextera Transposase Adapaters

## Quality Control tools: FastQC and MultiQC
MultiQC version: ewels/multiqc

## Trimming tools: Trimgalore
Docker Trimgalore version: 0.6.7--hdfd78af_0
Trimming mode: Paired-end

## Alignment: STAR 
Docker STAR version: 2.7.11b  
Annotation: genome.vM38.basic.annotation.gtf  
Genome: GRCm38.primary_assembly.genome.fa  
Ensure STAR version pulled in index generation is the same or compatible with read alignment. 

## Alignment: Salmon
Docker Salmon Biocontainer version: gffread:0.12.7--h9a82719_0

## Feature Counts
Docker FeatureCounts Biocontainer version: subread:2.0.3--h7132678_1
