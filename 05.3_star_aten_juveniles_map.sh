#!/bin/bash

#Aligning reads using STAR is a two step process: 1. Create a genome index 2. Map reads to the genome
#STAR tutorial https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html 
#RNA=seq tutorial https://biocorecrg.github.io/RNAseq_course_2019/alnpractical.html
#another very useful tutorial: https://sydney-informatics-hub.github.io/training-RNAseq/03-MapReads/index.html

###I map A. tenuis reads to the genome of A. tenuis

# Setting up input directory
INPUT_DIR=/home/jc828813/RNAseq/samples/Aten_juv_samples/trimmed_reads/paired
mkdir -p /home/jc828813/RNAseq/results/aten_juveniles/host/mapping


#for each item in the folder assign each item temporarly into a variable DATASET after which do ..
#the basename command prints the last element of a file path
#-d "_" is the delimiter _
#-f is the field I want to display
#example file name: P1_12_S25_R1_001.trim.fastq.gz
for DATASET in $(ls $INPUT_DIR/*_R1_001.trim.fastq.gz); do
PREFIX=$(basename $DATASET | cut -d "_" -f 1,2,3)
echo "for the file $DATASET, the prefix is $PREFIX"
echo "and the current file is $(basename $DATASET)"

#submitting pbs job - one job per sample
#-v list of variables is passed to the job - you also need to incude #PBS -V in the PBS script
qsub -v PBSDATASET="$DATASET" /home/jc828813/RNAseq/scripts/05.3_star_aten_juveniles_map.pbs
done