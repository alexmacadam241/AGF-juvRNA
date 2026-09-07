#!/bin/bash
#PBS -o aten_trial2.log
#PBS -e aten_trial2.err

#Aligning reads using STAR is a two step process: 1. Create a genome index 2. Map reads to the genome
#STAR tutorial https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html 
#RNA=seq tutorial https://biocorecrg.github.io/RNAseq_course_2019/alnpractical.html
#another very useful tutorial: https://sydney-informatics-hub.github.io/training-RNAseq/03-MapReads/index.html

###I map porites reads that did no map to the host, to D. are

# Setting up input directory
INPUT_DIR=/home/jc828813/RNAseq/results/aten_juveniles/host/mapping/D1
#mkdir -p /home/jc489418/RNA_seq_porites_cliona_heat/11_star_mapping_relaxed/unmapped_host_to_UQ_Cladocopium_goreaui


#for each item in the folder assign each item temporarly into a variable DATASET after which do ..
#the basename command prints the last element of a file path
#-d "_" is the delimiter _
#-f is the field I want to display
#example file name: P1_12_S25_combinedLanes_R1_001.fastq.gz
for DATASET in $(ls $INPUT_DIR/*_aten_juveniles_Unmapped.out.mate1); do
PREFIX=$(basename $DATASET | cut -d "_" -f 1,2,3)
echo "for the file $DATASET, the prefix is $PREFIX"
echo "and the current file is $(basename $DATASET)"


#submitting pbs job - one job per sample
#-v list of variables is passed to the job - you also need to incude #PBS -V in the PBS script
qsub -v PBSDATASET="$DATASET" /home/jc828813/RNAseq/scripts/05.4_star_aten_unmapped_host_Dtre2.pbs
done