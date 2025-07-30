#!/bin/bash

## Author(s): Lara Heckmann
## Contact: heckmann@ie-freiburg.mpg.de
## This software is distributed without any guarantee under the terms of the GNU General

## This bash script runs spikein normalization using the "01_normalise.sh" script
## for all 4 modifications (including H3 as control 
## using replicates as data

# Define the histone modifications and conditions
histones=("H3" "H3K27ac" "H3K27me3" "H3K9me3")   
conditions=("ctr" "cdk1KD")
replicates=("R1" "R2")

# Redirect stdout (1) and stderr (2) to a log file
# Including timestamps
exec > >(while IFS= read -r line; do echo "$(date '+%Y-%m-%d %H:%M:%S') $line"; done | tee batch_H3_spikein_norm.log) 2>&1

# load required tools
module load samtools
module load deeptools

## H3
sh ../01_normalise.sh -r filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R1.filtered.bam -q filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R1.filtered.bam --refh3 filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R1.filtered.bam --queryh3 filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R1.filtered.bam --genomeid J02459.1_spikein
sh ../01_normalise.sh -r filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R2.filtered.bam -q filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R2.filtered.bam --refh3 filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R2.filtered.bam --queryh3 filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R2.filtered.bam --genomeid J02459.1_spikein


## H3K27ac
sh ../01_normalise.sh -r filtered_bam/*_FRA_CUTTAG_H3K27ac_WD_ctr_R1.filtered.bam -q filtered_bam/*_FRA_CUTTAG_H3K27ac_WD_cdk1KD_R1.filtered.bam --refh3 filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R1.filtered.bam --queryh3 filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R1.filtered.bam --genomeid J02459.1_spikein
sh ../01_normalise.sh -r filtered_bam/*_FRA_CUTTAG_H3K27ac_WD_ctr_R2.filtered.bam -q filtered_bam/*_FRA_CUTTAG_H3K27ac_WD_cdk1KD_R2.filtered.bam --refh3 filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R2.filtered.bam --queryh3 filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R2.filtered.bam --genomeid J02459.1_spikein


## H3K9me3
sh ../01_normalise.sh -r filtered_bam/*_FRA_CUTTAG_H3K9me3_WD_ctr_R1.filtered.bam -q filtered_bam/*_FRA_CUTTAG_H3K9me3_WD_cdk1KD_R1.filtered.bam --refh3 filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R1.filtered.bam --queryh3 filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R1.filtered.bam --genomeid J02459.1_spikein
sh ../01_normalise.sh -r filtered_bam/*_FRA_CUTTAG_H3K9me3_WD_ctr_R2.filtered.bam -q filtered_bam/*_FRA_CUTTAG_H3K9me3_WD_cdk1KD_R2.filtered.bam --refh3 filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R2.filtered.bam --queryh3 filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R2.filtered.bam --genomeid J02459.1_spikein


# H3K27me3
sh ../01_normalise.sh -r filtered_bam/*_FRA_CUTTAG_H3K27me3_WD_ctr_R1.filtered.bam -q filtered_bam/*_FRA_CUTTAG_H3K27me3_WD_cdk1KD_R1.filtered.bam --refh3 filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R1.filtered.bam --queryh3 filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R1.filtered.bam --genomeid J02459.1_spikein
sh ../01_normalise.sh -r filtered_bam/*_FRA_CUTTAG_H3K27me3_WD_ctr_R2.filtered.bam -q filtered_bam/*_FRA_CUTTAG_H3K27me3_WD_cdk1KD_R2.filtered.bam --refh3 filtered_bam/*_FRA_CUTTAG_H3_WD_ctr_R2.filtered.bam --queryh3 filtered_bam/*_FRA_CUTTAG_H3_WD_cdk1KD_R2.filtered.bam --genomeid J02459.1_spikein

echo "H3 normalization + spike-in normalization is done for all 4 modifications and 2 replicates."
