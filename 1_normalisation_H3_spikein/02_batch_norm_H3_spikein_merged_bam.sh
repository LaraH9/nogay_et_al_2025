#!/bin/bash

## Author(s): Lara Heckmann
## Contact: heckmann@ie-freiburg.mpg.de
## This software is distributed without any guarantee under the terms of the GNU General

## This bash script runs spikein normalization using the "01_normalise_merged.sh" script
## for all 4 antibodies (including H3 norm as well) and all 2 replicates

# Define the histone modifications and conditions
histones=("H3", "H3K27ac" "H3K27me3" "H3K9me3") 
conditions=("ctr" "cdk1KD")
replicates=("R1" "R2")

# Redirect stdout (1) and stderr (2) to a log file
# Including timestamps
exec > >(while IFS= read -r line; do echo "$(date '+%Y-%m-%d %H:%M:%S') $line"; done | tee batch_merged_H3_spikein_norm_BAM.log) 2>&1

# load required tools
module load samtools
module load deeptools

## H3
sh ../01_normalise_merged.sh -r filtered_bam/*_H3_WD_ctr.filtered.bam -q filtered_bam/*_H3_WD_cdk1KD.filtered.bam --refh3 filtered_bam/*_H3_WD_ctr.filtered.bam --queryh3 filtered_bam/*_H3_WD_cdk1KD.filtered.bam --genomeid J02459.1_spikein
# subtract control bw file from RNAi bw file
bigwigCompare -b1 filtered_bam/*_H3_WD_cdk1KD.filtered.bw -b2 filtered_bam/*_H3_WD_ctr.filtered.bw -o filtered_bam/2norm_H3_subtracted.bw --operation subtract

echo "H3 + spike-in normalization of bam files done.\nBigwig files subtration of control from RNAi done."

## H3K27ac
sh ../01_normalise_merged.sh -r filtered_bam/*_H3K27ac_WD_ctr.filtered.bam -q filtered_bam/*_H3K27ac_WD_cdk1KD.filtered.bam --refh3 filtered_bam/*_H3_WD_ctr.filtered.bam --queryh3 filtered_bam/*_H3_WD_cdk1KD.filtered.bam --genomeid J02459.1_spikein
# subtract control bw file from RNAi bw file
bigwigCompare -b1 filtered_bam/*_H3K27ac_WD_cdk1KD.filtered.bw -b2 filtered_bam/*_H3K27ac_WD_ctr.filtered.bw -o filtered_bam/2norm_H3K27ac_subtracted.bw --operation subtract

echo "H3 + spike-in normalization of bam files done.\nBigwig files subtration of control from RNAi done."

## H3K9me3
sh ../01_normalise_merged.sh -r filtered_bam/*_H3K9me3_WD_ctr.filtered.bam -q filtered_bam/*_H3K9me3_WD_cdk1KD.filtered.bam --refh3 filtered_bam/*_H3_WD_ctr.filtered.bam --queryh3 filtered_bam/*_H3_WD_cdk1KD.filtered.bam --genomeid J02459.1_spikein
# subtract control bw file from RNAi bw file
bigwigCompare -b1 filtered_bam/*_H3K9me3_WD_cdk1KD.filtered.bw -b2 filtered_bam/*_H3K9me3_WD_ctr.filtered.bw -o filtered_bam/2norm_H3K9me3_subtracted.bw --operation subtract

echo "H3 + spike-in normalization of bam files done.\nBigwig files subtration of control from RNAi done."

## H3K27me3
sh ../01_normalise_merged.sh -r filtered_bam/*_H3K27me3_WD_ctr.filtered.bam -q filtered_bam/*_H3K27me3_WD_cdk1KD.filtered.bam --refh3 filtered_bam/*_H3_WD_ctr.filtered.bam --queryh3 filtered_bam/*_H3_WD_cdk1KD.filtered.bam --genomeid J02459.1_spikein
# subtract control bw file from RNAi bw file
bigwigCompare -b1 filtered_bam/*_H3K27me3_WD_cdk1KD.filtered.bw -b2 filtered_bam/*_H3K27me3_WD_ctr.filtered.bw -o filtered_bam/2norm_H3K27me3_subtracted.bw --operation subtract

echo "H3 + spike-in normalization of bam files done.\nBigwig files subtration of control from RNAi done."

echo "H3 + spike-in normalization is done for all 4 modifications with merged bam files from both replicates."
