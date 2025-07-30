#!/bin/bash

## Author(s): Lara Heckmann
## Contact: heckmann@ie-freiburg.mpg.de
## This software is distributed without any guarantee under the terms of the GNU General

## This bash script runs deeptools like computeMatrix, plotheatmap, plotProfile for all 4 modifications
## including TSS, center around unqiue or shared peaks of wt and cdk1KD.

## Change input (1) normalization [spike-in, RPKM,...] - according to whether bam is normalized to RPKM or spike-in
## USAGE: bash visualize_signals_H3spikein.sh spikein /PATH/TO/dm6_ensembl.gtf

# Check if the user provided a normalization method
if [ $# -eq 0 ]; then
    echo "Error: Please provide a normalization method (e.g., RPKM, spike-in)."
    exit 1
fi

# Assign the normalization method from the first command-line argument
normalization=$1
gtf_file=$2          # e.g. /data/repository/organisms/dm6_ensembl/Ensembl/release-96/genes.gtf

# Define the histone modifications and conditions
histones=("H3" "H3K27ac" "H3K27me3" "H3K9me3")
colors=("Blues" "Oranges" "Greens" "Reds")
profile_colors=("blue" "orange" "green" "red")
conditions=("ctr" "cdk1KD")

# Redirect stdout (1) and stderr (2) to a log file
# Including timestamps
exec > >(while IFS= read -r line; do echo "$(date '+%Y-%m-%d %H:%M:%S') $line"; done | tee logs/${normalization}_visualize_heatmaps.log) 2>&1

# If necessary create directory
mkdir -p deepTools_qc/computeMatrix
mkdir -p deepTools_qc/plotHeatmap
mkdir -p deepTools_qc/plotProfile

# specify TSS BED file with TSS annotations from /PATH/TO/dm6_ensembl.gtf
# Extract TSS coordinates for positive and negative strand genes and add them
awk '$3 == "gene" && $7 == "+" {print $1 "\t" $4 "\t" $4+1 "\t" $9}' ${gtf_file} | awk -F'gene_id "|";' '{print $1 "\t" $2 "\t" $3 "\t" $4}' > dm6_release-96_genes_TSS.bed
awk '$3 == "gene" && $7 == "-" {print $1 "\t" $5-1 "\t" $5 "\t" $9}' ${gtf_file} | awk -F'gene_id "|";' '{print $1 "\t" $2 "\t" $3 "\t" $4}' >> dm6_release-96_genes_TSS.bed

# load required tools
module purge
module load deeptools/3.5.6 
module load bedtools

# Loop through each histone and condition
for i in "${!histones[@]}"; do
    histone=${histones[$i]}
    color=${colors[$i]}
    profile_color=${profile_colors[$i]}

    ### Get SHARED or UNIQUE PEAKS
    
    # Get shared peak of WT (ctr) and cdk1KD
    echo " "
    echo "Compute shared peak of ${normalization}_${histone}_*${conditions[1]}_host.BAMPE*broadPeak and ${normalization}_${histone}_*${conditions[0]}.filtered.BAMPE*broadPeak"
    bedtools intersect -a MACS2/${normalization}_${histone}_*${conditions[1]}_host.BAMPE*broadPeak -b MACS2/${normalization}_${histone}_*${conditions[0]}_host.BAMPE*broadPeak > MACS2/${normalization}_${histone}_shared.filtered.BAMPE.broadPeak 
   
    # Get unique peaks for WT and cdk1KD, respectively
    echo "Compute shared peak of ${normalization}_${histone}_*${conditions[0]}.filtered.BAMPE*broadPeak"
    bedtools subtract -a MACS2/${normalization}_${histone}_*${conditions[0]}_host.BAMPE*broadPeak -b MACS2/${normalization}_${histone}_*${conditions[1]}_host.BAMPE*broadPeak > MACS2/${normalization}_${histone}_${conditions[0]}_UNIQUE.filtered.BAMPE.broadPeak 

    echo "Compute shared peak of ${normalization}_${histone}_*${conditions[1]}.filtered.BAMPE*broadPeak"
    bedtools subtract -a MACS2/${normalization}_${histone}_*${conditions[1]}_host.BAMPE*broadPeak -b MACS2/${normalization}_${histone}_*${conditions[0]}_host.BAMPE*broadPeak > MACS2/${normalization}_${histone}_${conditions[1]}_UNIQUE.filtered.BAMPE.broadPeak


    ### COMPUTE MATRICES
    
    # center aound SHARED peaks 
    echo "CHECK FILES:"
    ls -l bamCoverage/spikein_${histone}_WD_ctr.filtered.seq_depth_norm.bw
    ls -l bamCoverage/spikein_${histone}_WD_cdk1KD.filtered.seq_depth_norm.bw
    
    # center aound SHARED peaks 
    echo "Doing computeMatrix of shared peaks for ${histone} with normalization ${normalization} into deepTools_qc/computeMatrix..."
    computeMatrix reference-point --referencePoint center -S bamCoverage/spikein_${histone}_WD_ctr.filtered.seq_depth_norm.bw bamCoverage/spikein_${histone}_WD_cdk1KD.filtered.seq_depth_norm.bw -R MACS2/${normalization}_${histone}_shared.filtered.BAMPE.broadPeak -b 5000 -a 5000 -o deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_matrix.gz --missingDataAsZero

    # center aound UNIQUE peaks of CDK1KD
    echo "Doing computeMatrix of UNIQUE cdk1KD peaks for ${histone} with normalization ${normalization} into deepTools_qc/computeMatrix..."
    computeMatrix reference-point --referencePoint center -S bamCoverage/spikein_${histone}_WD_ctr.filtered.seq_depth_norm.bw bamCoverage/spikein_${histone}_WD_cdk1KD.filtered.seq_depth_norm.bw -R MACS2/${normalization}_${histone}_cdk1KD_UNIQUE.filtered.BAMPE.broadPeak -b 5000 -a 5000 -o deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_cdk1KD_UNIQUE_matrix.gz --missingDataAsZero

    # center aound UNIQUE peaks of CONTROL
    echo "Doing computeMatrix of UNIQUE CONTROL peaks for ${histone} with normalization ${normalization} into deepTools_qc/computeMatrix..."
    computeMatrix reference-point --referencePoint center -S bamCoverage/spikein_${histone}_WD_ctr.filtered.seq_depth_norm.bw bamCoverage/spikein_${histone}_WD_cdk1KD.filtered.seq_depth_norm.bw -R MACS2/${normalization}_${histone}_ctr_UNIQUE.filtered.BAMPE.broadPeak -b 5000 -a 5000 -o deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_ctr_UNIQUE_matrix.gz --missingDataAsZero
    
    # TSS around genes with created TSS annoation bed file
    echo "Doing computeMatrix of genes for ${histone} with normalization ${normalization} into deepTools_qc/computeMatrix..."
    computeMatrix reference-point --referencePoint TSS -S bamCoverage/spikein_${histone}_WD_ctr.filtered.seq_depth_norm.bw bamCoverage/spikein_${histone}_WD_cdk1KD.filtered.seq_depth_norm.bw -R dm6_release-96_genes_TSS.bed -b 5000 -a 5000 -o deepTools_qc/computeMatrix/${normalization}_${histone}_GENE_matrix.gz --missingDataAsZero


   ### Plot HEATMAP

   # Plot coverage heatmap around shared peaks
   echo "Plot peak heatmap for ${histone} with normalization ${normalization} into deepTools_qc/plotHeatmap..."
   plotHeatmap -m deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_matrix.gz -out deepTools_qc/plotHeatmap/${normalization}_${histone}_PEAK_heatmap.png --colorMap ${color} --samplesLabel WT cdk1KD --regionsLabel ${histone},\ shared\ peaks --xAxisLabel ' ' --whatToShow heatmap\ and\ colorbar
   
   # Plot coverage heatmap around UNIQUE cdk1KD peaks
   echo "Plot peak heatmap around UNIQUE cdk1KD peaks for ${histone} with normalization ${normalization} into deepTools_qc/plotHeatmap..."
   plotHeatmap -m deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_cdk1KD_UNIQUE_matrix.gz -out deepTools_qc/plotHeatmap/${normalization}_${histone}_PEAK_cdk1KD_UNIQUE_heatmap.png --colorMap ${color} --samplesLabel WT cdk1KD --regionsLabel ${histone},\ unique\ cdk1KD\ peaks --xAxisLabel ' ' --whatToShow heatmap\ and\ colorbar

   # Plot coverage heatmap around UNIQUE control peaks
   echo "Plot peak heatmap around UNIQUE CONTROL peaks for ${histone} with normalization ${normalization} into deepTools_qc/plotHeatmap..."
   plotHeatmap -m deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_ctr_UNIQUE_matrix.gz -out deepTools_qc/plotHeatmap/${normalization}_${histone}_PEAK_ctr_UNIQUE_heatmap.png --colorMap ${color} --samplesLabel WT cdk1KD --regionsLabel ${histone},\ unique\ WT\ peaks --xAxisLabel ' ' --whatToShow heatmap\ and\ colorbar

   # Plot coverage heatmap around TSS genes
   echo "Plot gene heatmap for ${histone} with normalization ${normalization} into deepTools_qc/plotHeatmap..."
   plotHeatmap -m deepTools_qc/computeMatrix/${normalization}_${histone}_GENE_matrix.gz -out deepTools_qc/plotHeatmap/${normalization}_${histone}_GENE_heatmap.png --colorMap ${color} --samplesLabel WT cdk1KD --regionsLabel ${histone}\ peaks\ around\ TSS --xAxisLabel ' ' --whatToShow heatmap\ and\ colorbar

   ### Plot PROFILE

   # Plot profile around peaks
   echo "Plot peak profile for ${histone} with normalization ${normalization} into deepTools_qc/plotProfile..."
   plotProfile -m deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_matrix.gz -out deepTools_qc/plotProfile/${normalization}_${histone}_PEAK_profile.png --colors black ${profile_color} --regionsLabel '' --samplesLabel WT cdk1KD --perGroup
 
   # Plot profile around UNIQUE cdk1KD peaks
   echo "Plot peak profile around UNIQUE cdk1KD peaks for ${histone} with normalization ${normalization} into deepTools_qc/plotProfile..."
   plotProfile -m deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_cdk1KD_UNIQUE_matrix.gz -out deepTools_qc/plotProfile/${normalization}_${histone}_PEAK_cdk1KD_UNIQUE_profile.png --colors black ${profile_color} --regionsLabel '' --samplesLabel WT cdk1KD --perGroup 

   # Plot profile around UNIQUE control peaks
   echo "Plot peak profile around UNIQUE CONTROL peaks for ${histone} with normalization ${normalization} into deepTools_qc/plotProfile..."
   plotProfile -m deepTools_qc/computeMatrix/${normalization}_${histone}_PEAK_ctr_UNIQUE_matrix.gz -out deepTools_qc/plotProfile/${normalization}_${histone}_PEAK_ctr_UNIQUE_profile.png --colors black ${profile_color} --regionsLabel '' --samplesLabel WT cdk1KD --perGroup 

   # Plot profile around TSS genes
   echo "Plot gene profile for ${histone} with normalization ${normalization} into deepTools_qc/plotProfile..."
   plotProfile -m deepTools_qc/computeMatrix/${normalization}_${histone}_GENE_matrix.gz -out deepTools_qc/plotProfile/${normalization}_${histone}_GENE_profile.png --colors black ${profile_color} --regionsLabel '' --samplesLabel WT cdk1KD --perGroup 

done

echo "All computeMatrix, plotheatmap, plotProfile done with ${normalization} + ${normalization_bamCoverage} data."
