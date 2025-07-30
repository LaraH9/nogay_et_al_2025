#!/bin/bash

## Author(s): Lara Heckmann
## Contact: heckmann@ie-freiburg.mpg.de
## This software is distributed without any guarantee under the terms of the GNU General

## This bash script merge 8x both replicates of WT or KD and all 4 modifications
## Change input (1) normalization and (2) directory in bash argruments - according to whether bam is normalized to RPKM or spike-in

## Use tmux and slurm module to run it
# For that you might first load the slurm module with "module load slurm".
# Then run 'sbatch merge_bam_files.sh NORMALIZATION NORM_BAMCOVERAGE DIRECTORY'
# You can check your job with "squeue" and cancel it with "scancel".

# Check if the user provided a normalization method
if [ $# -eq 0 ]; then
    echo "Error: Please provide a normalization method (e.g., RPKM, CPM, BPM)."
    exit 1
fi

# Assign the normalization method from the first command-line argument
normalization=$1
# Assign bamCoverage parameter either [RPKM, CPM, BPM, RPGC, None]
normalization_bamCoverage=$2
# Assign the input directory of bam files from second command-line argument
input_dir=$3

# Define the histone modifications and conditions
histones=("H3" "H3K27ac" "H3K27me3" "H3K9me3")
conditions=("ctr" "cdk1KD")

# Redirect stdout (1) and stderr (2) to a log file
# Including timestamps
exec > >(while IFS= read -r line; do echo "$(date '+%Y-%m-%d %H:%M:%S') $line"; done | tee logs/${normalization}_merge_bam_files.log) 2>&1

# load required tools
module load samtools
module load deeptools

# create output directory
mkdir -p merged_bam

# Loop through each histone and condition
for histone in "${histones[@]}"; do
    for condition in "${conditions[@]}"; do
        
        # Create the output file name
        output_bam="merged_bam/${normalization}_${histone}_WD_${condition}.filtered.bam"

        # Create the input file pattern
        input_bams="${input_dir}/*_FRA_CUTTAG_${histone}_WD_${condition}_R*.filtered.bam"

        # Run samtools merge
        echo "\nMerging files for ${histone} ${condition} with normalization ${normalization} into ${output_bam}..."
        
        samtools merge -f -@ 4 "$output_bam" $input_bams

        # Check if the command succeeded
        if [ $? -eq 0 ]; then
            echo "Successfully merged files for ${histone} ${condition} with normalization ${normalization} + ${normalization_bamCoverage}."
            
            # Index the merged BAM file
            echo "Indexing ${output_bam}..."
            samtools index "$output_bam"
            
            if [ $? -eq 0 ]; then
                echo "Successfully indexed ${output_bam}."
            
                # Generate normalized bigWig file
                output_bw="bamCoverage/${normalization}_${histone}_WD_${condition}.filtered.seq_depth_norm.bw"
                echo "Generating normalized bigWig file for ${output_bam}..."

                bamCoverage -b "$output_bam" -o "$output_bw" --normalizeUsing "$normalization_bamCoverage" --binSize 25  --numberOfProcessors 8

                if [ $? -eq 0 ]; then
                    echo "Successfully generated normalized bigWig file: ${output_bw}."
                else
                    echo "Error generating normalized bigWig file for ${output_bam}."
                fi
            
        else
                echo "Error indexing ${output_bam}."
            fi
            
        else
            echo "Error merging files for ${histone} ${condition} with normalization ${normalization}."
        fi

    done
done

echo "All merging, indexing, and bigWig generation operations completed with normalization ${normalization} + ${normalization_bamCoverage}."
