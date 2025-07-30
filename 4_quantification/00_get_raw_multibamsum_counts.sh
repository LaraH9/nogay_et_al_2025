#!/bin/bash

## Author(s): Lara Heckmann
## Contact: heckmann@ie-freiburg.mpg.de
## This software is distributed without any guarantee under the terms of the GNU General

### Calculated signal per replicate and condition (cdk1KD and ctr) for ...
## ... all genome bins with a bin size of 500 bp
## ... 3 different peaks from merged data peak calling (shared, 2x unique)

# Redirect stdout (1) and stderr (2) to a log file
# Including timestamps
exec > >(while IFS= read -r line; do echo "$(date '+%Y-%m-%d %H:%M:%S') $line"; done | tee logs/get_raw_multibamsum_counts.log) 2>&1

# load required tools
module purge
module load deeptools/3.5.6

## variables

modifications=("H3" "H3K27ac" "H3K27me3" "H3K9me3")
peak_types=("shared" "cdk1KD_UNIQUE" "ctr_UNIQUE")

## Required folders

BED_DIR="/data/processing/heckmann/Analysis_3498_Nogay_Classen/CUTandTag-seq_SPIKE-IN_hybrid/Output_spikein/MACS2"    # Folder with BED files
BAM_DIR="filtered_bam_original"   # Folder with raw BAMs
OUT_DIR="multiBamOutput_raw"

# Loop through mods and peak types

for mod in "${modifications[@]}"; do
  
  ## Run mulitbamsummary for all genome bins with a bin size of 500 bp
  
  echo "Running multiBamSummary on raw replicate bam files for ${mod} - 500 bp bins..."
  multiBamSummary bins \
    --bamfiles ${BAM_DIR}/*${mod}_WD_ctr_R1*.bam ${BAM_DIR}/*${mod}_WD_ctr_R2*.bam ${BAM_DIR}/*${mod}_WD_cdk1KD_R1*.bam ${BAM_DIR}/*${mod}_WD_cdk1KD_R2*.bam \
    --binSize 500 \
    --extendReads \
    --outRawCounts ${OUT_DIR}/${mod}_bin_counts.tsv \
    -o ${OUT_DIR}/${mod}_bin_summary.npz
   
  ### Run multibamsummary for all peak types (shared, 2x unique), which where calculated from merged normalized files
  
  for peak in "${peak_types[@]}"; do

    BED_FILE="${BED_DIR}/spikein_${mod}_${peak}.filtered.BAMPE.broadPeak"
    
    # Output files
    OUT_RAW="${OUT_DIR}/${mod}_${peak}_counts.tsv"
    OUT_NPZ="${OUT_DIR}/${mod}_${peak}_summary.npz"

    # Check BED exists
    if [[ ! -f "$BED_FILE" ]]; then
      echo "BED file missing: $BED_FILE"
      continue
    fi

    # Match BAMs
    BAMs=$(ls ${BAM_DIR}/*${mod}_WD_ctr_R1*.bam ${BAM_DIR}/*${mod}_WD_ctr_R2*.bam ${BAM_DIR}/*${mod}_WD_cdk1KD_R1*.bam ${BAM_DIR}/*${mod}_WD_cdk1KD_R2*.bam)

    # Skip if no BAMs found
    if [[ -z "$BAMs" ]]; then
      echo "No BAMs found for ${mod}"
      continue
    fi

    echo "Running multiBamSummary on raw replicate bam files for ${mod} - ${peak}..."
    multiBamSummary BED-file \
      --BED "$BED_FILE" \
      --bamfiles ${BAM_DIR}/*${mod}_WD_ctr_R1*.bam ${BAM_DIR}/*${mod}_WD_ctr_R2*.bam ${BAM_DIR}/*${mod}_WD_cdk1KD_R1*.bam ${BAM_DIR}/*${mod}_WD_cdk1KD_R2*.bam \
      --extendReads \
      --numberOfProcessors 8 \
      --outRawCounts "$OUT_RAW" \
      -o "$OUT_NPZ"

  done

done
