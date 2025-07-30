# Code for nogay_et_al_2025

Here we put the scripts used to analyse the CUT&Tag data for the paper Nogay et al. 2025. The scripts are example scripts that need to be adapted to the user path and data structure.

###  DNAmapping
We used snakePipes-v3.0.0 to align the fastq.gz files to a constructed hybrid dm6 and Lambda phage genome (Genbank: J02459.1). The adjusted config yaml file used to do the alignment is provided [here](./0_mapping/hybrid_adjusted.yaml).

###  Normalisation
Aligned replicates were [merged](./1_normalisation_H3_spikein/00_merge_bam_files.sh) before normalization to both H3 and spike-in signals. The code used for the H3 and Lambda normalization was based on work done by [Yinxiu Zhan](https://github.com/zhanyinx/atinbayeva_paper_2023) and modified according top [our purposes](./1_normalisation_H3_spikein/02_batch_norm_H3_spikein_merged_bam.sh).

###  Quantification
For the quantification, the number of PE reads were counted in each shared or unique peak for wt and cdk1KD (multiBamSummary BED-file --bamfiles <input files> --extendReads --outRawCounts <file>) as well as in 500 bp bins across the dm6 genome (multiBamSummary bins --bamfiles <input files> --binSize 500 --extendReads --outRawCounts <file>).
