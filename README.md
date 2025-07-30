# Code for nogay_et_al_2025

Here we provide the scripts used to analyse the CUT&Tag data for the paper Nogay et al. 2025. The scripts are example scripts that need to be adapted to the user path and data structure.

###  0. DNAmapping
We used snakePipes-v3.0.0 to align the fastq.gz files to a constructed hybrid dm6 and Lambda phage genome (Genbank: J02459.1). The adjusted config yaml file used to do the alignment is provided [here](./0_mapping/hybrid_adjusted.yaml).

###  1. Normalisation
Aligned replicates were [merged](./1_normalisation_H3_spikein/00_merge_bam_files.sh) before normalization to both H3 and spike-in signals. The code used for the H3 and Lambda normalization was based on work done by [Yinxiu Zhan](https://github.com/zhanyinx/atinbayeva_paper_2023) and modified according top [our purposes](./1_normalisation_H3_spikein/02_batch_norm_H3_spikein_merged_bam.sh).

###  2. Peak Calling
We used snakePipes-v3.0.0 with MACS2-v 2.2.9.1 to call CUT&Tag peaks. Used config files are provided [here](./2_peak_calling).

###  3. Heatmaps & Profile Plot 
We used deepTools-v3.5.6 to compute matrices of signal enrichment +/− 5kb around transcription start sites (TSS) or peak centers of shared or unique peaks for wt and cdk1KD. Coverage heatmaps and profiles were created using plotHeatmap or plotProfile from deepTools-v.2.5.7. Used bash script is provided [here](./3_signal_visualization/visualize_signals_H3spikein.sh).

###  4. Quantification
For the quantification, the number of PE reads were counted in each shared or unique peak for wt and cdk1KD (multiBamSummary BED-file --bamfiles <input files> --extendReads --outRawCounts <file>) as well as in 500 bp bins across the dm6 genome (multiBamSummary bins --bamfiles <input files> --binSize 500 --extendReads --outRawCounts <file>). Bash script is provided [here](./4_quantification/00_get_raw_multibamsum_counts.sh). The [DESeq2-v1.38.3 analysis](./4_quantification/01_DESeq_merged_norm_per_bin_and_peaks.R) was performed using the previously computed H3 and spike-in scaling factors from merged data for size factor. Bins and peaks with basemean (average read count across samples) of less than or equal to 25 read counts were discarded. The design matrix was setup to compare the samples by condition and correct for replicate effects (design = ~replicate + condition). Finally, we executed the DESeq2 shrinkage of log2 fold changes (type = ‘normal’).
