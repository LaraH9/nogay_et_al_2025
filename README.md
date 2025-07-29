# Code for nogay_et_al_2025

Here we put the scripts used to analyse the CUT&Tag data. The scripts are example scripts that need to be adapted to the user path and data structure.

###  DNAmapping
We used snakePipes-v3.0.0 to align the fastq.gz files to a constructed hybrid dm6 and Lambda phage genome (Genbank: J02459.1). The config yaml file used to do the alignment is provided here.

###  Normalisation
The code used for the H3 and Lambda normalization was based on work done by Yinxiu Zhan (https://github.com/zhanyinx/atinbayeva_paper_2023) and modified according top our purposes.

###  Quantification
For the quantification, the number of PE reads were counted in each shared or unique peak for wt and cdk1KD (multiBamSummary BED-file --bamfiles <input files> --extendReads --outRawCounts <file>) as well as in 500 bp bins across the dm6 genome (multiBamSummary bins --bamfiles <input files> --binSize 500 --extendReads --outRawCounts <file>).
