
DNAmapping --DAG --mapq 1 -i fastq/ -o Output_spikein --trim --fastqc --properPairs --dedup --cutntag --alignerOpts="--local --very-sensitive-local --no-discordant --no-mixed -I 10 -X 700" hybrid_adjusted.yaml 
