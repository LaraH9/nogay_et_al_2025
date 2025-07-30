#################################################################################
###
###  (A) diffbind_analysis.R to get log2FC values of Cut&Tag replicates
###
###   Per shared or unique (2x) peaks from MERGED bam files as well as 500 bp bin
###
#################################################################################

library(DESeq2)
library(tidyverse)
library(GenomicRanges)

### -----------------------------
### Parameters
count_threshold <- 25                                           ### -------- alter this accordingly  ------------
subfolder <- "norm_via_merged_scaling_flipped_canonical_chr"    ### -------- alter this accordingly  ------------

mods <- c("H3", "H3K27ac", "H3K27me3", "H3K9me3")
regions <- c("bin", "background", "shared", "ctr_UNIQUE", "cdk1KD_UNIQUE")
output_dir <- paste0(subfolder, "\\DESeq_",count_threshold, "_count_filter")
canonical_chr <- c("2L", "2R", "3L", "3R", "X", "Y")
### -----------------------------

# check if sub directory exists
if (file.exists(subfolder) == FALSE){
  print(paste0("Create directory: ", subfolder))
  dir.create(subfolder, showWarnings = FALSE)
}
if (file.exists(output_dir) == FALSE){
  print(paste0("Create directory: ", output_dir))
  dir.create(output_dir, showWarnings = FALSE)
}

### -----------------------------
# Provide spike-in factors in same order as your samples (rep1, rep2 for ctr and cdk1KD)

# from merged files       --------> Use this one, so replicates are normalized according to merged normalization factors
# [ctr_merged, ctr_merged, cdk1KD_merged, cdk1KD_merged]
spikein_factors_merged <- list(
  H3        = c( 0.999999, 0.999999, 1.03364,  1.03364),
  H3K27ac   = c( 1,               1, 0.937718, 0.937718),
  H3K27me3  = c( 0.999998, 0.999998, 1.31981,  1.31981),
  H3K9me3   = c( 1,               1, 1.3313,   1.3313)
)

#	Scaling factor from the normalization pipeline need to be flipped to obtain correct sizefactors 
for (mod in mods) {
    spikein_factors_merged[[mod]] <- 1 /  spikein_factors_merged[[mod]]
}
print(spikein_factors_merged)

# ---- Main loop ----
all_results <- list()

for (mod in mods) {
  
  message(glue::glue("\n--------------------------------------------\nProcessing: {mod}"))
  
  ### -----------------------------
  ### Identify bins that are outside of any peak region 
  ### (i.e., not overlapping any region from your broadPeak files for shared, ctr_UNIQUE, or cdk1KD_UNIQUE).
  # Load all bins
  all_bins <- read.delim(paste0("multiBamOutput_raw\\", mod, "_bin_counts.tsv"),
                       stringsAsFactors = FALSE)
  colnames(all_bins)[1:3] <- c("chr", "start", "end")
  
  # Filter for canonical chromosomes only
  all_bins <- all_bins %>% filter(chr %in% canonical_chr)
  
  # Load peak regions from counts files
  load_regions <- function(file, region_label) {
    df <- read.delim(file, stringsAsFactors = FALSE)
    colnames(df)[1:3] <- c("chr", "start", "end")
    df$region <- region_label
    # Filter for canonical chromosomes only
    df %>% filter(chr %in% canonical_chr)

  }
  df_shared     <- load_regions(paste0("multiBamOutput_raw\\", mod, "_shared_counts.tsv"), "shared")
  df_ctr_unique <- load_regions(paste0("multiBamOutput_raw\\", mod, "_ctr_UNIQUE_counts.tsv"), "ctr_UNIQUE")
  df_kd_unique  <- load_regions(paste0("multiBamOutput_raw\\", mod, "_cdk1KD_UNIQUE_counts.tsv"), "cdk1KD_UNIQUE")
  # Combine all peak bins
  peak_bins <- bind_rows(df_shared, df_ctr_unique, df_kd_unique)
  
  # Label bins as "background" only if they do NOT overlap 
  # with any of the peak intervals (shared or unique) on the same chromosome.
  # Convert bins to GRanges
  bin_gr <- GRanges(
    seqnames = all_bins$chr,
    ranges = IRanges(start = all_bins$start, end = all_bins$end)
  )
  # Convert peak regions to GRanges
  peak_gr <- GRanges(
    seqnames = peak_bins$chr,
    ranges = IRanges(start = peak_bins$start, end = peak_bins$end)
  )
  # Find overlaps
  overlaps <- findOverlaps(bin_gr, peak_gr)
  # Initialize region as "background"
  all_bins$region <- "background"
  # Assign overlapping regions from peak_bins
  matched_idx <- queryHits(overlaps)
  region_labels <- peak_bins$region[subjectHits(overlaps)]
  # If multiple peaks overlap one bin, choose the first (or use a rule like priority)
  all_bins$region[matched_idx] <- region_labels
  # extract background bins as single df
  background_df <- all_bins %>%  filter(region == "background") %>% select(-region)
  # Convert count columns to numeric (excluding chr)
  cols_to_convert <- colnames(background_df)[c(-1)]
  background_df[ , cols_to_convert] <- lapply(background_df[ , cols_to_convert], as.numeric)
  # save background bins a extra tsv file
  output_file <- paste0("multiBamOutput_raw\\", mod, "_background_counts.tsv")
  write.table(as.data.frame(background_df), file=output_file,
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  ### -----------------------------
  ### Go through all regions
  for (reg in regions) {

    message(glue::glue("\n----------------------\nProcessing: {mod} - {reg}"))
    
    # Load count data
    counts_file <- paste0("multiBamOutput_raw\\", mod, "_", reg, "_counts.tsv")
    if (reg == "background"){
      raw_counts <- read.table(counts_file, header=TRUE)
      raw_counts <- raw_counts %>% filter(chr %in% canonical_chr)
    } else {
      raw_counts <- read.table(counts_file, header=FALSE)
      raw_counts <- raw_counts %>% filter(V1 %in% canonical_chr)
    }
    counts_matrix <- as.matrix(raw_counts[ , -c(1:3)])  # remove first 3 columns if it's bin info
    colnames(counts_matrix) <- c("ctr_R1",  "ctr_R2", "cdk1KD_R1",  "cdk1KD_R2")
    print(type(counts_matrix))
    
    ##################################################
    ### Get REPLICATE MEANS
    
    # Calculate unnormalized means per row
    raw_counts$mean_ctr <- rowMeans(counts_matrix[, c("ctr_R1", "ctr_R2")], na.rm = TRUE)
    raw_counts$mean_cdk1KD <- rowMeans(counts_matrix[, c("cdk1KD_R1", "cdk1KD_R2")], na.rm = TRUE)
    
    # Normalize the matrix manually using flipped spike-in factors
    norm_factors <- spikein_factors_merged[[mod]]
    normalized_matrix <- sweep(counts_matrix, 2, norm_factors, FUN = "/")
    
    # Calculate normalized means per row
    raw_counts$norm_mean_ctr <- rowMeans(normalized_matrix[, c("ctr_R1", "ctr_R2")], na.rm = TRUE)
    raw_counts$norm_mean_cdk1KD <- rowMeans(normalized_matrix[, c("cdk1KD_R1", "cdk1KD_R2")], na.rm = TRUE)
    
    # Save these means to file
    raw_counts$modification <- mod
    raw_counts$region <- reg
    colnames(raw_counts) <- c("chr", "start", "end", "ctr_R1",  "ctr_R2", "cdk1KD_R1",  "cdk1KD_R2", 
                              "mean_ctr", "mean_cdk1KD", "norm_mean_ctr", "norm_mean_cdk1KD",
                              "Modification", "Region")
    write.table(raw_counts,paste0(subfolder, "/means_", mod, "_", reg, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    ##################################################
    ### DESeq PIPELINE
    # Create colData
    sample_info <- data.frame(
      row.names = colnames(counts_matrix),
      replicate = rep(c("R1", "R2"), 2),
      condition = rep(c("ctr", "cdk1KD"), each = 2)
    )
    sample_info$condition <- factor(sample_info$condition, levels = c("ctr", "cdk1KD"))
    
    # Build DESeq dataset
    dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                  colData = sample_info,
                                  design = ~ replicate + condition)
    
    # Apply custom spike-in normalization
    sizeFactors(dds) <- spikein_factors_merged[[mod]]

    # Filter low bins
    dds <- dds[rowMeans(counts(dds)) > count_threshold, ]
    
    # Run DESeq2
    dds <- DESeq(dds)
    res <- lfcShrink(dds, coef = "condition_cdk1KD_vs_ctr", type = "normal")
    
    # Annotate
    res_df <- as.data.frame(res)
    res_df$modification <- mod
    res_df$region <- reg
    # Save output
    output_file <- paste0(output_dir, "\\DESeq2_", mod,"_", reg, ".csv")
    write.csv(as.data.frame(res_df), file=output_file)
    
    all_results[[paste(mod, reg, sep = "_")]] <- res_df
  }
}

# ---- Save output ----
final_df <- bind_rows(all_results)
write.csv(final_df, paste0(output_dir, "\\All_log2FC_results.csv"), row.names = FALSE)
