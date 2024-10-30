#!/usr/bin/env Rscript

intersect_plink_frq <- function(freq = NULL, snp = NULL, bed = NULL) {
  frq <- read.table(freq, header = TRUE, sep = "")
  snp <- read.table(snp, header = FALSE, sep = "")
  bed <- read.table(bed, header = FALSE, sep = "\t")
  
  colnames(frq) <- c("CHR", "SNP", "A1", "A2", "MAF", "NCHROBS")
  colnames(snp) <- c("SNP", "chr", "length", "pos", "REF", "ALT")
  colnames(bed) <- c("chr", "start", "end", "REF", "ALT")

  # Get rsids which intersect the original .frq, the base .snp and
  # the provided .bed file
  rsid_intersect <- merge(
    x    = merge(x = frq, y = snp, by.x = "SNP", by.y = "SNP"),
    y    = bed,
    by.x = c("CHR", "pos"),
    by.y = c("chr", "end")
  )
  
  # Ensure the order remains the same by querying the original frq file.
  which(frq$SNP %in% rsid_intersect$SNP)
}

if (!interactive()) {
  args  <- commandArgs(trailingOnly = TRUE)
  rsids <- intersect_plink_frq(freq = args[1], snp = args[2], bed = args[3])
  write(readLines(args[1])[c(1, rsids + 1)], stdout())
}
