#!/usr/bin/env Rscript

# Usage:   plot-pedigree.R <INPUT> <OUTPUT>
#            - INPUT  path. File must be a PLINK .fam format
#            - OUTPUT path. File must end with one of the "accepted" file formats
# 
#  - Accepted file formats: pdf bmp svg png jpeg
# 
# Example: plot-pedigree.R ./EUR-pedigrees-everyone.fam ./EUR-pedigreees-everyone.pdf
library(kinship2)
library(tools)

args <- commandArgs(trailingOnly=TRUE)

fam.file <- args[1]
requested_output <- args[2]
requested_format <- tools::file_ext(requested_output)

fam.df <- read.table(fam.file, header=FALSE, sep=" ")
colnames(fam.df) <- c("Family.id", "Individual.id", "Father.id", "Mother.id", "Sex", "Phenotype")

# Convert 0 to <NA> values
fam.df$Father.id[fam.df$Father.id == "0"] <- NA
fam.df$Mother.id[fam.df$Mother.id == "0"] <- NA
pedigree <- kinship2::pedigree(
  id    = fam.df$Individual.id,
  dadid = fam.df$Father.id,
  momid = fam.df$Mother.id,
  sex   = fam.df$Sex,
)


acceptable_formats <- c("pdf", "bmp", "svg", "png", "jpeg")

# Match the user-requested format withing the "acceptable" list, evaluate and store in a variable for later use.
# Kind reminder to self: If the answer implies using parse() or eval(), then the question should be rethinked...
format_function <- eval(as.symbol(acceptable_formats[match(requested_format, acceptable_formats)]))

# Plot
format_function(file = requested_output)
kinship2::plot.pedigree(pedigree, packed=F, cex=0.5, mar=c(5,5,5,5), branch=0)
dev.off()

