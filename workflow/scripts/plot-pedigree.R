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
library(plot.matrix)

args <- commandArgs(trailingOnly=TRUE)

fam.file         <- args[1]
requested_output <- args[2]

requested_format   <- tools::file_ext(requested_output)
acceptable_formats <- c("pdf", "bmp", "svg", "png", "jpeg")

# Match the user-requested format withing the "acceptable" list, evaluate and store in a variable for later use.
# Kind reminder to self: If the answer implies using parse() or eval(), then the question should be rethinked...
format_function <- eval(as.symbol(acceptable_formats[match(requested_format, acceptable_formats)]))



fam.df <- read.table(fam.file, header=FALSE, sep=" ")
colnames(fam.df) <- c("Family.id", "Individual.id", "Father.id", "Mother.id", "Sex", "Phenotype")

# Convert 0 to <NA> values
fam.df$Father.id[fam.df$Father.id == "0"] <- NA
fam.df$Mother.id[fam.df$Mother.id == "0"] <- NA

fam.df$founder <- factor(
  ifelse(is.na(fam.df$Father.id) & is.na(fam.df$Mother.id), 0, 1),
  labels = c("founder", "simulated")
)

reps <- length(unique(fam.df$Family.id))


format_function(file = requested_output)
for (famid in unique(fam.df$Family.id)) {

  subset <- fam.df[fam.df$Family.id == famid, ]
  pedigree <- kinship2::pedigree(
    id    = subset$Individual.id,
    dadid = subset$Father.id,
    momid = subset$Mother.id,
    sex   = subset$Sex,
    affected = rep(1, nrow(subset)),
  )


  #par(mfrow = c(1,2))
  layout(matrix(c(1,1,1,2,2,3), nrow = 2, ncol = 3, byrow = TRUE))
  kinship2::plot.pedigree(
    pedigree,
    packed  = FALSE,
    cex     = 0.6,
    width   = 8,
    density = c(-1, 35, 65, 20),
    mar     = c(1, 1, 1, 1),
    symbolsize = 4,
    branch     = 1,
    angle      = c(90, 65, 40, 0),
    col        = c("cyan4", "darkorange3")[subset$founder],
    align      = F,
    lwd        = 2
  )

  mat <- kinship2::kinship(pedigree) * 2
  plot(mat, main = "Kinship matrix", las = 2)
  plot(factor(mat[lower.tri(mat)]), main="Relationship counts")
}
dev.off()

