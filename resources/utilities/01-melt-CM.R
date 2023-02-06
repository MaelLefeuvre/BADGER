#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

preds = read.table(args[1], header=TRUE, sep="\t")
truths = read.table("resources/ped-sim/ped-definition/pedigree-expected-results.tsv", header=FALSE, sep="\t")

truths$Pair_name_a = paste0("ped1_",truths$V1,"-","ped1_", truths$V2)
truths$Pair_name_b = paste0("ped1_",truths$V2,"-","ped1_", truths$V1)
merged_a = merge(preds, truths, by.x="Pair_name", by.y="Pair_name_a")
merged_b = merge(preds, truths, by.x="Pair_name", by.y="Pair_name_b")

output   = rbind(merged_a[,c("Most_Likely_rel", "V4")], merged_b[,c("Most_Likely_rel", "V4")])


cat(paste0(output$V4, " ", output$Most_Likely_rel,"\n"), sep = "")
