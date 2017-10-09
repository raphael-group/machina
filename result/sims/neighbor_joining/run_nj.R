#!/usr/bin/Rscript
suppressMessages(library(phangorn))

args <- commandArgs(trailingOnly = TRUE)
d <- read.table(args[1], header = T, row.names = "sample")
mat.nj = nj(dist.gene(d))
cat(paste0(write.tree(mat.nj), "\n"))
