#!/usr/bin/Rscript
library(sciClone)

args = commandArgs(TRUE)
dirname = args[1]
clusters = as.integer(args[2])
samples = gsub(".tsv", "", list.files(dirname, pattern = "*_.*.tsv"))

print(samples)

cn= lapply(paste0(dirname, "/", 
                   list.files(dirname, pattern = "*_.*.tsv")), 
            function(f) read.table(f, header=T))
writeClusterTable(sciClone(cn, 
                           sampleNames = samples, 
                           maximumClusters = clusters,
                           #clusterMethod = "bmm"
                           clusterParams = "no.min.items.condition,no.apply.overlapping.std.dev.condition"
                           ),
                  paste0(dirname, "/", "clustering.tsv"))
