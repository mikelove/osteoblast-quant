# Michael Love
# December 9, 2022

# Script to convert allelic (wide) data to total counts, including bootstraps

library(SummarizedExperiment)
load("data/gse_filtered.rda")
colData(gse)
cross <- "B6xCAST"
gse <- gse[,gse$cross == cross]
assaynms <- assayNames(gse)
for (a in assaynms) {
  assay(gse, a)[,1:9] <- assay(gse, a)[,1:9] + assay(gse, a)[,10:18]
}
gse <- gse[,1:9]
assay(gse, "length") <- 0.5 * assay(gse, "length") # average the lengths, not sum
gse$allele <- "total"
colSums(assay(gse, "abundance")) # should be close to 1e6
