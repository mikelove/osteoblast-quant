library(SummarizedExperiment)
devtools::load_all("../fishpond/fishpond")
load("data/gse_filtered_collapsed.rda")
y <- gse_coll
assays(y) <- assays(y)[1:3] # drop inf reps
summary(colSums(assay(y))/1e6)
tot <- y[,1:18]
# add the two alleles
assay(tot, "counts") <- assay(tot, "counts") + assay(y, "counts")[,19:36]
for (a in c("abundance","length")) {
  assay(tot, a) <- (assay(tot, a) + assay(y, a)[,19:36])/2
}
tot$allele <- NULL
colnames(tot) <- sub("-a2","",colnames(tot))
save(y, file="data/allelic_counts.rda")
save(tot, file="data/total_counts.rda")
