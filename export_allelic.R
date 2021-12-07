library(SummarizedExperiment)
#devtools::load_all("../../fishpond/fishpond")
load("data/gse_filtered.rda")
y <- gse
assays(y) <- assays(y)[1:3] # drop inf reps
summary(colSums(assay(y))/1e6)
nsamp <- 18 # 129xB6 + CASTxB6 together
tot <- y[,1:nsamp]
# add the two alleles
idx <- (nsamp+1):(2*nsamp)
assay(tot, "counts") <- assay(tot, "counts") +
                        assay(y, "counts")[,idx]
summary(colSums(assay(tot))/1e6)
for (a in c("abundance","length")) {
  assay(tot, a) <- (assay(tot, a) + assay(y, a)[,idx])/2
}
tot$allele <- NULL
colnames(tot) <- sub("-a2","",colnames(tot))
save(y, file="data/allelic_counts.rda")
save(tot, file="data/total_counts.rda")
