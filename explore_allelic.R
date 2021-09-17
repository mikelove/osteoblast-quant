library(SummarizedExperiment)
library(fishpond)
load("data/gse_filtered_collapsed.rda")

# the sample table
colData(gse_coll)

# one cross at a time
y <- gse_coll[,gse_coll$cross == "129xB6"]

# this is needed for plotting
y <- computeInfRV(y)

# no scaling! we look at ratios of alleles

# remove lowly expressed genes according to fishpond's function,
# note some filtering was already performed
y <- labelKeep(y)
table(mcols(y)$keep)
y <- y[mcols(y)$keep,]

# plot the inferential replicate data over time for one gene
gene <- "Runx2"
ensgene <- rownames(y)[which(mcols(y)$symbol == gene)]

# NOTE: this use of 'x' and 'cov' is not what you will use for inference
plotInfReps(y, idx=ensgene, x="day", cov="allele", main=gene)
