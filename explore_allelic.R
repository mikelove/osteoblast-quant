library(SummarizedExperiment)
devtools::load_all("../fishpond/fishpond")
load("data/gse_filtered.rda")
load("data/gse_filtered_collapsed.rda")

y <- gse_coll
y <- gse

# the sample table
colData(y)

# one cross at a time
y <- y[,y$cross == "129xB6"]

library(ggplot2)
gene <- "Cped1"
ensgene <- rownames(y)[which(mcols(y)$symbol == gene)]
getTrace(y, idx=ensgene, samp_idx=c("129xB6-14-a2","129xB6-14-a1")) %>%
  ggplot(aes(infRep, count, col=sample)) + geom_point() + geom_line() + ylim(0,5000)

getTrace(y, idx=8684, samp_idx=c("129xB6-14-a2","129xB6-14-a1")) %>% ggplot(aes(infRep, count, col=sample)) + geom_point() + geom_line()

# this is needed for plotting
y <- computeInfRV(y)

# no scaling! we look at ratios of alleles

# remove lowly expressed genes according to fishpond's function,
# note some filtering was already performed
y <- labelKeep(y)
table(mcols(y)$keep)
y <- y[mcols(y)$keep,]

# plot the inferential replicate data over time for one gene
gene <- "Cped1"
gene <- "Wdr26"
ensgene <- rownames(y)[which(mcols(y)$symbol == gene)]
ensgene

# NOTE: this use of 'x' and 'cov' is not what you will use for inference
plotInfReps(y, idx=ensgene, x="day", cov="allele",
            main=gene, legend=TRUE, legendPos="bottom", useMean=TRUE)
plotInfReps(y, idx=ensgene, x="day", cov="allele",
            main=gene, legend=TRUE, legendPos="bottom", useMean=FALSE)

