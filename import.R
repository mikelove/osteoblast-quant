samples <- list.files("ref_quants")
files <- file.path("ref_quants", samples, "quant.sf")
cross <- factor(rep(c("129xB6","CASTxB6"),each=27)) # cross
day <- rep(rep(1:9 * 2, each=3), times=2)
rep <- rep(1:3, times=18)
names <- paste0(cross, "-", day, "-", rep)
coldata <- data.frame(cross, day, files, names)

library(tximeta)
se <- tximeta(coldata)
gse <- summarizeToGene(se)

suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSet(gse, ~cross + day)
keep <- rowSums(counts(dds) >= 10) >= 6
table(keep)
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)

library(org.Mm.eg.db)
dds <- addIds(dds, "SYMBOL", gene=TRUE)

save(dds, file="ref_dds_filtered.rda")

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="cross")
plotPCA(vsd, intgroup="day")

dds_sub <- dds[!is.na(mcols(dds)$symbol),] # remove some genes
rownames(dds_sub) <- mcols(dds_sub)$symbol
plotCounts(dds_sub, "Runx2", intgroup="cross")

library(ggplot2)
plotGene <- function(gene) {
  suppressMessages({
    dat <- plotCounts(dds_sub, gene, intgroup=c("cross", "day"), returnData=TRUE)
  })
  maxcnt <- max(dat$count)
  ggplot(dat, aes(day, count, color=cross, group=cross)) +
    geom_point() + stat_smooth(se=FALSE, method="loess", formula=y~x) +
    ggtitle(gene) + ylim(0,1.1*maxcnt) + ylab("scaled count")
}

plotGene("Runx2")
plotGene("Cped1")
plotGene("Sparc")
plotGene("Col1a1")
plotGene("Col1a2")
plotGene("Lars2")

library(patchwork)
plotGene("Actb") + plotGene("B2m")
