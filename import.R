samples <- list.files("ref_quants")
library(tximport)
library(readr)
files <- file.path("ref_quants", samples, "quant.sf")
strain <- factor(rep(c("129xB6","CASTxB6"),each=27))
day <- rep(rep(1:9 * 2, each=3), times=2)
rep <- rep(1:3, times=54)
names <- paste0(strain, "-", day, "-", rep)
coldata <- data.frame(strain, day, files, names)

library(tximeta)
se <- tximeta(coldata)
gse <- summarizeToGene(se)

suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSet(gse, ~strain + day)
keep <- rowSums(counts(dds) >= 10) >= 6
table(keep)
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="strain")
plotPCA(vsd, intgroup="day")

library(org.Mm.eg.db)
dds <- addIds(dds, "SYMBOL", gene=TRUE)

dds_sub <- dds[!is.na(mcols(dds)$symbol),]
rownames(dds_sub) <- mcols(dds_sub)$symbol
plotCounts(dds_sub, "Runx2", intgroup="strain")

library(ggplot2)
plotGene <- function(gene) {
  suppressMessages({
    dat <- plotCounts(dds_sub, gene, intgroup=c("strain", "day"), returnData=TRUE)
  })
  maxcnt <- max(dat$count)
  ggplot(dat, aes(day, count, color=strain, group=strain)) +
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
