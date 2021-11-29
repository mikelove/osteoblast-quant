dir <- "boot_quants"
#dir <- "old_quants"
samples <- list.files(dir)
files <- file.path(dir, samples, "quant.sf")
cross <- factor(rep(c("129xB6","CASTxB6"),each=27))
day <- rep(rep(1:9 * 2, each=3), times=2)
rep <- rep(1:3, times=18)
names <- paste0(cross, "-", day, "-", rep)
coldata <- data.frame(cross, day, files, names)

#library(fishpond)
devtools::load_all("../fishpond/fishpond")
library(SummarizedExperiment)

library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()
#query(ah, c("EnsDb","102","Mus musculus"))
edb <- ah[["AH89211"]]

txps <- transcripts(edb, return.type="DataFrame")
tx2gene <- txps[,c("tx_id","gene_id")]

# attempt sub-gene resolution `Gnas` gene

isoform <- FALSE
if (isoform) {
  se <- importAllelicCounts(
    coldata, a1="alt", a2="ref", format="wide"
  )
  keep <- rowSums(assay(se) >= 10) >= 6
  table(keep)
  se <- se[keep,]
  all(rownames(se) %in% tx2gene[,1])
  library(org.Mm.eg.db)
  mcols(se)$gene <- tx2gene[match(rownames(se), tx2gene[,1]),2]
  mcols(se)$symbol <- mapIds(org.Mm.eg.db, mcols(se)$gene, "SYMBOL", "ENSEMBL")
  # collapse technical replicates
  idx <- 1:36 * 3
  se_coll <- se[,idx]
  nrep <- 30
  for (a in c("counts",paste0("infRep",1:nrep))) {
    cat(a,"")
    assay(se_coll,a) <- assay(se,a)[,idx] +
      assay(se,a)[,idx-1] + assay(se,a)[,idx-2]
  }
  for (a in c("abundance","length")) {
    cat(a,"")
    assay(se_coll,a) <- (assay(se,a)[,idx] +
                          assay(se,a)[,idx-1] + assay(se,a)[,idx-2])/3
  }
  colnames(se_coll) <- sub("^(.*-.*)-.*-(a.)$","\\1-\\2",colnames(se_coll))
  save(se_coll, file="data/se_filtered_collapsed.rda")
}

gse <- importAllelicCounts(
  coldata, a1="alt", a2="ref",
  format="wide", tx2gene=tx2gene
)
keep <- rowSums(assay(gse) >= 10) >= 6
table(keep)
gse <- gse[keep,]
library(org.Mm.eg.db)
mcols(gse)$symbol <- mapIds(org.Mm.eg.db, rownames(gse), "SYMBOL", "ENSEMBL")
save(gse, file="data/gse_filtered.rda")
#save(gse, file="old_data/gse_filtered.rda")

# collapse technical replicates
idx <- 1:36 * 3
gse_coll <- gse[,idx]
#nrep <- 30
nrep <- 20
for (a in c("counts",paste0("infRep",1:nrep))) {
  cat(a,"")
  assay(gse_coll,a) <- assay(gse,a)[,idx] +
    assay(gse,a)[,idx-1] + assay(gse,a)[,idx-2]
}
for (a in c("abundance","length")) {
  cat(a,"")
  assay(gse_coll,a) <- (assay(gse,a)[,idx] +
    assay(gse,a)[,idx-1] + assay(gse,a)[,idx-2])/3
}
colnames(gse_coll) <- sub("^(.*-.*)-.*-(a.)$","\\1-\\2",colnames(gse_coll))
#save(gse_coll, file="data/gse_filtered_collapsed.rda")
save(gse_coll, file="old_data/gse_filtered_collapsed.rda")

library(DESeq2)
gse_coll$fday <- factor(gse_coll$day)
gse_coll2 <- gse_coll
assays(gse_coll2) <- assays(gse_coll2)["counts"]
dds <- DESeqDataSet(gse_coll2, ~cross + cross:fday)
dds <- estimateSizeFactors(dds)
sf <- sizeFactors(dds)
plot(sf[1:18], sf[19:36]); abline(0,1)
sizeFactors(dds) <- rep((sf[1:18] + sf[19:36])/2, 2)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, "cross")
plotPCA(vsd, "day")
plotPCA(vsd, "allele")

dds_sub <- dds[!is.na(mcols(dds)$symbol),]
rownames(dds_sub) <- mcols(dds_sub)$symbol

library(ggplot2)
plotGene <- function(gene) {
  suppressMessages({
    dat <- plotCounts(dds_sub, gene, intgroup=c("cross", "day", "allele"), returnData=TRUE)
  })
  maxcnt <- max(dat$count)
  ggplot(dat, aes(day, count, color=allele, group=allele)) +
    geom_point() + stat_smooth(se=FALSE, method="loess", formula=y~x) +
    ggtitle(gene) + ylim(0,1.1*maxcnt) + ylab("scaled count") +
    facet_wrap(~cross)
}

plotGene("Runx2")
plotGene("Cped1")
plotGene("Sparc")
plotGene("Col1a1")
plotGene("Col1a2")
plotGene("Lars2")

dds129 <- dds[,dds$cross == "129xB6"]
ddsCast <- dds[,dds$cross == "CASTxB6"]

plotRatio <- function(dds, main="") {
  total <- counts(dds)[,10:18] + counts(dds)[,1:9]
  ratio <- (counts(dds)[,10:18] + 5) / (total + 10)
  plot(total[total > 100], ratio[total > 100],
       log="x", col=rgb(0,0,0,.05), cex=1, main=main,
       xlab="total", ylab="alt ratio (with pseudocount=5)")
  abline(h=0.5, col="red", lty=2)
}

plotRatio(dds129, main="129xB6")
plotRatio(ddsCast, main="CASTxB6")

histRatio <- function(dds, main="") {
  total <- counts(dds)[,10:18] + counts(dds)[,1:9]
  ratio <- (counts(dds)[,10:18] + 5) / (total + 10)
  #brks <- c(.45, seq(from=.45, to=.55, length=101) + .0005)
  #hist(ratio[total > 100 & ratio > .45 & ratio < .55], breaks=brks, xlab="ratio", main=main)
  hist(ratio[total > 100 & ratio > .2 & ratio < .8], breaks=1000, xlab="ratio", main=main)
}

histRatio(dds129, "129xB6")
histRatio(ddsCast, "CASTxB6")
