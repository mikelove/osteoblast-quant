dir <- "quants"
samples <- list.files(dir)
files <- file.path(dir, samples, "quant.sf")
cross <- factor(rep(c("B6x129","B6xCAST"),each=9))
day <- rep(rep(1:9 * 2), times=2)
names <- paste0(cross, "-d", ifelse(day < 10, paste0("0", day), day))
coldata <- data.frame(cross, day, files, names)

#library(fishpond)
devtools::load_all("../../fishpond/fishpond")
library(SummarizedExperiment)

library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()
#query(ah, c("EnsDb","102","Mus musculus"))
edb <- ah[["AH89211"]]

# new fishpond code for importing allelic with GRanges
#txps <- transcripts(edb, return.type="DataFrame")
#tx2gene <- txps[,c("tx_id","gene_id")]

g <- genes(edb)

library(plyranges)
library(org.Mm.eg.db)

# increasing levels of aggregation
for (iter in 1:4) {
  # iter 1 = gene
  isoform <- FALSE
  tss <- FALSE
  fuzzy_tss <- FALSE
  # other types of output are here:
  if (iter == 2) {
    isoform <- TRUE; tss <- FALSE; fuzzy_tss <- FALSE
  }
  if (iter == 3) {
    tss <- TRUE; isoform <- FALSE; fuzzy_tss <- FALSE
  }
  if (iter == 4) {
    fuzzy_tss <- TRUE; isoform <- FALSE; tss <- FALSE
  }
  if (tss) {
    #tx2gene$gene_id <- paste0(tx2gene$gene_id, "-", txps$tx_seq_start)
    txps <- makeTx2Tss(edb) %>%
      select(tx_id, gene_id, group_id)
  } else if (fuzzy_tss) {
    txps <- makeTx2Tss(edb, maxgap=50)  %>%    
      select(tx_id, gene_id, group_id, tss)
  } else if (isoform) {
    txps <- transcripts(edb) %>%
      select(tx_id, gene_id)
  } else {
    # gene-level
    txps <- transcripts(edb) %>%
      select(tx_id, group_id=gene_id)
  }
  # import and save
  if (isoform) {
    se <- importAllelicCounts(
      coldata, a1="alt", a2="ref", format="wide"
    )
    rowRanges(se) <- txps[rownames(se)]
    keep <- rowSums(assay(se) >= 10) >= 6
    table(keep)
    se <- se[keep,]
    mcols(se)$symbol <- mapIds(org.Mm.eg.db, mcols(se)$gene_id,
                               "SYMBOL", "ENSEMBL")
    save(se, file="data/se_filtered.rda")
  } else {
    gse <- importAllelicCounts(
      coldata, a1="alt", a2="ref",
      format="wide", tx2gene=txps)
    keep <- rowSums(assay(gse) >= 10) >= 6
    table(keep)
    gse <- gse[keep,]
    if (tss | fuzzy_tss) {
      # tss-level
      mcols(gse)$symbol <- mapIds(org.Mm.eg.db, mcols(gse)$gene_id,
                                  "SYMBOL", "ENSEMBL")
      if (tss) {
        save(gse, file="data/tss_se_filtered.rda")
      } else {
        save(gse, file="data/fuzzy_50bp_tss_se_filtered.rda")
      }
    } else {
      # normal gene-level
      rowRanges(gse) <- g[rownames(gse)]
      # mcols(gse)$symbol <- mapIds(org.Mm.eg.db, rownames(gse), "SYMBOL", "ENSEMBL")
      save(gse, file="data/gse_filtered.rda")
    }
  }
}

###

if (FALSE) {

  library(DESeq2)
  gse$fday <- factor(gse$day)
  gse2 <- gse
  assays(gse2) <- assays(gse2)["counts"]
  dds <- DESeqDataSet(gse2, ~cross + cross:fday)
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

  dds129 <- dds[,dds$cross == "B6x129"]
  ddsCast <- dds[,dds$cross == "B6xCAST"]

  plotRatio <- function(dds, main="") {
    total <- counts(dds)[,10:18] + counts(dds)[,1:9]
    ratio <- (counts(dds)[,10:18] + 5) / (total + 10)
    plot(total[total > 100], ratio[total > 100],
         log="x", col=rgb(0,0,0,.05), cex=1, main=main,
         xlab="total", ylab="alt ratio (with pseudocount=5)")
    abline(h=0.5, col="red", lty=2)
  }

  plotRatio(dds129, main="B6x129")
  plotRatio(ddsCast, main="B6xCAST")

  histRatio <- function(dds, main="") {
    total <- counts(dds)[,10:18] + counts(dds)[,1:9]
    ratio <- (counts(dds)[,10:18] + 5) / (total + 10)
    #brks <- c(.45, seq(from=.45, to=.55, length=101) + .0005)
    #hist(ratio[total > 100 & ratio > .45 & ratio < .55], breaks=brks, xlab="ratio", main=main)
    hist(ratio[total > 100 & ratio > .2 & ratio < .8], breaks=1000, xlab="ratio", main=main)
  }

  histRatio(dds129, "B6x129")
  histRatio(ddsCast, "B6xCAST")

}
