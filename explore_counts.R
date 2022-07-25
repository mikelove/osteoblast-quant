# Exploring counts script
# Michael Love
# Dec 6 2021

# if you don't have SummarizedExperiment, install it with these 3 lines:
if (FALSE) {
  if (!requireNamespace("BiocManager"), quietly=TRUE)
    install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
}

# load package and the allelic + total counts
library(SummarizedExperiment)
load("data/allelic_counts.rda")
load("data/total_counts.rda")

# FYI none of the counts are scaled yet, more notes on that below...

# Allelic counts: (reference a2 = B6 and alternate a1 = 129 or CAST)

# In the data, the B6 allelic count is first, then the other strain is next
# e.g. it looks like the following matrices stuck together (9 samples each)
# [B6x129 B6 counts, B6xCAST B6 counts, B6x129 129 counts, B6xCAST CAST counts]

# the metadata for the allelic data
colData(y)

# y = allelic counts
# tot = total counts

# the metadata for the total count data:
colData(tot)

# the counts are here:
assay(y, "counts")[1:3, 1:9]

# also abundance (TPM) and length (bp)
assay(y, "abundance")[1:3,1:3]
assay(y, "length")[1:3,1:3]

# example: plot Runx2 alleles
library(ggplot2)

# match the gene symbol to Ensembl ID (our row identifiers)
getEnsID <- function(symbol) {
  stopifnot(symbol %in% mcols(y)$symbol)
  ens <- rownames(y)[ which(mcols(y)$symbol == symbol) ]
  if (length(ens) > 1) warning("more than one symbol")
  ens
}

gene <- "Runx2"
ensgene <- getEnsID(gene)

# convert metadata from Bioc DataFrame to normal data.frame
dat <- as.data.frame( colData(y) )

# add the count for this gene to our little table 
dat$count <- assay(y, "counts")[ensgene,]

# remember, a2 = B6 in all plots
cols <- c(a2="dodgerblue", a1="goldenrod3")
ggplot(dat, aes(day, count, color=allele)) +
  geom_point() + geom_line() + facet_wrap(~cross) +
  ggtitle(gene) +
  scale_colour_manual(values=cols) +
  guides(color = guide_legend(reverse=TRUE))

# again with total count (still not scaled for seq depth)
gene <- "Sparc"
ensgene <- getEnsID(gene)
dat <- as.data.frame( colData(tot) )
dat$count <- assay(tot, "counts")[ensgene,]
ggplot(dat, aes(day, count, color=cross)) +
  geom_point(size=2) + geom_line() + ggtitle(gene) +
  scale_colour_manual(values=cols) +
  guides(color = guide_legend(reverse=TRUE))

# Scaling for seq depth...
# ok, first, the above plots are not so bad, bc seq depth is fairly flat,
# still we probably want to scale for total expression comparisons
dat$seq_depth <- colSums(assay(tot))/1e6
ggplot(dat, aes(day, seq_depth, color=cross)) +
  geom_point(size=2) + geom_line() + ylim(0,50) +
  scale_colour_manual(values=cols) +
  guides(color = guide_legend(reverse=TRUE))
  ylab("seq depth (millions)")
  
# Scaling: I will recommend to use the methods of `tximport`
# for scaling such that isoform changes are taken into account when
# comparing counts across samples, each gene x sample has a "length".
# I can explain this more on a call, but it helps to prevent mistaking
# splicing changes as differential 'gene' expression
# (e.g. changes to the total output of gene)

# use DESeq2 scaling
library(DESeq2)
tot$day.scale <- (tot$day - 10) / 8
dds <- DESeqDataSet(tot, ~cross + day.scale)

# correct for sequencing depth
dds <- estimateSizeFactors(dds)

# the corrections are here
normalizationFactors(dds)[1:3,1:5]

# change rownames to gene symbols, if we have one
rownames(dds) <- ifelse(is.na(mcols(dds)$symbol),
                        rownames(dds),
                        mcols(dds)$symbol)

# now plotting scaled counts:
gene <- "Runx2"
dat <- plotCounts(dds, gene, intgroup=c("cross","day"), returnData=TRUE)
ggplot(dat, aes(day, count, col=cross)) +
  geom_point(size=2) + geom_line() +
  scale_colour_manual(values=cols) +
  guides(color = guide_legend(reverse=TRUE))
  ylim(0, 1.1 * max(dat$count)) + ggtitle(gene)

# scaled counts are here
counts(dds, normalized=TRUE)["Runx2",1:5]

# log transformed scaled counts for heatmap / sample distance / PCA
vsd <- vst(dds, blind=FALSE)
library(pheatmap)
pheatmap(assay(vsd)[c("Runx2","Sparc","Col1a1"),],
         scale="row", cluster_cols=FALSE)

# log transformed scaled counts are here
assay(vsd)["Runx2",1:5]

# PCA plots
plotPCA(vsd, intgroup="day")
plotPCA(vsd, intgroup="cross") +
  scale_colour_manual(values=cols) +
  guides(color = guide_legend(reverse=TRUE))

############################################################

# alternatively, forget about gene lengths, and just
# do simple scaling per sample.
# Note: I don't recommend this, I think the above is a better approach
dds2 <- dds
assays(dds2) <- assays(dds2)[1] # throw out abundance and length info
dds2 <- estimateSizeFactors(dds2)
sizeFactors(dds2) # the size factors
# size factors here capturing the deviation from middle total count
plot(colSums(counts(dds2))/1e6,
     sizeFactors(dds2), xlim=c(0,50), ylim=c(0,2))
abline(h=1, v=median(colSums(counts(dds2))/1e6,), lty=2)


#############################################################

# assess the overdispersion for the simulation

library(DESeq2)
y$sample <- factor(paste0(y$cross, "_", y$day))
dds <- DESeqDataSet(y, ~sample + allele)
keep <- rowSums(counts(dds) >= 10) == ncol(dds)
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
nrow(dds)
system.time({
  dds <- estimateDispersions(dds)
})

library(dplyr)
library(ggplot2)
dat <- mcols(dds) %>%
  as.data.frame() %>%
  select(baseMean, dispGeneEst)
dat %>% filter(baseMean > 1e3, dispGeneEst > 1e-4) %>%
  ggplot(aes(dispGeneEst)) +
  geom_histogram(color="red",fill="mistyrose2") +
  scale_x_log10()
