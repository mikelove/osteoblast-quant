# Exploring counts script
# Michael Love
# Oct 6 2021

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
# [129xB6 B6 counts, CASTxB6 B6 counts, 129xB6 129 counts, CASTxB6 CAST counts]

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
ggplot(dat, aes(day, count, color=allele)) +
  geom_point() + geom_line() + facet_wrap(~cross) +
  ggtitle(gene) +
  scale_colour_brewer(palette = "Set1")

# again with total count (still not scaled for seq depth)
gene <- "Sparc"
ensgene <- getEnsID(gene)
dat <- as.data.frame( colData(tot) )
dat$count <- assay(tot, "counts")[ensgene,]
ggplot(dat, aes(day, count, color=cross)) +
  geom_point(size=2) + geom_line() + ggtitle(gene) +
  scale_colour_brewer(palette = "Dark2")

# Scaling for seq depth...
# ok, first, the above plots are not so bad, bc seq depth is fairly flat,
# still we probably want to scale for total expression comparisons
dat$seq_depth <- colSums(assay(tot))/1e6
ggplot(dat, aes(day, seq_depth, color=cross)) +
  geom_point(size=2) + geom_line() + ylim(0,45) + 
  scale_colour_brewer(palette = "Dark2") +
  ylab("seq depth (millions)")
  
# Scaling: I will recommend to use the methods of `tximport`
# for scaling such that isoform changes are taken into account when
# comparing counts across samples, each gene x sample has a "length".
# I can explain this more on a call, but it helps to prevent mistaking
# splicing changes as differential 'gene' expression
# (e.g. changes to the total output of gene)
