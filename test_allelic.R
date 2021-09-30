library(SummarizedExperiment)
devtools::load_all("../fishpond/fishpond")
load("data/gse_filtered_collapsed.rda")
y <- gse_coll

# the sample table
colData(y)

# one cross at a time
y <- y[,y$cross == "129xB6"]
y <- labelKeep(y)
table(mcols(y)$keep) # already filtered, so not many left
y <- y[mcols(y)$keep,]

# prefer symbols if we have them
symOrEns <- ifelse(is.na(mcols(y)$symbol), rownames(y), mcols(y)$symbol)
rownames(y) <- symOrEns

# assess the InfRV
y <- computeInfRV(y)
hist(log10(mcols(y)$meanInfRV + .001))

# remove genes with no information (this code specific to our data)
n <- ncol(y)/2
mcols(y)$someInfo <- rowSums(abs(assay(y, "infRep1")[,y$allele == "a2"] -
                  assay(y, "infRep1")[,y$allele == "a1"]) < 1) < n
table(mcols(y)$someInfo)
boxplot(log10(mcols(y)$meanInfRV + .001) ~ mcols(y)$someInfo)
y <- y[mcols(y)$someInfo,]

# no scaling, because we are comparing alleles within samples
y <- swish(y, x="allele", pair="day")

hist(mcols(y)$pvalue)

plotMASwish(y)

# worth looking at these
with(mcols(y), which.max(log2FC * as.numeric(qvalue < .05)))
with(mcols(y), which.min(log2FC * as.numeric(qvalue < .05)))

gene <- "Runx2"
plotInfReps(y, gene, x="day", cov="allele")
mcols(y)["Runx2",] # not global AI

# dynamic AI test
y <- swish(y, x="allele", pair="day", cov="day", cor="pearson")

hist(mcols(y)$stat)
plot(mcols(y)$stat, -log10(mcols(y)$pvalue))
hist(mcols(y)$pvalue)

with(mcols(y), plot(stat, log2FC, ylim=c(-1,1)))
abline(h=0, col="red")
with(mcols(y), identify(stat, log2FC))

# worth looking at these
with(mcols(y), which.max(log2FC * as.numeric(qvalue < .05))) # 3111
with(mcols(y), which.min(log2FC * as.numeric(qvalue < .05))) # 495

gene <- "Msx2"
mcols(y)[gene,]
plotInfReps(y, gene, x="day", cov="allele")
# symbol      keep meanInfRV  someInfo log10mean       stat    log2FC    pvalue    locfdr    qvalue
# <character> <logical> <numeric> <logical> <numeric>  <numeric> <numeric> <numeric> <numeric> <numeric>
#   Msx2        Msx2      TRUE   4.11872      TRUE   3.01567 0.00780183  0.603582  0.640518         1  0.653212

plotInfReps(y, 3111, x="day", cov="allele")
plotInfReps(y, 495, x="day", cov="allele")

dat <- data.frame(count=assay(y)["Msx2",], allele=y$allele, day=y$day)
library(ggplot2)
ggplot(dat, aes(day, count, color=allele, group=allele)) + geom_point() + geom_line() + ggtitle("Msx2")
ggplot(dat, aes(day, count, color=allele, group=allele)) + geom_point() + geom_line() + scale_y_log10() + ggtitle("Msx2")

# 467 1036 1560 2960

plotInfReps(y, 1036, x="day", cov="allele")
plotInfReps(y, 1560, x="day", cov="allele")

plotInfReps(y, 261, x="day", cov="allele")
plotInfReps(y, 1527, x="day", cov="allele")
plotInfReps(y, 2222, x="day", cov="allele")
plotInfReps(y, 2383, x="day", cov="allele")
plotInfReps(y, 3942, x="day", cov="allele")
