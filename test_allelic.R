library(SummarizedExperiment)
devtools::load_all("../../fishpond/fishpond")
#load("data/gse_filtered.rda")
#load("data/se_filtered.rda")
#load("data/tss_se_filtered.rda")
load("data/fuzzy_50bp_tss_se_filtered.rda")
y <- gse

# the sample table
colData(y)

# one cross at a time
y <- y[,y$cross == "CASTxB6"]
y <- labelKeep(y)
table(mcols(y)$keep) # already filtered, so not many left
y <- y[mcols(y)$keep,]

# prefer symbols if we have them
gene <- FALSE
tss <- FALSE
if (gene) {
  mcols(y)$gene <- rownames(y)
  symOrEns <- ifelse(is.na(mcols(y)$symbol), rownames(y), mcols(y)$symbol)
  rownames(y) <- symOrEns
} else if (tss) {
  symOrEns <- ifelse(is.na(mcols(y)$symbol), mcols(y)$gene, mcols(y)$symbol)
  tss <- sapply(strsplit(rownames(y), "-"), `[`, 2)
  rownames(y) <- paste0(symOrEns,"-",tss)
} else {
  symOrEns <- ifelse(is.na(mcols(y)$symbol), mcols(y)$gene, mcols(y)$symbol)
  rownames(y) <- paste0(symOrEns,"-",rownames(y))
}

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
y <- swish(y, x="allele", cov="day", pair="day", cor="pearson")

hist(mcols(y)$pvalue)

# worth looking at these
with(mcols(y), which.max(stat))
with(mcols(y), which.min(stat))

plotInfReps(y, 5953, x="day", cov="allele", shiftX=.15)
plotInfReps(y, 16127, x="day", cov="allele", shiftX=.15)
plotInfReps(y, 32850, x="day", cov="allele", shiftX=.15)

library(tibble)
tss <- mcols(y) %>% as.data.frame() %>% rownames_to_column("id") %>% tibble()
save(tss, file="tss.rda")
gene <- mcols(y) %>% as.data.frame() %>% rownames_to_column("id") %>% tibble()
save(gene, file="gene.rda")

# try dynamic AI for up/down pattern
plot(y$day[1:9])
y$day.sq <- (y$day - 10)^2
plot(y$day.sq[1:9])
y$day.sin <- sin((y$day - 2) / 16 * (2 * pi))
plot(y$day.sin[1:9])

y <- swish(y, x="allele", pair="day", cov="day.sq", cor="pearson")
y <- swish(y, x="allele", pair="day", cov="day.sin", cor="pearson")

hist(mcols(y)$stat)
plot(mcols(y)$stat, -log10(mcols(y)$pvalue))
hist(mcols(y)$pvalue)

with(mcols(y), plot(stat, log2FC, ylim=c(-1,1)))
abline(h=0, col="red")
#with(mcols(y), identify(stat, log2FC))

# worth looking at these
with(mcols(y), which.max(stat))
with(mcols(y), which.min(stat))

head(mcols(y)[order(mcols(y)$pvalue),],10)

plotInfReps(y, "Stc1", x="day", cov="allele")

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
