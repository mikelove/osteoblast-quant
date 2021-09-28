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
mcols(y)$someInfo <- rowSums(abs(assay(y)[,y$allele == "a2"] -
                  assay(y)[,y$allele == "a1"]) < 1e-3) < n
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
