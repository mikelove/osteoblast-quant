library(SummarizedExperiment)
devtools::load_all("../fishpond/fishpond")
load("data/gse_filtered_collapsed.rda")
y <- gse_coll
# the sample table
colData(y)
# one cross at a time
y <- y[,y$cross == "129xB6"]
