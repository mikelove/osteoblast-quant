samples <- list.files("ref_quants_b6")
files <- file.path("ref_quants_b6", samples, "quant.sf")
day <- 1:9 * 2
names <- paste0("B6-", day)
coldata <- data.frame(day, files, names)

library(tximeta)
se <- tximeta(coldata)
gse <- summarizeToGene(se)

library(org.Mm.eg.db)
se <- addIds(se, "SYMBOL", gene=TRUE)
gse <- addIds(gse, "SYMBOL", gene=TRUE)

library(SummarizedExperiment)
keep <- rowSums(assay(se) >= 10) >= 6
table(keep)
se <- se[keep,]

keep <- rowSums(assay(gse) >= 10) >= 6
table(keep)
gse <- gse[keep,]

save(se, file="data/b6_se_filtered.rda")
save(gse, file="data/b6_gse_filtered.rda")
