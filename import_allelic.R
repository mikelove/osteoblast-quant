samples <- list.files("quants")
files <- file.path("quants", samples, "quant.sf")
strain <- factor(rep(c("129xB6","CASTxB6"),each=27))
day <- rep(rep(1:9 * 2, each=3), times=2)
rep <- rep(1:3, times=54)
names <- paste0(strain, "-", day, "-", rep)
coldata <- data.frame(strain, day, files, names)

library(fishpond)
library(SummarizedExperiment)
se <- importAllelicCounts(coldata, a1="alt", a2="ref", format="wide")
# oops not enough memory...
