files <- commandArgs(TRUE)
#files <- c(list.files("txps", "transcripts.fa", full.names=TRUE),
#           "txps/Mus_musculus.GRCm38.102.fa")

library(Biostrings)
cdna <- lapply(files, readDNAStringSet)
nstrains <- length(files) - 1
ref <- length(files)
for (i in 1:nstrains) {
  names(cdna[[i]]) <- sub(" .*","",names(cdna[[i]]))
}
# check equal:
stopifnot(all.equal(names(cdna[[1]]), names(cdna[[2]])))

# re-order and subset the reference (some missing from g2gtools output)
cdna[[ref]] <- cdna[[ref]][names(cdna[[1]])]

# check equal:
stopifnot(all.equal(names(cdna[[1]]), names(cdna[[ref]])))

# now add a ref and strain ('alt') indicator on the txp names
names(cdna[[ref]]) <- paste0(names(cdna[[ref]]), "_ref")
writeXStringSet(cdna[[ref]], filepath="txps/Mus_musculus.GRCm38.102.subset.fa")

for (i in 1:nstrains) {
  names(cdna[[i]]) <- paste0(names(cdna[[i]]), "_alt")
  out <- sub("transcripts", "tagged", files[i])
  writeXStringSet(cdna[[i]], filepath=out)
}

# check basepairs
## ws <- sapply(cdna, width)
## plot(ws[,1:2], log="xy")
## plot(ws[,c(1,3)], log="xy")
