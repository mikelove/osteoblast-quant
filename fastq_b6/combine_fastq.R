x <- read.csv("sra_result.csv")
day <- sub("GSM.*B6_(.*)_rep.*","\\1",x[,2])
dat <- data.frame(day, Experiment=x[,1])
y <- read.csv("SraRunTable.txt")
y <- y[,c("Run","Experiment")]
dat <- merge(dat, y)
days <- 1:9 * 2
for (d in days) {
  print(d)
  system(paste0("cat ", paste(paste0(dat$Run[dat$day == d],".fastq"), collapse=" "), " > B6_T", d/2, ".fastq"))
}
