files <- list.files("orig_fastq",pattern=".gz")
# simplify file structure
new <- sub("(.*)_Col3.6-(.*)_GES.*L00(.)_(R.)_ALL.fastq.gz", "\\1_T\\2_L\\3_\\4.fastq.gz", files)
# simplify the Snakemake workflow:
new <- sub("L4","L1",new)
new <- sub("L5","L2",new)
new <- sub("L6","L3",new)
for (i in seq_along(files)) {
  cat(i)
  system(paste0("cp orig_fastq/", files[i], " fastq/", new[i]))
}
