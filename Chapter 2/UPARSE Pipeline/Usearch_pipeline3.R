source("~/code/Usearch_parameters3.R")           ## Source all of the parameters for the run, note that these are used at different stages
source("~/code/PipelineFunctions_abacusV3.txt")  ## Huge file of functions that Ian has written
dir <- paste(projectDir, dataDir, sep="")             ## directory where your data is
paths <- setPaths()                                   ## This sets the paths to your various locations of softward
print(paste("Paths :", paths))
runVparse(dir, fPrimer, rPrimer, minseqlength = minseqlength, maxee = maxee, skipMerge = FALSE, stripL = stripL, stripR = stripR, maxdiffs = maxdiffs, pct_id = pct_id)
systemP(paste(paths$vsearch, " --fastq_stats ", dir, "/merged_reads.fq -log logStats.log", sep="", " --fastq_qmax 42"))
print("Usearch pipeline finished")
