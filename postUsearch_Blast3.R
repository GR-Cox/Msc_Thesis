source("~/code/Usearch_parameters3.R")
source("~/code/PipelineFunctions_abacusV3.txt")
dir <- paste(projectDir, dataDir, sep="")
paths <- setPaths()

OtusFasta <- readFAST(paste(dir, "otus.fa", sep=""))

uniquesfa <- readFAST(paste(dir, "uniques.fa", sep=""))
ZotusFasta <- readFAST(paste(dir, "zotus.fa", sep=""))


otuMatching <- read.table(paste(dir, "sequence_match_to_Otus.txt", sep=""))
zotuMatching <- read.table(paste(dir, "sequence_match_to_Zotus.txt", sep=""))


otus <- runBlastOnOtus(dir, OtusFasta, database, cores=coreN)

##THIS BIT BACTERIA SPECIFIC:
#fullHeaders <- readFAST(paste(paths$dbPath, "SILVA_138_SSURef_tax_silva.fasta", sep=""))
#fullHeaders$ID  <- unlist(sapply(strsplit(fullHeaders$descs, " "), "[", 1))
#fullHeaders$ID <- gsub("\\.", "", fullHeaders$ID)
#otus$fullHeader <- fullHeaders$descs[match(otus$otuMatch, fullHeaders$ID)]

save(file=paste(projectName, "OtuOutput.Rdat", sep="_"), list=c("otus", "OtusFasta", "uniquesfa", "ZotusFasta", "otuMatching", "zotuMatching"))