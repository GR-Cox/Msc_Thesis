projectName <- "George_Cox_Masters2"
codeDir <- "code"
projectDir <- "~/"
dataDir <- "data/run/processed_fastq/"
database <- "pr2_version_5.0.0_SSU_UTAX.fasta"
fPrimer <- "TTAAARVGYTCGTAGTYG"  
rPrimer <- "CCGTCAATTHCTTYAART"
minseqlength <- 270
maxee <- 3
maxdiffs <- 30                        ## Number of differences allowed when merging paired sequences. default is 5 but if using longer illumina reads with a long overlap region (eg: 616*F and 1132r for 18S euks) this needs to be increased. 
pct_id <- 65                          ## minimum percent identity match for alignments, default is 90 but should be decreased if working with long regions of overlap. 
stripL <- 0                           ## Bases to remove from start and end (24/26 are for Antonina's primers). Set both to 0 to not strip
stripR <- 0
coreN <- 16                          ## Number of cores to use for BLAST step
