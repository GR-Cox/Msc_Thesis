#DADA2 Pipeline with variable parameters
#each parameter needs to be defined at the top for easy editing, including pile paths and names.

#Parameters Block (CHECK DATES AND FILE PATHS)
#_______________________________________________________________________________
Parameter_Names <- c("Date", "datapath", "outpath", "taxpath", "FWD Primer", "FWD Primer Length", "REV Primer",
                     "REV Primer Length", "Amplicon Length", "Amplicon Length Minus Primers", "Contig Pipeline",
                     "Chop_F_contig", "Chop_R_contig", "EE_F_contig", "EE_R_contig", "qual_chop_contig", "Pool_contig", 
                     "Forward Only Pipeline", "Chop_F_forward_only", "EE_F_forward_only", "qual_chop_forward_only", "Pool_forward_only", 
                     "N Mash Pipeline", "Chop_F_N_mash", "Chop_R_N_mash", "EE_F_N_mash", "EE_R_N_mash", "qual_chop_N_mash", "Pool_N_mash", "Proceed to Formatting")
Date <- "20_9_23"     #Date when run will start
datapath <- "~/euk_data/processed_fastq" #Shows where the demultiplexed sequences are. 
dir.create('~/outputs/20_9_23')  #Directory where all outputs will go
outpath <- '~/outputs/20_9_23' #Where all outputs will end up
taxpath <- "~/euk_data/pr2_version_5.0.0_SSU_dada2.fasta.gz"   #Where the assign taxonomy function can find the sequence database, in this case the PR2 database.

FWD <- "TTAAARVGYTCGTAGTYG"     #Forward primer sequence, in this case its 616*F
FWD_length <- nchar(FWD) #Returns the length of the forward primer sequence, for 616*F it should be 18 bp
REV <- "CCGTCAATTHCTTYAART"     #Reverse primer sequence, in this case its 1132r
REV_length <- nchar(REV) #returns the length of the reverse primer sequence, for 1132r it should be 18bp
Amplicon_length <- 535 #Amplicon length including primer seqeunces, for the 616*F-1132r primer pair it should be 535bp
Amplicon_length_minus_primers <- (Amplicon_length-(FWD_length+REV_length)) #Records the expected length of the amplicon once primers have been removed, nessecary for the N_mash_pipeline and useful to know for optimising trimming parameters.

Contig_pipeline <- FALSE #Runs the dada2 pipeline as intended which merges both the forward and reverse reads to make contigs, requires less stringent trimming parameters or too many reads will be lost.
Chop_F_contig <- 270  #Sets the Trunclen parameter for the forward reads for contiging
Chop_R_contig <- 250  #Sets the Trunclen parameter for the reverse reads for contiging
EE_F_contig <- 3    #maxEE parameter 1 (Forward Reads) for contiging
EE_R_contig <- 6     #maxEE parameter 2 (Reverse Reads) for contiging
qual_chop_contig <- 1  #min quality score for contiging (truncQ)
pool_contig <- FALSE #Whether or not the dada2 denoising step should pool the samples (takes longer but can detect more rare sequence variants), if you want to pool change pool=Psuedo (intermediate) or pool=true (takes ages but super high sensitivity)

Forward_only_pipeline <- FALSE #Runs the dada2 pipeline without merging the reads, uses only the forward reads allowing for more stringent trimming parameters (Idea from Nick Dragone)
Chop_F_forward_only <- 270 #Sets the TruncLen parameter for the forward only pipeline.
EE_F_forward_only <- 3 #Sets the maxEE value for the forward only pipeline.
qual_chop_forward_only <- 1 #set the truncQ parameter for the forwarc only pipeline.
pool_forward_only <- FALSE

N_mash_pipeline <- TRUE #Runs the dada2 pipeline using both forward and reverse reads but instead of contiging, the reads are trimmed heavily without need for overlap, then they are merged togther using ambigious base pairs (N) in the middle. Allows for seqeuncing errors to be removed at the cost of losing some taxonomic information. (Idea from Craig Herbold)
Chop_F_N_mash <- 250  #Sets the Trunclen parameter for the forward reads for the N Mash pipeline.
Chop_R_N_mash <- 200  #Sets the Trunclen parameter for the reverse reads for the N Mash pipeline.
EE_F_N_mash <- 2     #maxEE parameter 1 (Forward Reads) for the N Mash pipeline.
EE_R_N_mash <- 3      #maxEE parameter 2 (Reverse Reads) for the N Mash pipeline.
qual_chop_N_mash <- 2  #min quality score for the N Mash pipeline. (truncQ)
pool_N_mash <- FALSE

proceed_to_formatting <- FALSE #Whether or not the pipeline should immediately run the formatting R script to generate fasta files and taxonomically filtered data.

Parameter_Values <- c(Date, datapath, outpath, taxpath, FWD, FWD_length, REV, REV_length, Amplicon_length, Amplicon_length_minus_primers, Contig_pipeline, Chop_F_contig, Chop_R_contig, EE_F_contig, EE_R_contig, qual_chop_contig, pool_contig, Forward_only_pipeline, Chop_F_forward_only, EE_F_forward_only, qual_chop_forward_only, pool_forward_only, N_mash_pipeline, Chop_F_N_mash, Chop_R_N_mash, EE_F_N_mash, EE_R_N_mash, qual_chop_N_mash, pool_N_mash, proceed_to_formatting) #used to make the parameter matrix which gives info on each run

#Printing Parameters for stdout
Parameters <- data.frame(Parameter_Names, Parameter_Values)
print(Parameters)


#Packages Block
#_______________________________________________________________________________

#Package Installation, Should only need to be done once per environment.
#this pipeline also requires cutadapt to be loaded in the environment (in my case the dada2env)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("dada2")
#BiocManager::install("Shortread")
#install.packages("stringr")

library(stringr) 
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
cutadapt <- "cutadapt" #finds cutadapt on the VM, but only seems to work if cutadapt is sitting in the home directory (/Home/rccuser/cutadapt).
system2(cutadapt, args = "--version") #Tests to make sure cutadapt has actually loaded.


#Functions Block
#_______________________________________________________________________________
set.seed(100) #Sets seed so that randomised steps are replicatable

allOrients <- function(primer) {
  # Create all orientations of the input sequence eg: primers
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


getN <- function(x) sum(getUniques(x)) #Gets the Number of a thing, used to track reads.



#Code to run on VM
#_________________________________________________________________________________________
#conda activate dada2env
#nohup Rscript ~/pipelines/DADA2_pipeline_DATE.R >~/outputs/stdout_pipeline_DATE.txt
#_________________________________________________________________________________________



#Pipeline Start
#________________________________________________________________________________________


list.files(datapath) #checks that the pipeline can find the processed sequences

#Cleansing the rundir of previous run things, if important will have already been moved and renamed at end of previous run
unlink('~/rundir', recursive = TRUE)

dir.create('~/rundir') #Recreates a new clean rundir to dump all files from this run into. Will only get deleted when a new run occurs
setwd("~/rundir")
rundir <- "~/rundir"

#Creating directories within the rundir for filtering.
dir.create('~/rundir/N_filtered')     #No N's in this directory
dir.create('~/rundir/P_filtered')     #No primers in this directory
dir.create('~/rundir/filtered_contig')       #Post trimming, ready for denoising for the normal contiging pipeline.
dir.create('~/rundir/filtered_forward_only') #Post trimming, ready for denoising for the forward only pipeline.
dir.create('~/rundir/filtered_N_mash') #Post trimming, ready for denoising for the N mash pipeline.



# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(datapath, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(datapath, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
#changed the code to match my file names as it just outputed the first word when file paths and _ were removed. First word was processed but due to changing the number on the last agrument it is now correct (2nd word)
head(sample_names)

pdf(file = '~/rundir/F_quality_plot_DATE.pdf', width = 10, height = 10)
plotQualityProfile(fnFs)
dev.off()

pdf(file = '~/rundir/R_quality_plot_DATE.pdf', width = 10, height = 10)
plotQualityProfile(fnRs)
dev.off()

#Removes any sequences with Ns as dada2 and cutadapt can't deal with ambigious base pairs. 
filtF_Ns <- file.path('~/rundir/N_filtered', paste0(sample_names, "_F_filt_N.fastq.gz"))
filtR_Ns <- file.path('~/rundir/N_filtered', paste0(sample_names, "_R_filt_N.fastq.gz"))
names(filtF_Ns) <- sample_names
names(filtR_Ns) <- sample_names
N_out <- filterAndTrim(fnFs, filtF_Ns, fnRs, filtR_Ns, maxN = 0, multithread = TRUE) 
head(N_out)

#Checking for primers in the sequences, primers must have been defined in the parameters section


FWD_orients <- allOrients(FWD)
REV_orients <- allOrients(REV)
FWD_orients

rbind(FWD.ForwardReads = sapply(FWD_orients, primerHits, fn = filtF_Ns[[1]]), 
      FWD.ReverseReads = sapply(FWD_orients, primerHits, fn = filtR_Ns[[1]]), 
      REV.ForwardReads = sapply(REV_orients, primerHits, fn = filtF_Ns[[1]]), 
      REV.ReverseReads = sapply(REV_orients, primerHits, fn = filtR_Ns[[1]]))

filtF_Ps <- file.path('~/rundir/P_filtered', paste0(sample_names, "_F_filt_P.fastq.gz")) #Sets up files for primer removal with cutadapt
filtR_Ps <- file.path('~/rundir/P_filtered', paste0(sample_names, "_R_filt_P.fastq.gz"))
names(filtF_Ps) <- sample_names
names(filtR_Ps) <- sample_names

FWD_RC <- dada2:::rc(FWD) #saves the reverse complements of the primers
REV_RC <- dada2:::rc(REV)

#Cutadapt flags which tell cutadapt to trim the primers and their complements from the reads

R1_flags <- paste("-g", FWD, "-a", REV_RC) #Parameters taken from Fierer lab github
R2_flags <- paste("-G", REV, "-A", FWD_RC) 

#Run cutadapt

for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1_flags, R2_flags, "-n", 2, # -n 2 required to remove FWD and REV primers from reads
                             "-o", filtF_Ps[i], "-p", filtR_Ps[i], # output files
                             filtF_Ns[i], filtR_Ns[i], "--quiet", "-q", 0)) # input files + quiet command as this bit of code will make the stdout very unwieldy. 
  #-q 0 command should stop cutadapt from doing any quality filtering as other filters are being employed to do that.
}

#Checking that primers have been removed, should all be zero
rbind(FWD.ForwardReads = sapply(FWD_orients, primerHits, fn = filtF_Ps[[1]]), 
      FWD.ReverseReads = sapply(FWD_orients, primerHits, fn = filtR_Ps[[1]]), 
      REV.ForwardReads = sapply(REV_orients, primerHits, fn = filtF_Ps[[1]]), 
      REV.ReverseReads = sapply(REV_orients, primerHits, fn = filtR_Ps[[1]]))

#End of pre trimming and start of pipeline proper

#Contig pipeline (normal dada2 pipeline)
#_______________________________________________________________________________
if (Contig_pipeline == TRUE) { #Turns on this pipeline if Contig_pipeline is set to TRUE
  
#Trimming

  filtFs_contig <- file.path('~/rundir/filtered_contig', paste0(sample_names, "_F_filt_contig.fastq.gz"))
  filtRs_contig <- file.path('~/rundir/filtered_contig', paste0(sample_names, "_R_filt_contig.fastq.gz"))
  names(filtFs_contig) <- sample_names
  names(filtRs_contig) <- sample_names

  trim_out_contig <- filterAndTrim(filtF_Ps, filtFs_contig, filtR_Ps, filtRs_contig, truncLen=c(Chop_F_contig,Chop_R_contig),
                     maxEE=c(EE_F_contig,EE_R_contig), truncQ=qual_chop_contig, maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

  print(head(trim_out_contig))

#Learning error rates

  errF_contig <- learnErrors(filtFs_contig, nbases = 1e8, multithread=TRUE, randomize = TRUE)
  errR_contig <- learnErrors(filtRs_contig, nbases = 1e8, multithread=TRUE, randomize = TRUE)

  pdf(file = '~/rundir/plotErrors_F_contig_DATE.pdf')
  print(plotErrors(errF_contig, nominalQ=TRUE))
  dev.off()

  pdf(file = '~/rundir/plotErrors_R_contig_DATE.pdf')
  print(plotErrors(errF_contig, nominalQ=TRUE))
  dev.off()

  dadaFs_contig <- dada(filtFs_contig, err=errF_contig, multithread=TRUE, pool = pool_contig)
  dadaRs_contig <- dada(filtRs_contig, err=errR_contig, multithread=TRUE, pool = pool_contig)
  dadaFs_contig[[1]]
  dadaRs_contig[[1]]

  mergers <- mergePairs(dadaFs_contig, filtFs_contig, dadaRs_contig, filtRs_contig, verbose=TRUE)
# Inspect the merger data.frame from the first sample
  print(head(mergers[[1]]))

  seqtab_contig <- makeSequenceTable(mergers)
  dim(seqtab_contig)
  table(nchar(getSequences(seqtab_contig)))

  saveRDS(seqtab_contig, file = '~/rundir/seqtab_contig_DATE.rds')


  seqtab_nochim_contig <- removeBimeraDenovo(seqtab_contig, method="consensus", multithread=TRUE, verbose=TRUE)
  dim(seqtab_nochim_contig)

  saveRDS(seqtab_nochim_contig, file = '~/rundir/seqtab_nochim_contig_DATE.rds')

  print(sum(seqtab_nochim_contig)/sum(seqtab_contig))

  track_contig <- cbind(N_out, trim_out_contig, sapply(dadaFs_contig, getN), sapply(dadaRs_contig, getN), sapply(mergers, getN), rowSums(seqtab_nochim_contig)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track_contig) <- c("input", "N_filtered", "P_filtered", "trimmed", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track_contig) <- sample_names
  print(head(track_contig))

  saveRDS(track_contig, file = '~/rundir/track_contig_DATE.rds')

  taxa_contig <- assignTaxonomy(seqtab_nochim_contig, taxpath, taxLevels = c("Domain","Supergroup","Division","Kingdom","Phylum","Order","Family","Genus","Species"), multithread=TRUE)

  saveRDS(taxa_contig, file = '~/rundir/taxa_contig_DATE.rds')

  taxa_print <- taxa_contig # Removing sequence rownames for display only
  rownames(taxa_print) <- NULL
  print(head(taxa_print))
  

  print("Contig Pipeline Finished")
} 
#End of the contig pipeline.



#Forward only pipeline
#_______________________________________________________________________________

#if
if (Forward_only_pipeline == TRUE) {#Turns on the forward read only pipeline if the object is set to TRUE in the parameters section.

#Trimming
  filtFs_forward_only <- file.path('~/rundir/filtered_forward_only', paste0(sample_names, "_F_filt_forward_only.fastq.gz"))
  names(filtFs_forward_only) <- sample_names
  
  
  trim_out_forward_only <- filterAndTrim(filtF_Ps, filtFs_forward_only, truncLen=(Chop_F_forward_only),
                                   maxEE=(EE_F_forward_only), truncQ=qual_chop_forward_only, maxN=0, rm.phix=TRUE,
                                   compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
  
  print(head(trim_out_forward_only))
  
#Learn errors

  errF_forward_only <- learnErrors(filtFs_forward_only, nbases = 1e8, multithread=TRUE, randomize = TRUE)
  
  pdf(file = '~/rundir/plotErrors_F_forward_only_DATE.pdf')
  print(plotErrors(errF_forward_only, nominalQ=TRUE))
  dev.off()
  

#DADA

  dadaFs_forward_only <- dada(filtFs_forward_only, err=errF_forward_only, multithread=TRUE, pool = pool_forward_only)
 
  dadaFs_forward_only[[1]]

  
#Seqtab

  seqtab_forward_only <- makeSequenceTable(dadaFs_forward_only)
  dim(seqtab_forward_only)
  table(nchar(getSequences(seqtab_forward_only)))
  
  saveRDS(seqtab_forward_only, file = '~/rundir/seqtab_forward_only_DATE.rds')
  
#Seqtab no chim  
  seqtab_nochim_forward_only <- removeBimeraDenovo(seqtab_forward_only, method="consensus", multithread=TRUE, verbose=TRUE)
  dim(seqtab_nochim_forward_only)
  
  saveRDS(seqtab_nochim_forward_only, file = '~/rundir/seqtab_nochim_forward_only_DATE.rds')
  
  print(sum(seqtab_nochim_forward_only)/sum(seqtab_forward_only))
  
#Tracking reads

  track_forward_only <- cbind(N_out, trim_out_forward_only, sapply(dadaFs_forward_only, getN), rowSums(seqtab_nochim_forward_only)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track_forward_only) <- c("input", "N_filtered", "P_filtered", "trimmed", "denoisedF","nonchim")
  rownames(track_forward_only) <- sample_names
  print(head(track_forward_only))
  
  saveRDS(track_forward_only, file = '~/rundir/track_forward_only_DATE.rds')
 
#Assigning Taxonomy
  
  taxa_forward_only <- assignTaxonomy(seqtab_nochim_forward_only, taxpath, taxLevels = c("Domain","Supergroup","Division","Kingdom","Phylum","Order","Family","Genus","Species"), multithread=TRUE)
  
  saveRDS(taxa_forward_only, file = '~/rundir/taxa_forward_only_DATE.rds')
  
  taxa_print <- taxa_forward_only# Removing sequence rownames for display only
  rownames(taxa_print) <- NULL
  print(head(taxa_print))
  

  print("Forward Only Pipeline Finished")
} 

#N Mash pipeline
#_______________________________________________________________________________
if (N_mash_pipeline == TRUE) { #Turns on this pipeline if Contig_pipeline is set to TRUE
  
  #Trimming
  
  filtFs_N_mash <- file.path('~/rundir/filtered_N_mash', paste0(sample_names, "_F_filt_N_mash.fastq.gz"))
  filtRs_N_mash <- file.path('~/rundir/filtered_N_mash', paste0(sample_names, "_R_filt_N_mash.fastq.gz"))
  names(filtFs_N_mash) <- sample_names
  names(filtRs_N_mash) <- sample_names
  
  trim_out_N_mash <- filterAndTrim(filtF_Ps, filtFs_N_mash, filtR_Ps, filtRs_N_mash, truncLen=c(Chop_F_N_mash,Chop_R_N_mash),
                                   maxEE=c(EE_F_N_mash,EE_R_N_mash), truncQ=qual_chop_N_mash, maxN=0, rm.phix=TRUE,
                                   compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
  
  print(head(trim_out_N_mash))
  
  #Learning error rates
  
  errF_N_mash <- learnErrors(filtFs_N_mash, nbases = 1e8, multithread=TRUE, randomize = TRUE)
  errR_N_mash <- learnErrors(filtRs_N_mash, nbases = 1e8, multithread=TRUE, randomize = TRUE)
  
  pdf(file = '~/rundir/plotErrors_F_N_mash_DATE.pdf')
  print(plotErrors(errF_N_mash, nominalQ=TRUE))
  dev.off()
  
  pdf(file = '~/rundir/plotErrors_R_N_mash_DATE.pdf')
  print(plotErrors(errF_N_mash, nominalQ=TRUE))
  dev.off()
  
  dadaFs_N_mash <- dada(filtFs_N_mash, err=errF_N_mash, multithread=TRUE, pool = pool_N_mash)
  dadaRs_N_mash <- dada(filtRs_N_mash, err=errR_N_mash, multithread=TRUE, pool = pool_N_mash)
  dadaFs_N_mash[[1]]
  dadaRs_N_mash[[1]]
  
  mergers_N_mash <- mergePairs(dadaFs_N_mash, filtFs_N_mash, dadaRs_N_mash, filtRs_N_mash, verbose=TRUE, justConcatenate = TRUE) #justConcatenate will mash the two sequences togther with a string of 10 Ns in the middle.
  # Inspect the merger data.frame from the first sample
  print(head(mergers_N_mash[[1]]))
  
  seqtab_N_mash <- makeSequenceTable(mergers_N_mash)
  dim(seqtab_N_mash)
  table(nchar(getSequences(seqtab_N_mash)))
  
  saveRDS(seqtab_N_mash, file = '~/rundir/seqtab_N_mash_DATE.rds')
  
  
  seqtab_nochim_N_mash <- removeBimeraDenovo(seqtab_N_mash, method="consensus", multithread=TRUE, verbose=TRUE)
  dim(seqtab_nochim_N_mash)
  
  saveRDS(seqtab_N_mash, file = '~/rundir/seqtab_nochim_N_mash_DATE.rds')
  
  print(sum(seqtab_nochim_N_mash)/sum(seqtab_N_mash))
  
  track_N_mash <- cbind(N_out, trim_out_N_mash, sapply(dadaFs_N_mash, getN), sapply(dadaRs_N_mash, getN), sapply(mergers_N_mash, getN), rowSums(seqtab_nochim_N_mash)) # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track_N_mash) <- c("input", "N_filtered", "P_filtered", "trimmed", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track_N_mash) <- sample_names
  print(head(track_N_mash))
  
  saveRDS(track_N_mash, file = '~/rundir/track_N_mash_DATE.rds')
  
  taxa_N_mash <- assignTaxonomy(seqtab_nochim_N_mash, taxpath, taxLevels = c("Domain","Supergroup","Division","Kingdom","Phylum","Order","Family","Genus","Species"), multithread=TRUE)
  
  saveRDS(taxa_N_mash, file = '~/rundir/taxa_N_mash_DATE.rds')
  print("IMPORTANT NUMBERS BELOW")
  print("                         ")
  print(ncol(seqtab_nochim_N_mash))
  print(nrow(taxa_N_mash))
  
  taxa_print <- taxa_N_mash # Removing sequence rownames for display only
  rownames(taxa_print) <- NULL
  print(head(taxa_print))
  
  print("Contig Pipeline Finished")
}



print("Pipeline(s) Finished")


#______________________________________________________________________________________

#End of pipeline, all following lines of code are used to move outputs to the correct directory and name them. All have been tested and work. 

saveRDS(Parameters, file = "~/rundir/Parameters_DATE.rds") #Saves the run parameters so that each run can be replicated
file.copy(from = '~/pipelines/DADA2_pipeline_DATE.R', to = '~/rundir/DADA2_pipeline_DATE.R') #Moves the R code used to run the pipeline to the rundir so that it can be renamed and moved to the output. 

#Moving the stdout to the rundir to keep all run information together. Very janky method but it works. 
stdout <- readLines('~/outputs/stdout_pipeline_DATE.txt')
write.table(stdout, file = '~/rundir/stdout_pipeline_run_DATE.txt', row.names = FALSE, col.names = FALSE)

#Below code adds the date set at the start to all file names containing !, which should be all important.rds outputs and the stdout.

file.rename(list.files(pattern = "DATE"), str_replace(list.files(pattern = "DATE"),pattern = "DATE", Date))

#moving output files from rundir to dated output directory

outputs_to_be_moved <- list.files(rundir, pattern = Date) #Finds files with the correct date and prepares to move them

outputs_to_be_moved #Prints what outputs are about to be moved to the output directory

file.copy(from = outputs_to_be_moved, to = outpath, overwrite = TRUE) #Moves files from rundir to the dated output directory

if(proceed_to_formatting == TRUE){ #Proceeds on to immediately format the outputs generated (requires the formatting r script to be prepared for the data)
  system2("Rscript", args = ("~/pipelines/Pipeline_formatting_DATE.R"))
}

