remove(list=ls())
setwd("C:/Users/darth/WorkDir/Vsearch_pipeline")

source("Usearch_parameters3.R")
source("PipelineFunctions_abacusV3.txt")

library("ShortRead")
library("vegan")

load("C:/Users/darth/WorkDir/outputs/29_11_23/George_Cox_Masters2_OtuOutput.Rdat") #Change file path to your most recent pipeline output

joinedReads <- readFastq("C:/Users/darth/WorkDir/outputs/29_11_23/merged_reads.fq", sep="")

uniquesMatches <- match(as.character(sread(joinedReads)), uniquesfa$sequences)
seqData <- data.frame(ID = as.character(id(joinedReads)), uniqueSeq = uniquesfa$descs[uniquesMatches])

seqData$otu <- otuMatching$V2[match(seqData$ID, otuMatching$V1)] #something going on in here thats making the otu table lose seqs
seqData$zotu <- otuMatching$V2[match(seqData$ID, zotuMatching$V1)]


seqData$cleanID <- unlist(lapply(strsplit(seqData$ID, "\\."), function(x) x[1]))

resM <- table(seqData$cleanID[!is.na(seqData$otu)], seqData$otu[!is.na(seqData$otu)])

table2matrix <- function(table) { matrix(table, nrow=nrow(table),dimnames=list(rownames(table), colnames(table)))}

resM <- table2matrix(resM)

#mds1 <- metaMDS(resM[rowSums(resM) > 1000,])
#plot(mds1)

otus$otuMatch <- gsub("_sp$", "_sp\\.", otus$otuMatch) #Needed because for some reason the names would differ between otu table and database due to missing . after _sp
otus$otuMatch <- gsub("-sp$", "-sp\\.", otus$otuMatch)

otus_before <- otus


readDescriptions <- function(file)
    { ##just read description lines from database
    in.file <- readLines(file,  encoding="latin1")
    #in.file <- iconv(in.file, from="latin1", to="ASCII", ".")
    gsub(">", "", in.file[grep(">", in.file)])
    }

blastDbDescs <- strsplit(readDescriptions("C:/Users/darth/WorkDir/pr2_version_5.0.0_SSU_UTAX.fasta"), " ")


ID <- unlist(lapply(blastDbDescs,"[", 1))

blastDbDescs <- blastDbDescs[ID%in%otus$otuMatch]


DescsDF <- data.frame(ID = unlist(lapply(blastDbDescs,"[", 1)), taxonomy = unlist(lapply(blastDbDescs,function(x) paste(x[2:(length(x)-1)], collapse="y"))),
        SH = unlist(lapply(blastDbDescs,function(x) x[length(x)])))
otus <- merge(otus, DescsDF, by.x = 'otuMatch', by.y = 'ID', sort=FALSE) 


otus$identity <- as.numeric(otus$identity)
otus$length <- as.numeric(otus$length)

rownames(otus) <- otus$otu #makes the rownames of the otu table be the name of the otu. 

otus$Domain <- gsub(".*\\k\\:", "", otus$taxonomy) #Adds a new column with domain level taxonomy
otus$Domain <- gsub("\\,d.*", "", otus$Domain)

otus$Supergroup <- gsub(".*d\\:", "", otus$taxonomy)
otus$Supergroup <- gsub("\\,p.*", "", otus$Supergroup)

otus$Division <- gsub(".*p\\:", "", otus$taxonomy)
otus$Division <- gsub("\\,c.*", "", otus$Division) 


otus$Subdivision <- gsub(".*\\-", "", otus$Division)
otus$Division <- gsub("\\-.*", "", otus$Division)

otus$Class <- gsub(".*c\\:", "", otus$taxonomy)
otus$Class <- gsub("\\,o.*", "", otus$Class)

otus$Order <- gsub(".*o\\:", "", otus$taxonomy)
otus$Order <- gsub("\\,f.*", "", otus$Order)

otus$Family <- gsub(".*f\\:", "", otus$taxonomy)
otus$Family <- gsub("\\,g.*", "", otus$Family)

otus$Genus <- gsub(".*g\\:", "", otus$taxonomy)
otus$Genus <- gsub("\\,s.*", "", otus$Genus)

otus$Species <- gsub(".*s\\:", "", otus$taxonomy)

print(nrow(otus_before)) #Make sure these are the same
print(nrow(otus))




save(list=c("otus", "resM"), file="C:/Users/darth/WorkDir/OTU_folder/GC_OtuOutput2_tabulated.Rdat")

load(file = "C:/Users/darth/WorkDir/OTU_folder/GC_OtuOutput2_tabulated.Rdat")
