#Figure and stats on Forward only sequencing data


#TO DO:

#Fix the taxa colours but other than that everything else important works
#Key things, it doesnt represent the communtiy very well, also over estimates richness as it makes multiple ASVs of the same species (eg: phytophora and badhamia)
#Explains why the OTU data is better, more conservative. link to the long reads also potentially increasing the illumina error rates

#Venn diagram stuff

#Make pres look good.

#Reading in the data
#_______________________________________________________________________________

Inpath <- ("C:/Users/darth/Workdir/outputs/DADA2/12_9_23")
setwd(Inpath) 


Poscon_data <- read.csv("C:/USers/darth/Workdir/Positive_Control_Biomass_Sheet.csv")

Veges <- c("Beech", "Shrubland", "Pine")

resM <- readRDS("seqtab_nochim_forward_only_12_9_23.rds") #reads in the forward only pipeline

otus <- readRDS("taxa_forward_only_12_9_23.rds")

metadata <- read.csv("C:/USers/darth/Workdir/Field_Sampling_Metadata.csv")

setwd("C:/Users/darth/Workdir/Analysis_Folder_DADA2_FO")

set.seed(100)

#Packages
#_______________________________________________________________________________

library(vegan)
library(dplyr)
library(tidyr)

library(VennDiagram)
#Functions
#_______________________________________________________________________________
set.seed(100)

writeRepSetFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

rich <- function(x){
  rowSums(x>0)
}

#Data Preparation and filtering
#_______________________________________________________________________________

rep_set_ASVs <- as.data.frame(colnames(resM))
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:n())
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID) #renames the asvs to be the ASV_ID, the sub function uses regex to find and replace the name using the ^ wildcard.
rep_set_ASVs$ASV <- rep_set_ASVs$`colnames(resM)` #Adds the site/sample names to the ASV rep set
rep_set_ASVs$`colnames(resM)` <- NULL
colnames(resM) <- rep_set_ASVs$ASV_ID

colnames(otus) <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")


otus <- as.data.frame(otus)
otus$ASV <- as.factor(rownames(otus))
otus <- merge(rep_set_ASVs, otus, by = "ASV")
rownames(otus) <- otus$ASV_ID
otus$Species <- gsub("\\.", "", otus$Species) #removes all . (full stops) in the code so that downstream programs dont have issues with the fasta file. 
otus$Genus <- gsub("\\.", "", otus$Genus) #Same as above but for Genus



otus[is.na(otus)] <- "Unknown" #Next lien of code will remove all NAs so NAs are now known as Unknowns.
otus_filt <- subset(otus, Domain != "Bacteria" & Domain != "Archaea" & Subdivision != "Fungi" & Subdivision != "Metazoa" & Class != "Embryophyceae" & Domain != "Eukaryota:mito" & Domain != "Eukaryota:plas" & Domain != "Unknown")  #Filters out non-protists and mitochondrial/plastid sequences


to_remove <- anti_join(otus, otus_filt) #Finds out what is different between filtered and not filtered otus dfs.
to_remove <- rownames(to_remove) #Pulls out only the rownames (ASV_IDs)

resM_filt <- resM[ , !(colnames(resM) %in% to_remove)] #removes all the ASVs that did not pass the filter



#Bleedover checking code to look at data (could use sweep function to go through and remove things below a threshold)
#_______________________________________________________________________________

par(mfrow=c(5,4))
par(mar=c(1,1,1,1))
for(otu.i in names(sort(colSums(resM_filt), decreasing = TRUE))[1:20])
{
  plot(resM_filt[,otu.i], log = "y", main = otu.i)
  abline(h=max(resM_filt[,otu.i])*0.001, col="red")
}

#Remove bleedover for loop

resM_filt_bl <- resM_filt
bleedover <- resM_filt[1,]
for(otu.i in names(colSums(resM_filt_bl))){
  bleedover[otu.i] <- as.matrix(max(resM_filt_bl[,otu.i])*0.001)
  bleedover <- round(bleedover, digits = 0)
  resM_filt_bl[,otu.i] <- as.matrix(pmax(0, resM_filt_bl[,otu.i]-bleedover[otu.i]))
} #This loop removes reads that are likely to be bleedover

resM_filt_bl[which(apply(resM_filt_bl, 2, function(x) x == max(x, na.rm = TRUE)))] <- resM_filt[which(apply(resM_filt, 2, function(x) x == max(x, na.rm = TRUE)))]
#This line of code restores the max values to their original values. 


for(otu.i in names(sort(colSums(resM_filt_bl), decreasing = TRUE))[1:20]){
  plot(resM_filt_bl[,otu.i], log = "y", main = otu.i)
  abline(h=max(resM_filt[,otu.i])*0.001, col="red")
}




#Subsetting by vegetation type
row_names_fix <- rownames(resM_filt[c(1:33),])
row_names_fix <- append(row_names_fix, c("GCFLDCON", "GCLABCON", "GCNEGCON", "GCPOSCON"))#Fixes the rownames so that the controls are easier to pull out in the vegtype
rownames(resM_filt_bl) <- row_names_fix




#Back up subsetting 

pinesites <- c(3, 6, 8, 12, 13, 18, 20, 23, 26, 29, 31)
shrubsites <- c(1, 5, 9, 11, 14, 17, 21, 24, 25, 30, 33)
nothosites <- c(2, 4, 7, 10, 15, 16, 19, 22, 27, 28, 32)
controls <- c(34, 35, 36 ,37)


#Removing controls for later analysis and making resM_prop
resM_nc <- resM_filt[1:33, ] #No controls

resM_nc <- resM_nc[,colSums(resM_nc) > 0]

#saveRDS(resM_fin, "GC_OTU_community_matrix.RDS")




#Colours
palette.colors(palette = "Okabe-Ito")

pinecol <-  "#009E73"
nothocol <- "#D55E00" 
shrubcol <- "#56B4E9"  #Change these colours around

pinepch <- 21
nothopch <- 22
shrubpch <- 23

vegpch <- c(nothopch, shrubpch, pinepch)

vegcol <- c(nothocol,  shrubcol, pinecol)

Amoebocol <- "#F0E442"
cercocol <- "#56B4E9"
apicol <- "#CC79A7"
cilicol <- "#E69F00"
othercol <- "black"
unkcol <- "white"

amoebpch <- 21
cercopch <- 25
apipch <- 24
cilipch <- 23
otherpch <- 22
unkpch <- 13

otus_filt$colour <- othercol


otus_filt$colour[otus_filt$Subdivision == "Cercozoa"] <- cercocol
otus_filt$colour[otus_filt$Subdivision == "Apicomplexa"] <- apicol
otus_filt$colour[otus_filt$Subdivision == "Ciliophora"] <- cilicol
otus_filt$colour[otus_filt$Subdivision == "Unknown"] <- unkcol
otus_filt$colour[otus_filt$Supergroup == "Amoebozoa"] <- Amoebocol

otus_filt$pch <- otherpch

otus_filt$pch[otus_filt$Subdivision == "Cercozoa"] <- cercopch
otus_filt$pch[otus_filt$Subdivision == "Apicomplexa"] <- apipch
otus_filt$pch[otus_filt$Subdivision == "Ciliophora"] <- cilipch
otus_filt$pch[otus_filt$Subdivision == "Unknown"] <- unkpch
otus_filt$pch[otus_filt$Supergroup == "Amoebozoa"] <- amoebpch


key_tax <- c("Amoebozoa","Apicomplexa", "Cercozoa", "Ciliophora", "Other Division")
tax_col <- c(Amoebocol, apicol, cercocol, cilicol, othercol)
tax_pch <- c(amoebpch, apipch, cercopch, cilipch, otherpch)

col_df <- data.frame(key_tax, tax_col)


#Taxonomic assignment plots and filtering by database matches


nrow(otus_filt[otus_filt$Family == "Unknown",])



resM_Cons <- resM_filt_bl[ , !(colnames(resM_filt_bl)) %in% rownames(otus_filt[otus_filt$Family == "Unknown",])]

resM_fin <- resM_Cons[1:33, ] #No controls

resM_fin <- resM_fin[,colSums(resM_fin) > 0]

otus_fin <- otus_filt[otus_filt$Family != "Unknown",]

nrow(otus_filt)
nrow(otus_fin)
ncol(resM_fin)

setdiff(rownames(otus_fin), colnames(resM_fin)) 

to_remove3 <- setdiff(rownames(otus_fin), colnames(resM_fin)) #removes controls from otus data frame
otus_fin <- otus_fin[!(rownames(otus_fin) %in% to_remove3),]

nrow(otus_fin)
ncol(resM_fin)

resM_prop <- sweep(resM_fin, 1, rowSums(resM_fin), "/")*100 #Makes a table of proportions. 



#rarefaction

site_reads <- as.data.frame(rowSums(resM_fin))
colSums(site_reads)
mean(site_reads$`rowSums(resM_fin)`)

n_read <- site_reads[nothosites,]
s_read <- site_reads[shrubsites,]
p_read <- site_reads[pinesites,]

mean_reads <- c(mean(n_read), mean(s_read), mean(p_read))



rared <- rarefy(resM_fin, sample = 3464)
rared


#rarecurve(resM_fin, sample = 4164)
resM_rare <- rrarefy(resM_fin, sample = 3464) #Used for richness

reps <- 100
# Create an empty array (3 dimensional matrix) of the same size as data, but with multiple copies (or layers):
rareStack <- array(NA, dim=c(nrow(resM_fin), ncol(resM_fin), reps))

## Put a rarefied version of the data into each layer of the matrix:
for(stack.i in 1:reps){
  rareStack[,,stack.i] <- rrarefy(resM_fin, 3464)
} #USed for effective diversity.









#visualising controls
par(mfrow =c(1,1))
options(scipen = 999)

controls <- resM_Cons[34:37,]
controls <- controls[, colSums(controls) > 0]
colnames(controls) <- paste0(otus_filt[colnames(controls),]$Family, " ", otus_filt[colnames(controls),]$Species)

poscon_test <- as.data.frame(controls[4,], colnames(controls))
sort(controls[4,])
pie(sort(controls[4,]))
sum(controls[4,])

control_prop <- controls[4,]/sum(controls[4,])*100
control_prop <- control_prop[control_prop>0]
sort(control_prop)

length(control_prop)



#Read and OTU counting

ncol(resM_fin)
nrow(otus_fin) 

site_reads <- as.data.frame(rowSums(resM_fin))
colSums(site_reads)
mean(site_reads$`rowSums(resM_fin)`)

site_OTUs <- as.data.frame(rich(resM_fin))
mean(site_OTUs$`rich(resM_fin)`)

site_reads_OTUS <- as.data.frame(site_reads)
site_reads_OTUS$OTUs <- site_OTUs
colnames(site_reads_OTUS) <- c("Reads", "OTUs")
site_reads_OTUS

read_cor <- lm(unlist(site_OTUs) ~ unlist(site_reads))
summary(read_cor)



#Do I include this or not?


vegtype <- substr(row_names_fix[1:33],6, 6) #Subsets the data by the vegetation type
vegtype
vegtype <- factor(vegtype, levels = c("S", "N", "P")) #Adds the letters as levels so that each vegetype and controls can be pulled out.
vegtype
resM_fin[vegtype == "S"]





#Plots
#_______________________________________________________________________________
par(mar=c(7,7,7,7))
par(mfrow =c(1,1))




#Community description
#Total

length(unique(otus_fin$Supergroup))


subs <- otus_fin[colnames(resM_fin),]$Subdivision

sds <- resM_fin
colnames(sds) <- subs
sds <- sds %*% sapply(unique(subs), "==", subs)


groups <- as.data.frame(sort(colSums(sds), decreasing = TRUE))
group_prop <- ((groups/colSums(groups))*100)

nrow(group_prop)

#Beech
subs_n <- otus_fin[colnames(resM_fin[nothosites,]),]$Subdivision

sds_n <- resM_fin[nothosites,]
colnames(sds_n) <- subs
sds_n <- sds_n %*% sapply(unique(subs_n), "==", subs_n)


groups_n <- as.data.frame(sort(colSums(sds_n), decreasing = TRUE))
group_prop_n <- ((groups_n/colSums(groups_n))*100)



#Shrubland
subs_s <- otus_fin[colnames(resM_fin[shrubsites,]),]$Subdivision

sds_s <- resM_fin[shrubsites,]
colnames(sds_s) <- subs_s
sds_s <- sds_s %*% sapply(unique(subs_s), "==", subs_s)


groups_s <- as.data.frame(sort(colSums(sds_s), decreasing = TRUE))
group_prop_s <- ((groups_s/colSums(groups_s))*100)


#Pine
subs_p <- otus_fin[colnames(resM_fin[pinesites,]),]$Subdivision
sds_p <- resM_fin[pinesites,]
colnames(sds_p) <- subs_p
sds_p <- sds_p %*% sapply(unique(subs_p), "==", subs_p)


groups_p <- as.data.frame(sort(colSums(sds_p), decreasing = TRUE))
group_prop_p <- ((groups_p/colSums(groups_p))*100)


#groups within cercozoa

otus_cerconly <- otus_fin[otus_fin$Subdivision == "Cercozoa",]

to_remove4 <- anti_join(otus_fin, otus_cerconly)
to_remove4 <- rownames(to_remove4)
resM_cercs <- resM_prop[ , !(colnames(resM_fin) %in% to_remove4)]
ncol(resM_cercs)

cercs <- otus_cerconly[colnames(resM_cercs),]$Class

colnames(resM_cercs) <- cercs
resM_cercs <- resM_cercs %*% sapply(unique(cercs), "==", cercs)

colMeans(resM_cercs[pinesites,])
colMeans(resM_cercs[nothosites,])
colMeans(resM_cercs[shrubsites,])






options(scipen = 0)

#OTUs by group

nrow(otus_fin[otus_fin$Subdivision == "Cercozoa",])
nrow(otus_fin[otus_fin$Subdivision == "Ciliophora",])
nrow(otus_fin[otus_fin$Subdivision == "Apicomplexa",])
nrow(otus_fin[otus_fin$Supergroup == "Amoebozoa",])












#shannon/effective diversity (needs rarestack) (Done)

shannon <- apply(rareStack, 3, function(x) exp(diversity(x, index = "shannon")))

rownames(shannon) <- rownames(resM_fin)

shannon_mean <- rowMeans(shannon)

p_shannon <- shannon_mean[vegtype == "P"]
n_shannon <- shannon_mean[vegtype == "N"]
s_shannon <- shannon_mean[vegtype == "S"]

shannon_list <- c(n_shannon, s_shannon, p_shannon) #FIX THIs
shannon_avg<- c(mean(n_shannon), mean(s_shannon), mean(p_shannon))

sitelist <- names(shannon_list)
sitelist <- gsub("GC...", "", sitelist)
shan_aov_df <- data.frame(shannon_list, sitelist)
colnames(shan_aov_df) <- c("EffectiveDiversity", "Vegetation")
shan_aov <- aov(EffectiveDiversity ~ Vegetation, data = shan_aov_df)
summary(shan_aov)
TukeyHSD(shan_aov)
shannon_avg
summary(lm(shan_aov))

shannon_se <-c(23.802, 33.661, 33.661)

pdf(file = "ASV Mean Effective Diversity plot.pdf")

shannon_barplot <- barplot(shannon_avg, main = "Effective diversity of each vegetation type", names = Veges, ylab = "Effective ASV Diversity", ylim = c(0,500), col = "white", cex.names = 1.5)
arrows(x0 = shannon_barplot,
       y0 = shannon_avg + shannon_se,
       y1 =shannon_avg - shannon_se,
       angle = 90, code = 3, length = 0.1)
points(x = lapply(rep(shannon_barplot, each=11), jitter, amount=0.05),
       y = shannon_list,
       col = "black", pch=21)

dev.off()





#richness (rarefied)


richness <- rich(resM_rare)

p_rich_mean <- mean(richness[vegtype == "P"])
n_rich_mean <- mean(richness[vegtype == "N"])
s_rich_mean <- mean(richness[vegtype == "S"])

rich_list <- c(rich(resM_rare[nothosites,]), rich(resM_rare[shrubsites, ]), rich(resM_rare[pinesites,]))
mean_rich <- c(n_rich_mean, s_rich_mean, p_rich_mean)



sitelist <- names(rich_list)
sitelist <- gsub("GC...", "", sitelist)
rich_aov_df <- data.frame(rich_list, sitelist)
colnames(rich_aov_df) <- c("Richness", "Vegetation")
rich_aov <- aov(Richness ~ Vegetation, data = rich_aov_df)
summary(rich_aov)
TukeyHSD(rich_aov)
summary(lm(rich_aov))

mean_rich
rich_se <- c(43.45, 61.45, 61.45)

pdf(file = "ASV Mean Richness plot.pdf")

rich_barplot <- barplot(mean_rich, names.arg = Veges, ylab = "Rarefied ASV Richness", main = "Rarefied richness of each vegetation type", ylim = c(0,1000), col = "white", cex.names =1.5)
arrows(x0 = rich_barplot,
       y0 = mean_rich + rich_se,
       y1 = mean_rich - rich_se,
       angle = 90, code = 3, length = 0.1)
points(x = lapply(rep(rich_barplot, each=11), jitter, amount=0.05),
       y = rich_list,
       col = "black", pch=21)
dev.off()





#beta diversity
#beta_dist <- vegdist(resM_fin, method = "jaccard")
beta_dist2 <- raupcrick(resM_fin, nsimul = 100) # TO DO: should be at least 99, read chase paper (if still relevant)
#Raup Crick makes a dissimilarity matrix that is less impacted by differences in alpha diversity (which my data is very effected by) when compared to other distances measurements. 

#Removing site 33 does make the raup crick work and repeat the best solution. 

par(mfrow=c(1,1))


beta_mds <- metaMDS(beta_dist2, k =2, trymax = 1000)



mds_data <- as.data.frame(beta_mds$points)

#pdf(file = "ASV raupcrick_beta_diversity.pdf")

plot(mds_data, main = "Community dissimilarity of the three vegetation types")

vegtype_MDS <- substr(rownames(beta_mds$points),6,6)
vegtype_MDS <- factor(vegtype_MDS, levels = c("N", "S", "P"))


#plot(beta_mds$points, vegtype_MDS, col = vegcol, pch = vegpch, main = "Beta Diversity of the Three Vegetation Types") #Colour points by vegtype. 
ordihull(beta_mds, vegtype_MDS, col = vegcol, lwd = 3)
#ordiellipse(beta_mds, vegtype_MDS, col = vegcol)
points(beta_mds$points, bg = vegcol[as.factor(vegtype_MDS)], pch = vegpch[as.factor(vegtype_MDS)])
legend(-1.5,-0.255, legend = Veges, col = vegcol, pch = vegpch, cex = 1.4, pt.bg = vegcol)

#ordispider(beta_mds, substr(rownames(beta_mds$points),6,6))
#orditorp(beta_mds,display="sites",cex=0.5,air=0.01)



#dev.off()
beta_mds
stressplot(beta_mds)

beta_disp <- betadisper(beta_dist2, substr(rownames(beta_mds$points),6,6))
anova(beta_disp)
TukeyHSD(beta_disp)



#Composition plot

resM_dom <- resM_prop[, colSums(resM_prop) > 0.1] #

#pdf(file = 'compositionplot.pdf')




bar_order <- order(colSums(resM_dom[nothosites, ]),
                   colSums(resM_dom[shrubsites,]),
                   colSums(resM_dom[pinesites,]), decreasing = FALSE)
par(mfrow=c(1,3))

par(mar=c(0,1,0,0))
par(mai=c(1,0.5,1,0))

barplot(colSums(resM_dom[nothosites,bar_order])/11, main = "Beech", log = "x", offset = 0.1, ylab = "Operational Taxonomic Units (OTUs)",
        horiz = TRUE, col = otus_fin$colour[match(colnames(resM_dom[nothosites,bar_order]),otus_fin$ASV_ID)], border = NA,
        cex.names = 2, names.arg = "")



barplot(colSums(resM_dom[shrubsites,bar_order])/11, main = "Shrubland", log = "x", offset = 0.1, xlab = "Dominance (%)",
        horiz = TRUE, col = otus_fin$colour[match(colnames(resM_dom[nothosites,bar_order],),otus_fin$ASV_ID)], border = NA,
        names.arg = "")

barplot(colSums(resM_dom[pinesites,bar_order])/11, main = "Pine", log = "x", offset = 0.1, xlab = "",
        horiz = TRUE, col = otus_fin$colour[match(colnames(resM_dom[nothosites,bar_order]),otus_fin$ASV_ID)], border = NA, 
        names.arg = "") 


#dev.off()










#Rarity stuff
#_____________________________________________________________________________-



#Occupancy
#Ian had low occupancy being present in no more than 5% of sites (4 out of 81), so if I use a similar metric and round up I get 2 sites. 
frequency <- colSums(resM_prop > 0)

pine_freq <- colSums(resM_prop[pinesites, ] > 0)
pine_freq <- pine_freq[pine_freq > 0] #removes ones not present in the pine sites

notho_freq <- colSums(resM_prop[nothosites, ] > 0)
notho_freq <- notho_freq[notho_freq > 0]

shrub_freq <- colSums(resM_prop[shrubsites, ] > 0)
shrub_freq <- shrub_freq[shrub_freq > 0] 


#Dominance
#low dominance = less than 1% of all reads maximum across all occurrences


dominance <- colSums(resM_prop)/frequency

pine_dom <- colSums(resM_prop[pinesites,])/frequency
pine_dom <- pine_dom[pine_dom > 0]


notho_dom <- colSums(resM_prop[nothosites,])/frequency
notho_dom <- notho_dom[notho_dom > 0]

shrub_dom <- colSums(resM_prop[shrubsites,])/frequency
shrub_dom <- shrub_dom[shrub_dom > 0]




#how to do some stats on this? By taxa or by site?
dom_freq <- lm(dominance ~ frequency)
summary(dom_freq) #Maybe needs revisting

cor.test(dominance, frequency, formula = "pearsons")

dom_freq_p <- lm(pine_dom ~ pine_freq)
summary(dom_freq_p)


#domanincae by frequency for different taxa

#table (rowname = OTU), dominance, frequency, taxa.

rare_df <- as.data.frame(dominance)
rare_df$frequency <- frequency
rare_df$group <- otus_fin$Subdivision[match(rownames(rare_df),otus_fin$otu)]
rare_df$species <- otus_fin$Species[match(rownames(rare_df),otus_fin$otu)]

rare_notho <- as.data.frame(notho_dom)
nfreq <- colSums(resM_prop[nothosites,] > 0)
rare_notho$frequency <- nfreq[nfreq > 0]
rare_notho$group <- otus_fin$Subdivision[match(rownames(rare_notho),otus_fin$otu)]
rare_notho$species <- otus_fin$Species[match(rownames(rare_notho),otus_fin$otu)]

rare_shrub <- as.data.frame(shrub_dom)
sfreq <- colSums(resM_prop[shrubsites,] > 0)
rare_shrub$frequency <- sfreq[sfreq > 0]
rare_shrub$group <- otus_fin$Subdivision[match(rownames(rare_shrub),otus_fin$otu)]
rare_shrub$species <- otus_fin$Species[match(rownames(rare_shrub),otus_fin$otu)]

rare_pine <- as.data.frame(pine_dom)
pfreq <- colSums(resM_prop[pinesites,] > 0)
rare_pine$frequency <- pfreq[pfreq > 0]
rare_pine$group <- otus_fin$Subdivision[match(rownames(rare_pine),otus_fin$otu)]
rare_pine$species <- otus_fin$Species[match(rownames(rare_pine),otus_fin$otu)]

dom_taxa <- rare_df[rare_df$dominance > 1,]
View(dom_taxa)
nrow(dom_taxa)
nrow(dom_taxa[dom_taxa$group == "Cercozoa",])
unique(dom_taxa$group)
dom_taxa[dom_taxa$frequency < 3,]




#Plot of frequency (x axis) vs dominance (y axis) (make a plot for each vegtype). 
#pdf(file = 'rarityplot.pdf', )


par(mar=c(4.3,4.5,2,1))
par(cex.lab = 1.7)
par(mfrow = c(1,3))


#overall_freq_vs_dom <- plot(dominance ~ jitter(frequency), log = "xy", col = otus_filt$colour[match(colnames(resM_prop),
#otus_filt$otu)], pch = 16, main = "All Sites",  xlab = "Frequency", ylab = "Dominance")


#abline(h = 1, col = "red")



notho_freq_vs_dom <- plot(notho_dom ~ jitter(notho_freq), log = "xy", bg = otus_fin$colour[match(names(notho_dom), otus_fin$otu)],
                          pch = otus_fin$pch[match(names(notho_dom), otus_fin$otu)], main = "Beech", ,  xlab = "", ylab = "Dominance (%)", cex = 2, cex.main = 2)

abline(h = 1, col = "red", lwd =  1.5)



shrub_freq_vs_dom <- plot(shrub_dom ~ jitter(shrub_freq), log = "xy", bg = otus_fin$colour[match(names(shrub_dom),otus_fin$otu)], 
                          pch = otus_fin$pch[match(names(shrub_dom), otus_fin$otu)], main = "Shrubland",  xlab = "Frequency", ylab = "", cex = 2, cex.main = 2)

abline(h = 1, col = "red", lwd =  1.5)


pine_freq_vs_dom <- plot(pine_dom ~ jitter(pine_freq), log = "xy", bg = otus_fin$colour[match(names(pine_dom),otus_fin$otu)], 
                         pch = otus_fin$pch[match(names(pine_dom),otus_fin$otu)], main = "Pine",  xlab = "", ylab = "", cex =2, cex.main = 2)

abline(h = 1, col = "red", lwd =  1.5)





#plot.new()
par(mar=c(1,1,1,1))

legend('center', legend = key_tax, pt.bg = tax_col, pch = tax_pch, cex = 2)


dev.off()





dom_taxa <- sort(dominance, decreasing = TRUE)
dom_taxa <- dom_taxa[dom_taxa > 1]
names(dom_taxa) <- otus_fin[names(dom_taxa),]$Division
dom_taxa

pine_dom_top <- sort(pine_dom, decreasing = TRUE)
pine_dom_top <- pine_dom_top[pine_dom_top > 0.5]
names(pine_dom_top) <- otus_fin[names(pine_dom_top),]$Division
pine_dom_top

#ASVs unique to vegtypes


resM_mushed_pine <- as.matrix(colSums(resM_fin[pinesites,]))
resM_mushed_notho <- as.matrix(colSums(resM_fin[nothosites,]))
resM_mushed_shrub <- as.matrix(colSums(resM_fin[shrubsites,]))

resM_mushed <- cbind(colnames(resM_fin), resM_mushed_notho, resM_mushed_shrub, resM_mushed_pine)
colnames(resM_mushed) <- c("otu", Veges) #Names needed their own actual column in matrix for for loops to read the names to a list



#Make empty list
pine_uniques <- list()
notho_uniques <- list()
shrub_uniques <- list()

#For loop that finds sequences that are only present in one vegtype and adds to a list, all 3 need to run (one for each vegtype).
for(otu.i in rownames(resM_mushed)){
  if(resM_mushed[otu.i, 2] > 0){
    if(resM_mushed[otu.i, 3] == 0){
      if(resM_mushed[otu.i, 4] == 0){
        notho_uniques <- append(notho_uniques, resM_mushed[otu.i,1])
      }
    }
  }
}

for(otu.i in rownames(resM_mushed)){
  if(resM_mushed[otu.i, 2]== 0){
    if(resM_mushed[otu.i, 3] > 0){
      if(resM_mushed[otu.i, 4] == 0){
        shrub_uniques <- append(shrub_uniques, resM_mushed[otu.i,1])
      }
    }
  } 
}

for(otu.i in rownames(resM_mushed)){
  if(resM_mushed[otu.i, 2] == 0){
    if(resM_mushed[otu.i, 3] == 0){
      if(resM_mushed[otu.i, 4] > 0){
        pine_uniques <- append(pine_uniques, resM_mushed[otu.i,1])
      }
    }
  }
}

OTUs_unique_to_veg <- c(length(notho_uniques), length(shrub_uniques), length(pine_uniques), sum(c(length(pine_uniques), length(notho_uniques), length(shrub_uniques))))
barplot(OTUs_unique_to_veg, names.arg = c(Veges, "Total"), main = "Number of OTUs unique to each vegetation type")



#Shared ASVs
all_shared <- list()

pn_shared <- list()

ps_shared <- list()

ns_shared <- list()

for(otu.i in rownames(resM_mushed)){
  if(resM_mushed[otu.i, 2] > 0){
    if(resM_mushed[otu.i, 3] > 0){
      if(resM_mushed[otu.i, 4] > 0){
        all_shared <- append(all_shared, resM_mushed[otu.i,1])
      }
    }
  }
}

for(otu.i in rownames(resM_mushed)){
  if(resM_mushed[otu.i, 2] > 0){
    if(resM_mushed[otu.i, 3] > 0){
      if(resM_mushed[otu.i, 4] == 0){
        pn_shared <- append(pn_shared, resM_mushed[otu.i,1])
      }
    }
  }
}

for(otu.i in rownames(resM_mushed)){
  if(resM_mushed[otu.i, 2] > 0){
    if(resM_mushed[otu.i, 3] == 0){
      if(resM_mushed[otu.i, 4] > 0){
        ps_shared <- append(ps_shared, resM_mushed[otu.i,1])
      }
    }
  }
}

for(otu.i in rownames(resM_mushed)){
  if(resM_mushed[otu.i, 2] == 0){
    if(resM_mushed[otu.i, 3] > 0){
      if(resM_mushed[otu.i, 4] > 0){
        ns_shared <- append(ns_shared, resM_mushed[otu.i,1])
      }
    }
  }
}

shared_OTUs <- c(length(pn_shared), length(ps_shared), length(ns_shared),length(all_shared))
shared_OTUs

ncol(resM_fin)
c(length(pine_uniques), length(notho_uniques), length(shrub_uniques))
#MAKE A VENN DIAGRAM



shared_OTUs
OTUs_unique_to_veg[1:3]

# create Venn diagram and display all sets 



library(VennDiagram)



pine_freq <- colSums(resM_prop[pinesites, ] > 0)
pine_freq <- pine_freq[pine_freq > 0] #removes ones not present in the pine sites
pine_OTUs <- names(pine_freq)

notho_freq <- colSums(resM_prop[nothosites, ] > 0)
notho_freq <- notho_freq[notho_freq > 0]
notho_OTUs <- names(notho_freq)

shrub_freq <- colSums(resM_prop[shrubsites, ] > 0)
shrub_freq <- shrub_freq[shrub_freq > 0] 
shrub_OTUs <- names(shrub_freq)


# Chart
venn.diagram(
  x = list(notho_OTUs, shrub_OTUs, pine_OTUs),
  category.names = c("Beech" , "Shrubland" , "Pine"),
  filename = 'test_venn.png',
  output=TRUE,
  disable.logging = TRUE,
  print.mode = c('raw', 'percent')
)




#COpying Ians table from Deadwood fungal rarity paper
colnames_for_table <- c(Veges, "Total")
sites <- c(length(nothosites), length(shrubsites), length(pinesites), 33)
total_OTUs <- c(length(notho_OTUs),length(shrub_OTUs), length(pine_OTUs), ncol(resM_fin))
low_occupancy <- c(length(notho_freq[notho_freq < 3]),length(shrub_freq[shrub_freq < 3]),length(pine_freq[pine_freq < 3]), length(frequency[frequency < 3]))
low_dominance <- c(length(notho_dom[notho_dom < 1]), length(shrub_dom[shrub_dom < 1]), length(pine_dom[pine_dom < 1]), length(dominance))

rarity_table <- data.frame(sites, total_OTUs, low_occupancy, low_dominance, OTUs_unique_to_veg, shared_OTUs)
rarity_table <- t(rarity_table)
colnames(rarity_table) <- colnames_for_table

rarity_table

#COuld do some stats on number of unique OTUs to shrubland


#NPMANOVA


test_NPMANOVA <- adonis2(resM_fin ~ as.factor(metadata$Type[1:33]))
test_NPMANOVA