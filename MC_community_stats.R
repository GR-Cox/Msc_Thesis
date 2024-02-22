#MC community Analysis Script

#Loading data

setwd("C:/Users/darth/Workdir/MC_stats/")

resM <- read.csv("MC_resM.csv")
rownames(resM) <- resM$X
resM$X <- NULL



metadata <- read.csv("C:/USers/darth/Workdir/Field_Sampling_Metadata.csv")


  
mc_sheet <- read.csv("mc_sheet.csv", row.names = 1)


OTUs <- read.csv("Seqs_to_OTUs.csv")





#packages
library(vegan)

#Functions

set.seed(100)

rich <- function(x){
  rowSums(x>0)
}

#Setup code

vegtype <- substr(rownames(resM),6, 6) #Subsets the data by the vegetation type
vegtype <- factor(vegtype, levels = c("S", "N", "P")) #Adds the letters as levels so that each vegetype and controls can be pulled out.

pinesites <- c(3, 6, 8, 12, 13, 18, 20, 23, 26, 29, 31)
shrubsites <- c(1, 5, 9, 11, 14, 17, 21, 24, 25, 30, 33)
nothosites <- c(2, 4, 7, 10, 15, 16, 19, 22, 27, 28, 32)

Veges <- c("Beech", "Shrubland", "Pine")

pinecol <-  "#009E73"
nothocol <- "#D55E00" 
shrubcol <- "#56B4E9"  #Change these colours around

pinepch <- 21
nothopch <- 22
shrubpch <- 23

vegpch <- c(nothopch, shrubpch, pinepch)

vegcol <- c(nothocol,  shrubcol, pinecol)

#Analysis code






#Richness

richness <- rich(resM)

p_rich_mean <- mean(richness[vegtype == "P"])
n_rich_mean <- mean(richness[vegtype == "N"])
s_rich_mean <- mean(richness[vegtype == "S"])

rich_list <- c(rich(resM[vegtype == "N",]), rich(resM[vegtype == "S",]), rich(resM[vegtype == "P",]))

sitelist <- names(rich_list)
sitelist <- gsub("GC...", "", sitelist)
rich_aov_df <- data.frame(rich_list, sitelist)
colnames(rich_aov_df) <- c("Richness", "Vegetation")

glm_rich <- glm(Richness ~ Vegetation, family = "poisson", data = rich_aov_df)
summary(aov(glm_rich))
summary(glm_rich)


beech_log <- 0.5978
beech_log_se <- 0.2236

shrub_log <- 0.3001+0.5978
shrub_log_se <- 0.2950

pine_log <- -0.2231+0.5978
pine_log_se <- 0.3354

rich_mean_log <- c(beech_log, shrub_log, pine_log)
rich_mean <- exp(rich_mean_log)

rich_se_log <- c(beech_log_se, shrub_log_se, pine_log_se)

#pdf(file = "MC Mean Richness plot.pdf")

rich_barplot <- barplot(rich_mean, names.arg = Veges, ylab = "OTU Richness", main = "Richness of each Vegetation Type", ylim = c(0,6), col = "white", cex.names =1.5)
arrows(x0 = rich_barplot,
       y0 = exp(rich_mean_log + rich_se_log),
       y1 = exp(rich_mean_log - rich_se_log),
       angle = 90, code = 3, length = 0.1)
points(x = lapply(rep(rich_barplot, each=11), jitter, amount=0.1),
       y = rich_list,
       col = "black", pch=21)
#dev.off()



#Effective diversity

eff_div <- exp(diversity(resM, index = "shannon"))


p_eff_div <- eff_div[vegtype == "P"]
n_eff_div <- eff_div[vegtype == "N"]
s_eff_div <- eff_div[vegtype == "S"]

eff_div_list <- c(n_eff_div, s_eff_div, p_eff_div) 
eff_div_avg<- c(mean(n_eff_div), mean(s_eff_div), mean(p_eff_div))
eff_div_avg
eff_div_se <- c(sd(n_eff_div)/sqrt(11), sd(s_eff_div)/sqrt(11), sd(p_eff_div)/sqrt(11))

sitelist <- names(eff_div_list)
sitelist <- gsub("GC...", "", sitelist)
shan_aov_df <- data.frame(eff_div_list, sitelist)
colnames(shan_aov_df) <- c("EffectiveDiversity", "Vegetation")


hist(eff_div) #Not normally distributed and can't use poisson as there are non-integers

#Test of which type of transformation to do from BIOL 309 lab manual
shan_aov <- aov(EffectiveDiversity ~ Vegetation, data = shan_aov_df)
plot(shan_aov)
test <- lm(abs(shan_aov$residuals)~shan_aov$fitted.values)
summary(test) #2 x the slope is closest to 1 so should log transform.

eff_div_log <- log(eff_div_list +1)
hist(eff_div_log) #more normal but still not great
eff_div_df_log <- data.frame(eff_div_log, sitelist)
colnames(eff_div_df_log) <- c("LogEffectiveDiversity", "Vegetation")
eff_div_log_aov <- aov(LogEffectiveDiversity ~ Vegetation, data = eff_div_df_log)
plot(eff_div_log_aov)

summary(eff_div_log_aov)
TukeyHSD(eff_div_log_aov)

effdivlm <- lm(eff_div_log_aov)
summary(effdivlm) #Used to see intercepts to build plot via back transformation

beech_eff_log <- 1.0436
beech_eff_log_se <- 0.1018

shrub_eff_log <- 1.0436+0.1154
shrub_eff_log_se <- 0.1440

pine_eff_log <- 1.0436-0.1682
pine_eff_log_se <- 0.1440

plot_eff_div_log <- c(beech_eff_log, shrub_eff_log, pine_eff_log)
plot_eff_div_log_se <- c(beech_eff_log_se, shrub_eff_log_se, pine_eff_log_se)

#pdf(file = "Mean MC Effective Diversity plot.pdf")

eff_div_barplot <- barplot(exp(plot_eff_div_log)-1, main = "Effective Diversity of OTUs of each Vegetation Type", names = Veges, ylab = "Effective Diversity of OTUs", ylim = c(0,5.5), col = "white", cex.names = 1.5)
arrows(x0 = eff_div_barplot,
       y0 = exp(plot_eff_div_log + plot_eff_div_log_se)-1,
       y1 = exp(plot_eff_div_log - plot_eff_div_log_se)-1,
       angle = 90, code = 3, length = 0.1)
points(x = lapply(rep(eff_div_barplot, each=11), jitter, amount=0.1),
       y = eff_div_list,
       col = "black", pch=21)


#dev.off()






#Composition plot

resM_prop <- sweep(resM, 1, rowSums(resM), "/")*100
resM_prop[c(2, 12, 23),] <- 0

resM_dom <- resM[, colSums(resM) > 0.1] 

#pdf(file = 'compositionplot.pdf')




bar_order <- order(colSums(resM[vegtype == "N", ]),
                   colSums(resM[vegtype == "S",]),
                   colSums(resM[vegtype == "P",]), decreasing = FALSE)
par(mfrow=c(1,3))

par(mar=c(0,1,0,0))
par(mai=c(1,0.5,1,0))

barplot(colSums(resM[vegtype == "N",bar_order]), axisnames = FALSE, main = "Beech", ylab = "Slime Mould OTUs",
        horiz = TRUE)



barplot(colSums(resM[vegtype == "S",bar_order]), axisnames = FALSE, main = "Shrubland", xlab = "Frequency",
        horiz = TRUE)

barplot(colSums(resM[vegtype == "P",bar_order]), axisnames = FALSE, main = "Pine", xlab = "",
        horiz = TRUE)
#plot.new()


#dev.off()



#Beta diversity


#Removing site 33 does make the raup crick work and repeat the best solution. 

par(mfrow=c(1,1))
rowSums(resM)

MC_dist <- vegdist(resM[-c(2,12,23,14,28,31,33),])


beta_mds <- metaMDS(MC_dist, k =2, trymax = 1000)



mds_data <- as.data.frame(beta_mds$points)

#pdf(file = "raupcrick_beta_diversity.pdf")

plot(mds_data, main = "GC-OTU Community Dissimilarity of the Three Vegetation Types")

vegtype_MDS <- substr(rownames(beta_mds$points),6,6)
vegtype_MDS <- factor(vegtype_MDS, levels = c("N", "S", "P"))


#plot(beta_mds$points, vegtype_MDS, col = vegcol, pch = vegpch, main = "Beta Diversity of the Three Vegetation Types") #Colour points by vegtype. 
ordihull(beta_mds, vegtype_MDS, col = vegcol, lwd = 3)
#ordiellipse(beta_mds, vegtype_MDS, col = vegcol)
points(beta_mds$points, bg = vegcol[as.factor(vegtype_MDS)], pch = vegpch[as.factor(vegtype_MDS)])
legend(-2,-1.3, legend = Veges, col = vegcol, pch = vegpch, cex = 0.8, pt.bg = vegcol)

#ordispider(beta_mds, substr(rownames(beta_mds$points),6,6))
#orditorp(beta_mds,display="sites",cex=0.5,air=0.01)


#dev.off()


beta_disp <- betadisper(MC_dist, substr(rownames(beta_mds$points),6,6))
anova(beta_disp)
TukeyHSD(beta_disp)

test_NPMANOVA <- adonis2(MC_dist ~ as.factor(metadata$Type[-c(2,12,23,14,28,31,33,34,35,36,37)]))
test_NPMANOVA












#Rarity

resM_prop <- sweep(resM, 1, rowSums(resM_fin), "/")*100 #Makes a table of proportions. 


#Occupancy
#Ian had low occupancy being present in no more than 5% of sites (4 out of 81), so if I use a similar metric and round up I get 2 sites. 
frequency <- colSums(resM_prop > 0)

pine_freq <- colSums(resM_prop[pinesites, ] > 0)
pine_freq <- pine_freq[pine_freq > 0] #removes ones not present in the pine sites

notho_freq <- colSums(resM_prop[nothosites, ] > 0)
notho_freq <- notho_freq[notho_freq > 0]

shrub_freq <- colSums(resM_prop[shrubsites, ] > 0)
shrub_freq <- shrub_freq[shrub_freq > 0] 

frequency

#Dominance
#low dominance = less than 1% of all reads maximum across all occurrences


dominance <- colSums(resM_prop)/frequency

pine_dom <- colSums(resM_prop[pinesites,])/frequency
pine_dom <- pine_dom[pine_dom > 0]


notho_dom <- colSums(resM_prop[nothosites,])/frequency
notho_dom <- notho_dom[notho_dom > 0]

shrub_dom <- colSums(resM_prop[shrubsites,])/frequency
shrub_dom <- shrub_dom[shrub_dom > 0]

dominance


#how to do some stats on this? By taxa or by site?
dom_freq <- lm(dominance ~ frequency)
summary(dom_freq) #Maybe needs revisting

cor.test(dominance, frequency, formula = "pearsons")

dom_freq_p <- lm(pine_dom ~ pine_freq)
summary(dom_freq_p)





#ASVs unique to vegtypes


resM_mushed_pine <- as.matrix(colSums(resM[pinesites,]))
resM_mushed_notho <- as.matrix(colSums(resM[nothosites,]))
resM_mushed_shrub <- as.matrix(colSums(resM[shrubsites,]))

resM_mushed <- cbind(colnames(resM), resM_mushed_notho, resM_mushed_shrub, resM_mushed_pine)
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

ncol(resM)
c(length(pine_uniques), length(notho_uniques), length(shrub_uniques))
#MAKE A VENN DIAGRAM



shared_OTUs
OTUs_unique_to_veg[1:3]

# create Venn diagram and display all sets 



library(VennDiagram)



pine_freq <- colSums(resM[pinesites, ] > 0)
pine_freq <- pine_freq[pine_freq > 0] #removes ones not present in the pine sites
pine_OTUs <- names(pine_freq)

notho_freq <- colSums(resM[nothosites, ] > 0)
notho_freq <- notho_freq[notho_freq > 0]
notho_OTUs <- names(notho_freq)

shrub_freq <- colSums(resM[shrubsites, ] > 0)
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
total_OTUs <- c(length(notho_OTUs),length(shrub_OTUs), length(pine_OTUs), ncol(resM))
low_occupancy <- c(length(notho_freq[notho_freq < 3]),length(shrub_freq[shrub_freq < 3]),length(pine_freq[pine_freq < 3]), length(frequency[frequency < 3]))
low_dominance <- c(length(notho_dom[notho_dom < 1]), length(shrub_dom[shrub_dom < 1]), length(pine_dom[pine_dom < 1]), length(dominance))

rarity_table <- data.frame(sites, total_OTUs, low_occupancy, low_dominance, OTUs_unique_to_veg, shared_OTUs)
rarity_table <- t(rarity_table)
colnames(rarity_table) <- colnames_for_table

rarity_table

