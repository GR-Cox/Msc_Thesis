#MC success rate stats

#oo <- options(repos = "https://cran.r-project.org/")
#install.packages("Matrix")
#install.packages("lme4")
#options(oo)
library(lme4)
library(glmmTMB)
library(boot)

setwd("C:/Users/darth/Workdir/MC_stats/")

OTUs <- read.csv("Seqs_to_OTUs.csv")


seqs <- read.csv("list_of_seqs.csv", header = FALSE)

seqs <- as.vector(seqs$V1)



mc_data <- read.csv("mc_sheet.csv", row.names = 1)

t1 <- mc_data[match(seqs,mc_data$Slime.mould.1.ID), ]
t2 <- unique(mc_data[match(seqs,mc_data$Slime.mould.2.ID), ])


sites_with_seqs1 <- as.data.frame(rownames(t1), t1$Slime.mould.1.ID)
sites_with_OTUs1 <- as.data.frame(sites_with_seqs1, OTUs$OTU[match(t1$Slime.mould.1.ID, OTUs$SEQ)])

sites_with_OTUs1[rownames(sites_with_OTUs1) == "GC-OTU3",]



mc_data <- mc_data[-c(1:33),] #removes the rubbish original results for these tests as they were done a lil differently

vegtype <- substr(mc_data$Vegetation, 1, 1)
vegtype
vegtype <- factor(vegtype, levels=c("S", "N", "P"))
vegtype




mc_data$Slime.mould.Presence[mc_data$Slime.mould.Presence == "no"] <- 0
mc_data$Slime.mould.Presence[mc_data$Slime.mould.Presence == "yes"] <- 1

mc_data$Slime.mould.Presence <- as.numeric(mc_data$Slime.mould.Presence)
sum(mc_data$Slime.mould.Presence)
sum(mc_data$Number.of.Slime.moulds)
nrow(mc_data)

mc_data$Vegetation <- as.factor(mc_data$Vegetation)
mc_data$Season <- as.factor(mc_data$Season)
mc_data$Substrate <- as.factor(mc_data$Substrate)
mc_data$Site <- as.factor(mc_data$Site)

unique(mc_data$Slime.mould.2.life.cycle.stage)
plas <- sum(length(mc_data$Slime.mould.1.life.cycle.stage[mc_data$Slime.mould.1.life.cycle.stage == "plasmodium"]) + length(mc_data$Slime.mould.2.life.cycle.stage[mc_data$Slime.mould.2.life.cycle.stage == "plasmodium"]) + length(mc_data$Slime.mould.3.life.cycle.stage[mc_data$Slime.mould.3.life.cycle.stage == "plasmodium"]))
fruit <- sum(length(mc_data$Slime.mould.1.life.cycle.stage[mc_data$Slime.mould.1.life.cycle.stage == "fruiting body"]) + length(mc_data$Slime.mould.2.life.cycle.stage[mc_data$Slime.mould.2.life.cycle.stage == "fruiting body"])+ length(mc_data$Slime.mould.3.life.cycle.stage[mc_data$Slime.mould.3.life.cycle.stage == "fruiting body"]))
trail <- sum(length(mc_data$Slime.mould.1.life.cycle.stage[mc_data$Slime.mould.1.life.cycle.stage == "trail"]) + length(mc_data$Slime.mould.2.life.cycle.stage[mc_data$Slime.mould.2.life.cycle.stage == "trail"])+ length(mc_data$Slime.mould.3.life.cycle.stage[mc_data$Slime.mould.3.life.cycle.stage == "trail"]))
scler <- sum(length(mc_data$Slime.mould.1.life.cycle.stage[mc_data$Slime.mould.1.life.cycle.stage == "sclerotia"]) + length(mc_data$Slime.mould.2.life.cycle.stage[mc_data$Slime.mould.2.life.cycle.stage == "sclerotia"])+ length(mc_data$Slime.mould.3.life.cycle.stage[mc_data$Slime.mould.3.life.cycle.stage == "sclerotia"]))

plas
fruit
trail
scler
sum(plas,
    fruit,
    trail,
    scler)

pos_rate <- sum(as.numeric(mc_data$Slime.mould.Presence))/length(mc_data$Number.of.Slime.moulds)
pos_rate

pine_pos <- sum(as.numeric(mc_data$Slime.mould.Presence[vegtype == "P"]))/length(mc_data$Number.of.Slime.moulds[vegtype == "P"])
shrub_pos <- sum(as.numeric(mc_data$Slime.mould.Presence[vegtype == "S"]))/length(mc_data$Number.of.Slime.moulds[vegtype == "S"])
notho_pos <- sum(as.numeric(mc_data$Slime.mould.Presence[vegtype == "N"]))/length(mc_data$Number.of.Slime.moulds[vegtype == "N"])

pos_df <- c(pine_pos, shrub_pos, notho_pos)


pine_df <- mc_data[mc_data$Vegetation == "Pine",]
beech_df <- mc_data[mc_data$Vegetation == "Nothofagus",]
shrub_df <- mc_data[mc_data$Vegetation == "Shrubbery",]


summer <- mc_data[mc_data$Season == "Summer",]
autumn <- mc_data[mc_data$Season == "Autumn",]
spring <- mc_data[mc_data$Season == "Spring",]

summer_pos <- sum(summer$Slime.mould.Presence)/length(summer$Slime.mould.Presence)
autumn_pos <- sum(autumn$Slime.mould.Presence)/length(autumn$Slime.mould.Presence)
spring_pos <- sum(spring$Slime.mould.Presence)/length(spring$Slime.mould.Presence)

pinesum <- pine_df[pine_df$Season == "Summer",]
beechsum <- beech_df[beech_df$Season == "Summer",]
shrubsum <- shrub_df[shrub_df$Season == "Summer",]

pineaut <- pine_df[pine_df$Season == "Autumn",]
beechaut <- beech_df[beech_df$Season == "Autumn",]
shrubaut <- shrub_df[shrub_df$Season == "Autumn",]

pinespr <- pine_df[pine_df$Season == "Spring",]
beechspr <- beech_df[beech_df$Season == "Spring",]
shrubspr <- shrub_df[shrub_df$Season == "Spring",]


pinesum_pos <- sum(pinesum$Slime.mould.Presence/length(pinesum$Slime.mould.Presence))
beechsum_pos <- sum(beechsum$Slime.mould.Presence/length(beechsum$Slime.mould.Presence))
shrubsum_pos <- sum(shrubsum$Slime.mould.Presence/length(shrubsum$Slime.mould.Presence))
pineaut_pos <- sum(pineaut$Slime.mould.Presence/length(pineaut$Slime.mould.Presence))
beechaut_pos <- sum(beechaut$Slime.mould.Presence/length(beechaut$Slime.mould.Presence))
shrubaut_pos <- sum(shrubaut$Slime.mould.Presence/length(shrubaut$Slime.mould.Presence))
pinespr_pos <- sum(pinespr$Slime.mould.Presence/length(pinespr$Slime.mould.Presence))
beechspr_pos <- sum(beechspr$Slime.mould.Presence/length(beechspr$Slime.mould.Presence))
shrubspr_pos <- sum(shrubspr$Slime.mould.Presence/length(shrubspr$Slime.mould.Presence))


pinesum_se <- sd(pinesum$Slime.mould.Presence)/sqrt(length(pinesum$Slime.mould.Presence))
beechsum_se <- sd(beechsum$Slime.mould.Presence)/sqrt(length(beechsum$Slime.mould.Presence))
shrubsum_se <- sd(shrubsum$Slime.mould.Presence)/sqrt(length(shrubsum$Slime.mould.Presence))
pineaut_se <- sd(pineaut$Slime.mould.Presence)/sqrt(length(pineaut$Slime.mould.Presence))
beechaut_se <- sd(beechaut$Slime.mould.Presence)/sqrt(length(beechaut$Slime.mould.Presence))
shrubaut_se <- sd(shrubaut$Slime.mould.Presence)/sqrt(length(shrubaut$Slime.mould.Presence))
pinespr_se <- sd(pinespr$Slime.mould.Presence)/sqrt(length(pinespr$Slime.mould.Presence))
beechspr_se <- sd(beechspr$Slime.mould.Presence)/sqrt(length(beechspr$Slime.mould.Presence))
shrubspr_se <- sd(shrubspr$Slime.mould.Presence)/sqrt(length(shrubspr$Slime.mould.Presence))




veg_season_se <- c(beechsum_se, shrubsum_se, pinesum_se, beechaut_se, shrubaut_se, pineaut_se, beechspr_se, shrubspr_se, pinespr_se)
#Fix standard errors for this.

veg_season <- c(beechsum_pos, shrubsum_pos, pinesum_pos, beechaut_pos, shrubaut_pos, pineaut_pos, beechspr_pos, shrubspr_pos, pinespr_pos)
veg_seas_names <- c("Beech-Summer", "Shrubland-Summer", "Pine-Summer", "Beech-Autumn", "Shrubland-Autumn", "Pine-Autumn", "Beech-Spring", "Shrubland-Spring", "Pine-Spring")




test <- glmmTMB(as.numeric(Slime.mould.Presence) ~ Vegetation*Season+Substrate + 
                (1|Site), family = binomial, data = mc_data, REML=FALSE) #REML = FALSE because you need to use maximum lieklihood to test models
              
MASS::stepAIC(test)


glm1 <- glmmTMB(Slime.mould.Presence ~ Vegetation*Season+Substrate + 
                  (1|Site), family = binomial, data = mc_data, REML=TRUE) #reml = TRE because its better, I think


summary(glm1)
summary(aov(glm1))

anova(glm1, )


glmspring <- glmmTMB(Slime.mould.Presence ~ Vegetation + 
                  (1|Site), family = binomial, data = mc_data[mc_data$Season == "Spring",], REML=TRUE)


summary(glmspring)



glmsummer <- glmmTMB(Slime.mould.Presence ~ Vegetation + 
                       (1|Site), family = binomial, data = mc_data[mc_data$Season == "Summer",], REML=TRUE)


summary(glmsummer)

glmautumn <- glmmTMB(Slime.mould.Presence ~ Vegetation + 
                       (1|Site), family = binomial, data = mc_data[mc_data$Season == "Autumn",], REML=TRUE)


summary(glmautumn)


glm_substrate <- glmmTMB(Slime.mould.Presence ~ Substrate + (1|Site), family = binomial, data = mc_data, REML=TRUE)

summary(glm_substrate)

glm_pine <- glmmTMB(Slime.mould.Presence ~ Season + 
                      (1|Site), family = binomial, data = mc_data[mc_data$Vegetation== "Pine",], REML=TRUE)

summary(glm_pine)

glm_beech <- glmmTMB(Slime.mould.Presence ~ Season + 
                      (1|Site), family = binomial, data = mc_data[mc_data$Vegetation== "Nothofagus",], REML=TRUE)

summary(glm_beech)

glm_shrub <- glmmTMB(Slime.mould.Presence ~ Season + 
                      (1|Site), family = binomial, data = mc_data[mc_data$Vegetation== "Shrubbery",], REML=TRUE)

summary(glm_shrub)


#Get paramater estimates and std.errors from model then back transform, use for plots, should fix standard. 

#Beech in summer

beech_sum_logit <- -1.3358 # These values come from the sumumn GLM and match the 
beech_sum_se_logit <- 0.3393

beech_sum_plot_logit <- c((beech_sum_logit+beech_sum_se_logit), beech_sum_logit, (beech_sum_logit-beech_sum_se_logit))

inv.logit(beech_sum_plot_logit)
beechsum_pos

#Shrub in summer

shrub_sum_logit <- -0.9866+-1.3358 # These values come from the sumumn GLM and match the 
shrub_sum_se_logit <- 0.5433

shrub_sum_plot_logit <- c((shrub_sum_logit+shrub_sum_se_logit), shrub_sum_logit, (shrub_sum_logit-shrub_sum_se_logit))

inv.logit(shrub_sum_plot_logit)
shrubsum_pos

#pine in summer

pine_sum_logit <- -0.6618+-1.3358 # These values come from the sumumn GLM and match the 
pine_sum_se_logit <- 0.5140

pine_sum_plot_logit <- c((pine_sum_logit+pine_sum_se_logit), pine_sum_logit, (pine_sum_logit-pine_sum_se_logit))

inv.logit(pine_sum_plot_logit)
pinesum_pos

#Beech in autumn

beech_aut_logit <- -1.846e+00 # These values come from the autumn GLM and match the 
beech_aut_se_logit <- 3.106e-01

beech_aut_plot_logit <- c((beech_aut_logit+beech_aut_se_logit), beech_aut_logit, (beech_aut_logit-beech_aut_se_logit))

inv.logit(beech_aut_plot_logit)
beechaut_pos

#Shrub in Autumn

shrub_aut_logit <- 4.164e-01+-1.846e+00 # These values come from the autumn GLM and match the 
shrub_aut_se_logit <- 4.116e-01

shrub_aut_plot_logit <- c((shrub_aut_logit+shrub_aut_se_logit), shrub_aut_logit, (shrub_aut_logit-shrub_aut_se_logit))

inv.logit(shrub_aut_plot_logit)
shrubaut_pos

#Pine in Autumn

pine_aut_logit <- -1.005e-16+-1.846e+00 # These values come from the autumn GLM and match the 
pine_aut_se_logit <- 4.393e-01

pine_aut_plot_logit <- c((pine_aut_logit+pine_aut_se_logit), pine_aut_logit, (pine_aut_logit-pine_aut_se_logit))

inv.logit(pine_aut_plot_logit)
pineaut_pos

#beech in Spring

beech_spr_logit <- -0.9447 # These values come from the sprumn GLM and match the 
beech_spr_se_logit <- 0.3039

beech_spr_plot_logit <- c((beech_spr_logit+beech_spr_se_logit), beech_spr_logit, (beech_spr_logit-beech_spr_se_logit))

inv.logit(beech_spr_plot_logit)
beechspr_pos


#shrub in Spring

shrub_spr_logit <- -0.3030+-0.9447 # These values come from the sprumn GLM and match the 
shrub_spr_se_logit <- 0.4390

shrub_spr_plot_logit <- c((shrub_spr_logit+shrub_spr_se_logit), shrub_spr_logit, (shrub_spr_logit-shrub_spr_se_logit))

inv.logit(shrub_spr_plot_logit)
shrubspr_pos

#Pine in Spring

pine_spr_logit <- -0.1352+-0.9447 # These values come from the sprumn GLM and match the 
pine_spr_se_logit <- 0.4344

pine_spr_plot_logit <- c((pine_spr_logit+pine_spr_se_logit), pine_spr_logit, (pine_spr_logit-pine_spr_se_logit))

inv.logit(pine_spr_plot_logit)
pinespr_pos


summer_pos <- mean(inv.logit(mean(c(beech_sum_logit, shrub_sum_logit, pine_sum_logit))))
summer_pos

autumn_pos <- mean(inv.logit(mean(c(beech_aut_logit, shrub_aut_logit, pine_aut_logit))))
autumn_pos

spring_pos <- mean(inv.logit(mean(c(beech_spr_logit, shrub_spr_logit, pine_spr_logit))))
spring_pos




plot_means_logit <- c(beech_sum_logit, shrub_sum_logit, pine_sum_logit, beech_aut_logit, shrub_aut_logit, pine_aut_logit, beech_spr_logit, shrub_spr_logit, pine_spr_logit)

plot_means <- inv.logit(plot_means_logit)

plot_se_logit <- c(beech_sum_se_logit, shrub_sum_se_logit, pine_sum_se_logit, beech_aut_se_logit, shrub_aut_se_logit, pine_aut_se_logit, beech_spr_se_logit, shrub_spr_se_logit, pine_spr_se_logit)


pdf(file = "MC Success Plot Main.pdf", width = 16.5, height = 11.7)
par(mfrow=c(1,1))
par(mar = c(5, 4, 4, 8),                                 
    xpd = TRUE)

glm_cols <- c("#F0E44280", "#F0E44280", "#F0E44280", "#0072B280", "#0072B280", "#0072B280","#009E7380", "#009E7380", "#009E7380")
glm_names <- c("Beech", "Shrubland", "Pine", "Beech", "Shrubland", "Pine", "Beech", "Shrubland", "Pine")

glm_barplot <- barplot(plot_means, ylim = c(0, 0.35), col = glm_cols, names.arg = glm_names, cex.names = 1.5,
                       main = "Proportion of Successful Moist Chambers in Beech, Shrubland, and Pine across Three Seasons", cex.axis = 1.5)
arrows(x0 = glm_barplot,
       y0 = inv.logit(plot_means_logit + plot_se_logit),
       y1 = inv.logit(plot_means_logit - plot_se_logit),
       angle = 90, code = 3, length = 0.1)
legend('topright', legend = c("Summer", "Autumn", "Spring"), pt.bg = c("#F0E44280", "#0072B280", "#009E7380"), pch = 21, cex = 2, pt.cex = 2.5, inset = c(-0.1, 0))

dev.off()

#Plot should have each sites success rate as the point, then the inv logit se and mean as a bar.


#litter and bark

litter <- mc_data[mc_data$Substrate == "Litter",]
bark <- mc_data[mc_data$Substrate == "Bark",]

lit_pos <- sum(litter$Slime.mould.Presence)/nrow(litter)
bark_pos <- sum(bark$Slime.mould.Presence)/nrow(bark)


litterlogit <- 0.2572+-1.6387
litterlogit_se <- 0.1850

barklogit <- -1.6387
barklogit_se <- 0.1592

sub_means <- c(litterlogit, barklogit)
sub_se <- c(litterlogit_se, barklogit_se)

pdf(file = "MC Success Plot Substrate.pdf", width = 8.25, height = 6.5)

par(mfrow=c(1,1))
sub_barplot <- barplot(inv.logit(sub_means), ylim = c(0, 0.3), col = c("#56B4E980", "#D55E0080"), names.arg = c("Litter", "Bark"), cex.axis = 2, cex.names = 2, main = "Proportion of Successful Moist Chambers from each Substrate")
arrows(x0 = sub_barplot,
       y0 = inv.logit(sub_means + sub_se),
       y1 = inv.logit(sub_means - sub_se),
       angle = 90, code = 3, length = 0.1)

dev.off()

#Sequencing success


collections <- c(unique(mc_data$Slime.mould.1.ID),unique(mc_data$Slime.mould.2.ID), unique(mc_data$Slime.mould.3.ID))
collections <- collections[-c(1, 149, 152)]
collections

length(collections)
length(seqs)
base::setdiff(collections, seqs)


