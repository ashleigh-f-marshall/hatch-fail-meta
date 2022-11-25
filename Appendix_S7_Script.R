####################
# Authors: Ashleigh Fleming Marshall
# ashleigh.marshall@ioz.ac.uk | ashleigh.marshall.16@ucl.ac.uk
####################

### SYSTEMATIC REVIEW OF AVIAN HATCHING FAILURE AND IMPLICATIONS FOR CONSERVATION

##Clear workspace
rm(list=ls())

#### INSTALL AND LOAD PACKAGES ####
library("ape")
library("dplyr")
library("forcats")
library("forestplot")
library("ggplot2")
library("metafor")
library("phangorn")
library("phytools")

#### LOAD AND PREPARE DATA ####

## Read in data
data <- read.csv("Data/REVISED_Appendix_S8_Dataset.csv")
dim(data) #483
head(data$Study.Name.numbered)
tail(data$Study.Name.numbered)

## Rename rows from Management Level column
data$Wild.Managed.Captive<-gsub("Wild managed","ManagedWild",data$Wild.Managed.Captive)

#### COMPUTE PHYLOGENETIC RELATEDNESS CORRELATION MATRIX ####
## Use 1,000 trees downloaded from www.birdtree.org based on our dataset, compute the maximum clade credibility tree, compute branch lengths, compute the correlation matrix

## Remove records of species not available on BirdTree
data_final <- data[-which(is.na(data$BirdTree.Tip)),]
final_species <- sort(unique(as.character(data_final$BirdTree.Tip)))
tips <- data_final$BirdTree.Tip

## Load in tree downloaded from BirdTree
trees <- read.nexus("Data/BirdTree-DownloadedTrees/output.nex",tree.names=T)

pruned.trees <- lapply(trees,keep.tip,tip=tips)
class(pruned.trees) <- "multiPhylo" 
mcc_pruned <- maxCladeCred(pruned.trees)
mcc_pruned <- ladderize(mcc_pruned)
mcc_pruned$tip.label #241
mcc_ult <- compute.brlen(mcc_pruned, power = 1)
plot.phylo(mcc_ult, cex=0.2) #Figure S1
bird_phylo_cor <- vcv(mcc_ult, cor=T)

summary(bird_phylo_cor[row(bird_phylo_cor)!=col(bird_phylo_cor)])

## Add additional columns
data_final$phylo<-data_final$BirdTree.Tip
data_final$effect_size_id<-seq(nrow(data_final))  

## Tree for just wild and wild managed populations
data_final_no.captive <- data_final[data_final$Wild.Managed.Captive!="Captive",]
tips_no.captive <- data_final_no.captive$BirdTree.Tip
pruned.trees_no.captive <- lapply(trees,keep.tip,tip=tips_no.captive)
class(pruned.trees_no.captive)<-"multiPhylo"
mcc_pruned_no.captive <- maxCladeCred(pruned.trees_no.captive)
mcc_pruned_no.captive <- ladderize(mcc_pruned_no.captive)
mcc_pruned_no.captive$tip.label #233 species
mcc_multi_no.captive <- mcc_pruned_no.captive$tip.label
mcc_ult_no.captive <- compute.brlen(mcc_pruned_no.captive, power = 1)
bird_phylo_cor_no.captive <- vcv(mcc_ult_no.captive, cor=T) 

ies.da=escalc(xi=Total-Cases, ni=Total, data=data_final, measure="PFT", add=0)

#### CHECK PHYLOGENETIC SIGNAL ####

hatchfailure.es <- ies.da$yi
propHF <- ((ies.da$Total-ies.da$Cases)/(ies.da$Total))
names(hatchfailure.es) <- ies.da$BirdTree.Tip
names(propHF) <- ies.da$BirdTree.Tip

phylosig(mcc_ult,hatchfailure.es, method="K", test=TRUE, nsim=1000)
phylosig(mcc_ult,hatchfailure.es, method="lambda", test=TRUE)
phylosig(mcc_ult,propHF, method="K", test=TRUE, nsim=1000)
phylosig(mcc_ult,propHF, method="lambda", test=TRUE)

#### CHECK FOR MULTI-COLLINEARITY ####
## Testing for multi-collinearity of variables being used for moderators using pairwise Chi-square Tests of Independence and through calculation of VIF and GVIG values

## Chi-square Tests (significance and Cramer's V values) - Table S11

## Threat Status vs. Management Level
tbl_1<-table(data_final$Threat_Status_SUB,data_final$Wild.Managed.Captive)
chi2_1 <- chisq.test(tbl_1,correct=F)
c(chi2_1$statistic,chi2_1$p.value)
sqrt(chi2_1$statistic/sum(tbl_1)) #Cramer's V
## Threat Status vs. Incubation
tbl_2 = table(data_final$Threat_Status_SUB,data_final$Incubation)
chi2_2 = chisq.test(tbl_2,correct=F)
c(chi2_2$statistic,chi2_2$p.value)
sqrt(chi2_2$statistic/sum(tbl_2))
## Threat Status vs. Supplementary Feeding
tbl_3 = table(data_final$Threat_Status_SUB,data_final$Supplementary.Feeding)
chi2_3 = chisq.test(tbl_3,correct=F)
c(chi2_3$statistic,chi2_3$p.value)
sqrt(chi2_3$statistic/sum(tbl_3))
## Threat Status vs. Artificial Nest Provision
tbl_4 = table(data_final$Threat_Status_SUB,data_final$Nestboxes.Burrows.Sites)
chi2_4 = chisq.test(tbl_4,correct=F)
c(chi2_4$statistic,chi2_4$p.value)
sqrt(chi2_4$statistic/sum(tbl_4))
## Management Level vs. Incubation
tbl_5 = table(data_final$Wild.Managed.Captive,data_final$Incubation)
chi2_5 = chisq.test(tbl_5,correct=F)
c(chi2_5$statistic,chi2_5$p.value)
sqrt(chi2_5$statistic/sum(tbl_5))
## Management Level vs. Supplementary Feeding
tbl_6 = table(data_final$Wild.Managed.Captive,data_final$Supplementary.Feeding)
chi2_6 = chisq.test(tbl_6,correct=F)
c(chi2_6$statistic,chi2_6$p.value)
sqrt(chi2_6$statistic/sum(tbl_6))
## Management Level vs. Artificial Nest Provision
tbl_7 = table(data_final$Wild.Managed.Captive,data_final$Nestboxes.Burrows.Sites)
chi2_7 = chisq.test(tbl_7,correct=F)
c(chi2_7$statistic,chi2_7$p.value)
sqrt(chi2_7$statistic/sum(tbl_7))
## Incubation vs. Supplementary Feeding
tbl_8 = table(data_final$Incubation,data_final$Supplementary.Feeding)
chi2_8 = chisq.test(tbl_8,correct=F)
c(chi2_8$statistic,chi2_8$p.value)
sqrt(chi2_8$statistic/sum(tbl_8))
## Incubation vs. Artificial Nest Provision
tbl_9 = table(data_final$Incubation,data_final$Nestboxes.Burrows.Sites)
chi2_9 = chisq.test(tbl_9,correct=F)
c(chi2_9$statistic,chi2_9$p.value)
sqrt(chi2_9$statistic/sum(tbl_9))
## Supplementary Feeding vs. Artificial Nest Provision
tbl_10 = table(data_final$Supplementary.Feeding,data_final$Nestboxes.Burrows.Sites)
chi2_10 = chisq.test(tbl_10,correct=F)
c(chi2_10$statistic,chi2_10$p.value)
sqrt(chi2_10$statistic/sum(tbl_10))

## Variance Inflation Factors (VIF) and Generalised Variance Inflation Factors (GVIF) - Table S12
mods_test<-rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Threat_Status_SUB)+factor(Wild.Managed.Captive)+factor(Nestboxes.Burrows.Sites)+factor(Incubation)+factor(Supplementary.Feeding), data=ies.da)
vif.rma(mods_test)
vif.rma(mods_test,table=TRUE)
vif.rma(mods_test,btt=3:4) #GVIF for management level = 12.0335
(12.0335^(1/(2*2)))^2

#### CALCULATING OVERALL EFFECT AND FITTING MODELS ####

## Fit a hierarchical three-level (random effects for effect size and study) intercept-only model to calculate overall effect
overall <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name),test="t",method="REML",data=ies.da)
summary(overall)
predict(overall, transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) 

## Fit a two-level model without within-study variance
overall_novar2 <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name),sigma2=c(0,NA),test="t",method="REML",data=ies.da)
summary(overall_novar2)
predict(overall_novar2, transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) 

## Fit a two-level model without between-study variance
overall_novar3 <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name),
sigma2=c(NA,0),test="t",method="REML",data=ies.da)
summary(overall_novar3)
predict(overall_novar3, transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Use ANOVA to test the fit of the models
anova(overall,overall_novar2) #overall better than overall_novar2
anova(overall,overall_novar3) #overall better than overall_novar3

## There is significant within-study variance and between-study variance - therefore can apply moderator analyses

## Fit a multilevel meta-analytic model (random effects for effect size, study, species, and phylogeny)
overall_phylo_sp <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),test="t",method="REML",data=ies.da)
summary(overall_phylo_sp)
predict(overall_phylo_sp, transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Use ANOVA to test the fit of the models
anova(overall,overall_phylo_sp) #The model with phylogeny and species included is a better fit

## Use Profile Likelihood Plots to check the variance components are identifiable - Figure S2
par(mfrow=c(4,1))
sigma2=1
profile(overall_phylo_sp,sigma2=1,ylab="Restricted log-likelihood",xlab=bquote(sigma[.(sigma2)]^2 ~ "value"),main=bquote("Profile plot for" ~ sigma[.(sigma2)]^2))
sigma2=2
profile(overall_phylo_sp,sigma2=2,ylab="Restricted log-likelihood",xlab=bquote(sigma[.(sigma2)]^2 ~ "value"),main=bquote("Profile plot for" ~ sigma[.(sigma2)]^2))
sigma2=3
profile(overall_phylo_sp,sigma2=3,ylab="Restricted log-likelihood",xlab=bquote(sigma[.(sigma2)]^2 ~ "value"),main=bquote("Profile plot for" ~ sigma[.(sigma2)]^2))
sigma2=4
profile(overall_phylo_sp,sigma2=4,ylab="Restricted log-likelihood",xlab=bquote(sigma[.(sigma2)]^2 ~ "value"),main=bquote("Profile plot for" ~ sigma[.(sigma2)]^2))
par(mfrow=c(1,1))

#### PARTITIONING HETEROGENEITY ####

W <- diag(1/ies.da$vi)
X <- model.matrix(overall_phylo_sp)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

## Total I^2
100 * sum(overall_phylo_sp$sigma2) / (sum(overall_phylo_sp$sigma2) + (overall_phylo_sp$k-overall_phylo_sp$p)/sum(diag(P)))
## I^2 of each level
100 * overall_phylo_sp$sigma2 / (sum(overall_phylo_sp$sigma2) + (overall_phylo_sp$k-overall_phylo_sp$p)/sum(diag(P))) # Table S9

#### IDENTIFYING OUTLIERS ####
## Identifying potential influential outliers by calculating Cook's Distance - Table S7
cooks<-cooks.distance.rma.mv(overall_phylo_sp, progbar=TRUE, reestimate = FALSE)
plot(cooks,type="o",pch=19,xlab="Observed Outcome",ylab="Cook's Distance")
cooks[(cooks>=0.5)]
cooks[(cooks>(3*mean(cooks,na.rn=TRUE)))] #25, 201, 207, 264, 291, 306, 314, 366, 401, 402, 404
cooks[(cooks>(4/478))] #207, 366, 401

## Identifying potential influential outliers by calculating the externally studentized residuals - Table S7
options(max.print=10000)
rstud<-rstudent(overall_phylo_sp, progbar=TRUE,
                #reestimate=FALSE #Takes time elapsed from ~2.5 hrs to <10 mins but only yields estimation (very similar)
                )
abs.z_mv=abs(rstud$z)
rstud[order(abs.z_mv)] #[>1.96 = 197, 214, 13, 366, 258, 401, 259, 76, 471, 201, 472, 187]

all_outliers<-data_final[c(13,25,76,187,197,201,207,214,258,259,264,291,306,314,366,401,402,404,471,472),]
all_outliers$Study.Name.numbered

#### MIXED EFFECTS MULTILEVEL META-ANALYTICAL MODELS ####
## Applying mixed effects multilevel meta-analytical models to conduct moderator analyses; potential moderators = threat status (threatened or non-threatened), Red List classification (LC, NT, VU, EN, CR), management level (wild, wild managed, or captive), incubation type (natural or artificial), supplementary feeding (fed or not fed), artificial nest provision (provided or not provided), and publication year.

## Threat Status - Table S8
overall_phylo_sp_threat <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~relevel(factor(Threat_Status_SUB),ref="Threatened"), test="t", method="REML", data=ies.da)
summary(overall_phylo_sp_threat)
predict(overall_phylo_sp_threat,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Red List - Table S8
overall_phylo_sp_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~relevel(factor(IUCN.Status.Current_SUB),ref="VU"),test="t",method="REML", data=ies.da)
summary(overall_phylo_sp_redlist)
predict(overall_phylo_sp_redlist,newmods=rbind(c(0,0,0,0), c(1,0,0,0), c(0,1,0,0),c(0,0,1,0),c(0,0,0,1)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Management Level - Table S8
overall_phylo_sp_manage <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~relevel(factor(Wild.Managed.Captive),ref="Wild"), test="t", method="REML", data=ies.da)
summary(overall_phylo_sp_manage)
predict(overall_phylo_sp_manage,newmods=rbind(c(0,0), c(1,0), c(0,1)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Incubation Type - Table S8
overall_phylo_sp_incubate <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor),mods = ~relevel(factor(Incubation), ref="Natural"),test="t",method="REML",data=ies.da)
summary(overall_phylo_sp_incubate)
predict(overall_phylo_sp_incubate,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Incubation Type - no captive - Table S23
ies.da_no.captive=escalc(xi=Total-Cases, ni=Total, data=data_final_no.captive, measure="PFT", add=0) 

overall_phylo_sp_incubate_no.captive <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_no.captive),mods = ~relevel(factor(Incubation),ref="Natural"), test="t",method="REML", data=ies.da_no.captive)
summary(overall_phylo_sp_incubate_no.captive)
predict(overall_phylo_sp_incubate_no.captive,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da_no.captive$Total))

## Supplementary Feeding - Table S8
overall_phylo_sp_suppfed <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~relevel(factor(Supplementary.Feeding),ref="No"), test="t", method="REML",data=ies.da)
summary(overall_phylo_sp_suppfed)
predict(overall_phylo_sp_suppfed,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Supplementary Feeding - no captive - Table S23
overall_phylo_sp_suppfed_no.captive <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_no.captive),mods = ~relevel(factor(Supplementary.Feeding),ref="No"), test="t", method="REML", data=ies.da_no.captive)
summary(overall_phylo_sp_suppfed_no.captive)
predict(overall_phylo_sp_suppfed_no.captive,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da_no.captive$Total))

## Artificial Nest Provision - Table S8
overall_phylo_sp_nest <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~relevel(factor(Nestboxes.Burrows.Sites),ref="No"), test="t", method="REML", data=ies.da)
summary(overall_phylo_sp_nest)
predict(overall_phylo_sp_nest,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Artificial Nest Provision - no captive - Table S23
overall_phylo_sp_nest_no.captive <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_no.captive), mods = ~relevel(factor(Nestboxes.Burrows.Sites),ref="No"), test="t", method="REML", data=ies.da_no.captive)
summary(overall_phylo_sp_nest_no.captive)
predict(overall_phylo_sp_nest_no.captive,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da_no.captive$Total))

## Publication Year - Figure S3
overall_phylo_sp_pubyear <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~Publication.Year, test="t", method="REML", data=ies.da)
summary(overall_phylo_sp_pubyear)
preds=predict(overall_phylo_sp_pubyear,newmods=c(2004:2020),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

## Check for unexplained heterogeneity if all moderators included
overall_phylo_sp_all <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Threat_Status_SUB)+factor(Wild.Managed.Captive)+factor(Incubation)+factor(Supplementary.Feeding)+factor(Nestboxes.Burrows.Sites), test="t", method="REML", data=ies.da)
summary(overall_phylo_sp_all)

#### MODELLING INTERACTIONS####
## Multivariate mixed effects meta-analytical models used o test for evidence of interactions between threat status and Red List classification with management level, artificial incubation, supplementary feeding, and artificial nest provision.

## THREAT STATUS ####
## Threat Status + Management Level - Table 2 and Table S13
meta_additive <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Wild.Managed.Captive)+relevel(factor(Threat_Status_SUB),ref="Non-threatened"), test="t", method="REML", data=ies.da)
summary(meta_additive) #Threat Status is significant
anova(meta_additive,btt=2:3) #Management level is significant
predict(meta_additive,newmods=rbind(c(0,0,0),c(1,0,0),c(0,1,0)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #Non-threatened
predict(meta_additive,newmods=rbind(c(0,0,1),c(1,0,1),c(0,1,1)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #Threatened

meta_additive_nointercept <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Wild.Managed.Captive)+relevel(factor(Threat_Status_SUB),ref="Threatened")-1, test="t", method="REML", data=ies.da)
summary(meta_additive_nointercept)

## Threat Status * Management Level - Table S13
meta_interactive_1 <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Wild.Managed.Captive)*relevel(factor(Threat_Status_SUB),ref="Non-threatened"), test="t", method="REML", data=ies.da)
summary(meta_interactive_1)
anova(meta_interactive_1,btt=5:6) #Not significant = no interaction present

meta_additive_ml <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Wild.Managed.Captive)+factor(Threat_Status_SUB),test="t",method="ML",data=ies.da)
meta_interactive_1_ml <-rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Wild.Managed.Captive)*factor(Threat_Status_SUB),test="t",method="ML",data=ies.da)

anova(meta_additive_ml,meta_interactive_1_ml) #Not significant = no interaction present

predict(meta_interactive_1,newmods=c(0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #captive x non-threatened
predict(meta_interactive_1,newmods=c(1,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #wild managed x non-threatened
predict(meta_interactive_1,newmods=c(0,1,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #wild x non-threatened
predict(meta_interactive_1,newmods=c(0,0,1,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #captive x threatened
predict(meta_interactive_1,newmods=c(1,0,1,1,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #wild managed x threatened
predict(meta_interactive_1,newmods=c(0,1,1,0,1),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #wild x threatened

## Threat Status : Management Level
meta_interactive_2 <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Wild.Managed.Captive):factor(Threat_Status_SUB)-1, test="t", method="REML", data=ies.da)
summary(meta_interactive_2)

## Threat Status + Incubation Type - Table 2 and Table S14
meta_additive_inc <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Incubation)+relevel(factor(Threat_Status_SUB),ref="Non-threatened"), test="t", method="REML", data=ies.da)
summary(meta_additive_inc) #Threat Status and Incubation type are significant
predict(meta_additive_inc,newmods=rbind(c(0,0),c(1,0),c(0,1),c(1,1)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

meta_additive_inc_nointercept <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Incubation)+relevel(factor(Threat_Status_SUB),ref="Non-threatened")-1, test="t", method="REML", data=ies.da)
summary(meta_additive_inc_nointercept)

## Threat Status * Incubation Type - Table S14
meta_interactive_inc_1 <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Incubation)*relevel(factor(Threat_Status_SUB),ref="Non-threatened"), test="t", method="REML", data=ies.da)
summary(meta_interactive_inc_1) #Natural:Threatened = not significant --> No interactive effect
anova(meta_interactive_inc_1, btt=4) #Not significant = no interaction present

meta_additive_inc_ml <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Incubation)+factor(Threat_Status_SUB),test="t",method="ML",data=ies.da)
meta_interactive_inc_1_ml <-rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Incubation)*factor(Threat_Status_SUB),test="t",method="ML",data=ies.da)

anova(meta_additive_inc_ml,meta_interactive_inc_1_ml) #Not significant = no interaction present

predict(meta_interactive_inc_1,newmods=c(0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #artificial x non-threatened
predict(meta_interactive_inc_1,newmods=c(1,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #natural x non-threatened
predict(meta_interactive_inc_1,newmods=c(0,1,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #artificial x threatened
predict(meta_interactive_inc_1,newmods=c(1,1,1),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #natural x threatened

## Threat Status : Incubation Type
meta_interactive_inc_2 <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Incubation):relevel(factor(Threat_Status_SUB),ref="Non-threatened")-1, test="t", method="REML", data=ies.da)
summary(meta_interactive_inc_2)

## Threat Status + Supplementary Feeding - Table 2 and Table S15
meta_additive_supp <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Supplementary.Feeding)+relevel(factor(Threat_Status_SUB),ref="Non-threatened"), test="t", method="REML", data=ies.da)
summary(meta_additive_supp) #Threat Status and Supplementary Feeding are significant
predict(meta_additive_supp,newmods=rbind(c(0,0),c(1,0),c(0,1),c(1,1)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

meta_additive_supp_nointercept <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Supplementary.Feeding)+relevel(factor(Threat_Status_SUB),ref="Non-threatened")-1, test="t", method="REML", data=ies.da)
summary(meta_additive_supp_nointercept)

## Threat Status * Supplementary Feeding - Table S15
meta_interactive_supp_1 <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Supplementary.Feeding)*relevel(factor(Threat_Status_SUB),ref="Non-threatened"), test="t", method="REML", data=ies.da)
summary(meta_interactive_supp_1) #Yes:Threatened = not significant --> No interactive effect
anova(meta_interactive_supp_1, btt=4) #Not significant = no interaction present

meta_additive_supp_ml <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Supplementary.Feeding)+factor(Threat_Status_SUB),test="t",method="ML",data=ies.da)
meta_interactive_supp_1_ml <-rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Supplementary.Feeding)*factor(Threat_Status_SUB),test="t",method="ML",data=ies.da)

anova(meta_additive_supp_ml,meta_interactive_supp_1_ml) #Not significant = no interaction present

predict(meta_interactive_supp_1,newmods=c(0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #not fed x non-threatened
predict(meta_interactive_supp_1,newmods=c(1,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #fed x non-threatened
predict(meta_interactive_supp_1,newmods=c(0,1,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #not fed x threatened
predict(meta_interactive_supp_1,newmods=c(1,1,1),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #fed x threatened

#Threat Status : Supplementary Feeding
meta_interactive_supp_2 <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Supplementary.Feeding):relevel(factor(Threat_Status_SUB),ref="Non-threatened")-1, test="t", method="REML", data=ies.da)
summary(meta_interactive_supp_2)

## Threat Status + Artificial Nest Provision - Table 2 and Table S16
meta_additive_nest <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Nestboxes.Burrows.Sites)+relevel(factor(Threat_Status_SUB),ref="Non-threatened"), test="t", method="REML", data=ies.da)
summary(meta_additive_nest) #Threat Status and Artificial Nest Provision are significant
predict(meta_additive_nest,newmods=rbind(c(0,0),c(1,0),c(0,1),c(1,1)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))

meta_additive_nest_nointercept <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Nestboxes.Burrows.Sites)+relevel(factor(Threat_Status_SUB),ref="Threatened")-1, test="t", method="REML", data=ies.da)
summary(meta_additive_nest_nointercept)

## Threat Status * Artificial Nest Provision - Table S16
meta_interactive_nest_1 <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Nestboxes.Burrows.Sites)*relevel(factor(Threat_Status_SUB),ref="Non-threatened"), test="t", 
method="REML", data=ies.da)
summary(meta_interactive_nest_1) #Yes:Threatened = not significant --> No interactive effect
anova(meta_interactive_nest_1, btt=4) #Not significant = no interaction present

meta_additive_nest_ml <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Nestboxes.Burrows.Sites)+factor(Threat_Status_SUB),test="t",method="ML",data=ies.da)
meta_interactive_nest_1_ml <-rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Nestboxes.Burrows.Sites)*factor(Threat_Status_SUB),test="t",method="ML",data=ies.da)

anova(meta_additive_nest_ml,meta_interactive_nest_1_ml) #Not significant = no interaction present

predict(meta_interactive_nest_1,newmods=c(0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #not provided x non-threatened
predict(meta_interactive_nest_1,newmods=c(1,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #provided x non-threatened
predict(meta_interactive_nest_1,newmods=c(0,1,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #not provided x threatened
predict(meta_interactive_nest_1,newmods=c(1,1,1),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #provided x threatened

## Threat Status : Artificial Nest Provision
meta_interactive_nest_2 <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Nestboxes.Burrows.Sites):relevel(factor(Threat_Status_SUB),ref="Non-threatened")-1, test="t", method="REML", data=ies.da)
summary(meta_interactive_nest_2)

## RED LIST CLASSIFICATION ####
## Red List Classification + Management Level - Table 2 and Table S17
meta_additive_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Wild.Managed.Captive)+relevel(factor(IUCN.Status.Current_SUB),ref="LC"), test="t", method="REML", data=ies.da)
summary(meta_additive_redlist)
anova(meta_additive_redlist,btt=2:3) #Management level is significant
anova(meta_additive_redlist,btt=4:7) #Red List Classification is significant
predict(meta_additive_redlist,newmods=rbind(c(0,0,0,0,0,0),c(1,0,0,0,0,0),c(0,1,0,0,0,0)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total))  #Change reference level to get values for different Red List Classifications

meta_additive_nointercept_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Wild.Managed.Captive)+relevel(factor(IUCN.Status.Current_SUB),ref="LC")-1, test="t", method="REML", data=ies.da)
summary(meta_additive_nointercept_redlist)

## Red List Classification * Management Level - Table S17
meta_interactive_1_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Wild.Managed.Captive)*relevel(factor(IUCN.Status.Current_SUB),ref="CR"), test="t",method="REML", data=ies.da)
summary(meta_interactive_1_redlist)
anova(meta_interactive_1_redlist,btt=8:15) #Not significant = no interaction present

meta_additive_ml_redlist <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Wild.Managed.Captive)+factor(IUCN.Status.Current_SUB),test="t",method="ML",data=ies.da)
meta_interactive_1_ml_redlist  <-rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Wild.Managed.Captive)*factor(IUCN.Status.Current_SUB),test="t",method="ML",data=ies.da)

anova(meta_additive_ml_redlist,meta_interactive_1_ml_redlist) #Not significant = no interaction present

predict(meta_interactive_1_redlist,newmods=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #captive x reference level
predict(meta_interactive_1_redlist,newmods=c(1,0,0,0,0,0,0,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #wild managed x reference level
predict(meta_interactive_1_redlist,newmods=c(0,1,0,0,0,0,0,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #wild x reference level

## Red List Classification : Management Level
meta_interactive_2_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Wild.Managed.Captive):relevel(factor(IUCN.Status.Current_SUB),ref="CR")-1, test="t", method="REML", data=ies.da)
summary(meta_interactive_2_redlist)

## Red List Classification + Incubation Type - Table 2 and Table S18
meta_additive_inc_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Incubation)+relevel(factor(IUCN.Status.Current_SUB),ref="CR"), test="t", method="REML", data=ies.da)
summary(meta_additive_inc_redlist) #Incubation Type is significant
anova(meta_additive_inc_redlist,btt=3:6) #Red List Classification is significant
predict(meta_additive_inc_redlist,newmods=rbind(c(0,0,0,0,0),c(1,0,0,0,0)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) 

meta_additive_nointercept_inc_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Incubation)+relevel(factor(IUCN.Status.Current_SUB),ref="LC")-1, test="t", method="REML", data=ies.da)
summary(meta_additive_nointercept_inc_redlist)

## Red List Classification * Incubation Type - Table S18
meta_interactive_1_inc_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Incubation)*relevel(factor(IUCN.Status.Current_SUB),ref="LC"), test="t", method="REML", data=ies.da)
summary(meta_interactive_1_inc_redlist)
anova(meta_interactive_1_inc_redlist,btt=7:10) #Not significant = no interaction present

meta_additive_ml_inc_redlist <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Incubation)+factor(IUCN.Status.Current_SUB),test="t",method="ML",data=ies.da)
meta_interactive_1_ml_inc_redlist  <-rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Incubation)*factor(IUCN.Status.Current_SUB),test="t",method="ML",data=ies.da)

anova(meta_additive_ml_inc_redlist,meta_interactive_1_ml_inc_redlist) #Not significant = no interaction present

predict(meta_interactive_1_inc_redlist,newmods=c(0,0,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #artificial x reference level
predict(meta_interactive_1_inc_redlist,newmods=c(1,0,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #natural x reference level

## Red List Classification : Incubation Type
meta_interactive_2_inc_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Incubation):relevel(factor(IUCN.Status.Current_SUB),ref="CR")-1, test="t", method="REML", data=ies.da)
summary(meta_interactive_2_inc_redlist)

## Red List Classification + Supplementary Feeding - Table 2 and Table S19
meta_additive_supp_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Supplementary.Feeding)+relevel(factor(IUCN.Status.Current_SUB),ref="LC"), test="t", method="REML", data=ies.da)
summary(meta_additive_supp_redlist) #Supplementary Feeding is significant
anova(meta_additive_supp_redlist,btt=3:6) #Red List Classification is significant
predict(meta_additive_supp_redlist,newmods=rbind(c(0,0,0,0,0),c(1,0,0,0,0)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) 

meta_additive_nointercept_supp_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Supplementary.Feeding)+relevel(factor(IUCN.Status.Current_SUB),ref="LC")-1, test="t", method="REML", data=ies.da)
summary(meta_additive_nointercept_supp_redlist)

## Red List Classification * Supplementary.Feeding - Table S19
meta_interactive_1_supp_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Supplementary.Feeding)*relevel(factor(IUCN.Status.Current_SUB),ref="LC"), test="t", method="REML", data=ies.da)
summary(meta_interactive_1_supp_redlist)
anova(meta_interactive_1_supp_redlist,btt=7:10) #Not significant = no interaction present

meta_additive_ml_supp_redlist <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Supplementary.Feeding)+factor(IUCN.Status.Current_SUB),test="t",method="ML",data=ies.da)
meta_interactive_1_ml_supp_redlist  <-rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Supplementary.Feeding)*factor(IUCN.Status.Current_SUB),test="t",method="ML",data=ies.da)

anova(meta_additive_ml_supp_redlist,meta_interactive_1_ml_supp_redlist) #Not significant = no interaction present

predict(meta_interactive_1_supp_redlist,newmods=c(0,0,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #not fed x reference level
predict(meta_interactive_1_supp_redlist,newmods=c(1,0,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #fed x reference level

## Red List Classification : Supplementary.Feeding
meta_interactive_2_supp_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Supplementary.Feeding):relevel(factor(IUCN.Status.Current_SUB),ref="CR")-1, test="t", method="REML", data=ies.da)
summary(meta_interactive_2_supp_redlist)

## Red List Classification + Artificial Nest Provision - Table 2 and Table S20
meta_additive_nest_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Nestboxes.Burrows.Sites)+relevel(factor(IUCN.Status.Current_SUB),ref="LC"), test="t", method="REML", data=ies.da)
summary(meta_additive_nest_redlist) #Artificial Nest Provision is significant
anova(meta_additive_nest_redlist,btt=3:6) #Red List Classification is significant
predict(meta_additive_nest_redlist,newmods=rbind(c(0,0,0,0,0),c(1,0,0,0,0)),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) 

meta_additive_nointercept_nest_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Nestboxes.Burrows.Sites)+relevel(factor(IUCN.Status.Current_SUB),ref="LC")-1, test="t", method="REML", data=ies.da)
summary(meta_additive_nointercept_nest_redlist)

## Red List Classification * Artificial Nest Provision - Table S20
meta_interactive_1_nest_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Nestboxes.Burrows.Sites)*relevel(factor(IUCN.Status.Current_SUB),ref="CR"), test="t", method="REML", data=ies.da)
summary(meta_interactive_1_nest_redlist)
anova(meta_interactive_1_nest_redlist,btt=7:9) #Not significant = no interaction present

meta_additive_ml_nest_redlist <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Nestboxes.Burrows.Sites)+factor(IUCN.Status.Current_SUB),test="t",method="ML",data=ies.da)
meta_interactive_1_ml_nest_redlist  <-rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor),mods = ~factor(Nestboxes.Burrows.Sites)*factor(IUCN.Status.Current_SUB),test="t",method="ML",data=ies.da)

anova(meta_additive_ml_nest_redlist,meta_interactive_1_ml_nest_redlist) #Not significant = no interaction present

predict(meta_interactive_1_nest_redlist,newmods=c(0,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #not provided x reference level
predict(meta_interactive_1_nest_redlist,newmods=c(1,0,0,0,0,0,0,0),transf=transf.ipft.hm, targ=list(ni=ies.da$Total)) #provided x reference level

## Red List Classification : Artificial Nest Provision
meta_interactive_2_nest_redlist <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor), mods = ~factor(Nestboxes.Burrows.Sites):relevel(factor(IUCN.Status.Current_SUB),ref="LC")-1, test="t", method="REML", data=ies.da)
summary(meta_interactive_2_nest_redlist)

#### ROBUSTNESS CHECKS ####

## One intervention applied at a time - Table S24
## Only artificial incubation
data_final_onlyinc <- data_final[data_final$Incubation=="Artificial" & data_final$Supplementary.Feeding=="No" & data_final$Nestboxes.Burrows.Sites=="No",]
## Only supplementary feeding
data_final_onlyfed <- data_final[data_final$Incubation=="Natural" & data_final$Supplementary.Feeding=="Yes" & data_final$Nestboxes.Burrows.Sites=="No",]
## Only Artificial Nest Provision
data_final_onlynest <- data_final[data_final$Incubation=="Natural" & data_final$Supplementary.Feeding=="No" & data_final$Nestboxes.Burrows.Sites=="Yes",]
## Only Wild populations
data_final_wild <- data_final[data_final$Wild.Managed.Captive=="Wild" & !is.na(data_final$Wild.Managed.Captive),]

data_final_onlyinc.v.wild <- rbind(data_final_onlyinc,data_final_wild)
data_final_onlyfed.v.wild <- rbind(data_final_onlyfed,data_final_wild)
data_final_onlynest.v.wild <- rbind(data_final_onlynest,data_final_wild)

## Trees and data for only artificially incubated vs. wild
tips_onlyinc.v.wild <- data_final_onlyinc.v.wild$BirdTree.Tip
pruned.trees_onlyinc.v.wild <- lapply(trees,keep.tip,tip=tips_onlyinc.v.wild)
class(pruned.trees_onlyinc.v.wild) <- "multiPhylo"
mcc_pruned_onlyinc.v.wild <- maxCladeCred(pruned.trees_onlyinc.v.wild)
mcc_pruned_onlyinc.v.wild <- ladderize(mcc_pruned_onlyinc.v.wild)
mcc_multi_onlyinc.v.wild <- mcc_pruned_onlyinc.v.wild$tip.label
mcc_ult_onlyinc.v.wild <- compute.brlen(mcc_pruned_onlyinc.v.wild, power = 1)
bird_phylo_cor_onlyinc.v.wild <- vcv(mcc_ult_onlyinc.v.wild, cor=T) 

ies.da_onlyinc.v.wild=escalc(xi=Total-Cases, ni=Total, data=data_final_onlyinc.v.wild, measure="PFT", add=0)

## Artificial Incubation - only artificially incubated compared to wild
overall_phylo_sp_onlyinc.v.wild <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_onlyinc.v.wild),
                                          mods = ~relevel(factor(Incubation),ref="Natural"), test="t", method="REML", data=ies.da_onlyinc.v.wild)
summary(overall_phylo_sp_onlyinc.v.wild) #Incubation Type is significant
predict(overall_phylo_sp_onlyinc.v.wild,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da_onlyinc.v.wild$Total))

## Tree and data for only supplementary feeding only vs. wild
tips_onlyfed.v.wild <- data_final_onlyfed.v.wild$BirdTree.Tip
pruned.trees_onlyfed.v.wild <- lapply(trees,keep.tip,tip=tips_onlyfed.v.wild)
class(pruned.trees_onlyfed.v.wild) <- "multiPhylo"
mcc_pruned_onlyfed.v.wild <- maxCladeCred(pruned.trees_onlyfed.v.wild)
mcc_pruned_onlyfed.v.wild <- ladderize(mcc_pruned_onlyfed.v.wild)
mcc_multi_onlyfed.v.wild <- mcc_pruned_onlyfed.v.wild$tip.label
mcc_ult_onlyfed.v.wild <- compute.brlen(mcc_pruned_onlyfed.v.wild, power = 1)
bird_phylo_cor_onlyfed.v.wild <- vcv(mcc_ult_onlyfed.v.wild, cor=T) 

ies.da_onlyfed.v.wild=escalc(xi=Total-Cases, ni=Total, data=data_final_onlyfed.v.wild, measure="PFT", add=0)

## Supplementary Feeding - only supplementary fed compared to wild
overall_phylo_sp_onlyfed.v.wild <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_onlyfed.v.wild), mods = ~relevel(factor(Supplementary.Feeding),ref="Yes"), test="t", method="REML", data=ies.da_onlyfed.v.wild)
summary(overall_phylo_sp_onlyfed.v.wild) #Supplementary Feeding is significant
predict(overall_phylo_sp_onlyfed.v.wild,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da_onlyfed.v.wild$Total))

## Tree and data for only artificial nest provision only vs. wild
tips_onlynest.v.wild <- data_final_onlynest.v.wild$BirdTree.Tip
pruned.trees_onlynest.v.wild <- lapply(trees,keep.tip,tip=tips_onlynest.v.wild)
class(pruned.trees_onlynest.v.wild) <- "multiPhylo"
mcc_pruned_onlynest.v.wild <- maxCladeCred(pruned.trees_onlynest.v.wild)
mcc_pruned_onlynest.v.wild <- ladderize(mcc_pruned_onlynest.v.wild)
mcc_multi_onlynest.v.wild <- mcc_pruned_onlynest.v.wild$tip.label
mcc_ult_onlynest.v.wild <- compute.brlen(mcc_pruned_onlynest.v.wild, power = 1)
bird_phylo_cor_onlynest.v.wild <- vcv(mcc_ult_onlynest.v.wild, cor=T) 

ies.da_onlynest.v.wild=escalc(xi=Total-Cases, ni=Total, data=data_final_onlynest.v.wild, measure="PFT", add=0)

#Artificial Nest Provision - only provided nests compared to wild
overall_phylo_sp_onlynest.v.wild <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_onlynest.v.wild), mods = ~relevel(factor(Nestboxes.Burrows.Sites),ref="No"), test="t", method="REML", data=ies.da_onlynest.v.wild)
summary(overall_phylo_sp_onlynest.v.wild) #Artificial Nest Provision is NOT significant
predict(overall_phylo_sp_onlynest.v.wild,newmods=c(0,1),transf=transf.ipft.hm, targ=list(ni=ies.da_onlynest.v.wild$Total))

## Only populations of threatened species
data_final_threatened <- data_final[data_final$Threat_Status_SUB=="Threatened",]
ies.da_threatened=escalc(xi=Total-Cases, ni=Total, data=data_final_threatened, measure="PFT", add=0)

tips_threatened <- data_final_threatened$BirdTree.Tip
pruned.trees_threatened <- lapply(trees,keep.tip,tip=tips_threatened)
class(pruned.trees_threatened) <- "multiPhylo"
mcc_pruned_threatened <- maxCladeCred(pruned.trees_threatened)
mcc_pruned_threatened <- ladderize(mcc_pruned_threatened)
mcc_multi_threatened <- mcc_pruned_threatened$tip.label
mcc_ult_threatened <- compute.brlen(mcc_pruned_threatened, power = 1)
bird_phylo_cor_threatened <- vcv(mcc_ult_threatened, cor=T) 

overall_phylo_sp_manage_THREAT <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_threatened), mods = ~relevel(factor(Wild.Managed.Captive),ref="Wild"), test="t", method="REML", data=ies.da_threatened)
summary(overall_phylo_sp_manage_THREAT) #Management Level is significant
predict(overall_phylo_sp_manage_THREAT,newmods=rbind(c(0,0), c(1,0), c(0,1)),transf=transf.ipft.hm, targ=list(ni=ies.da_threatened$Total))

## Only populations of non-threatened species
data_final_nonthreatened <- data_final[data_final$Threat_Status_SUB=="Non-threatened",]
ies.da_nonthreatened=escalc(xi=Total-Cases, ni=Total, data=data_final_nonthreatened, measure="PFT", add=0)

tips_nonthreatened <- data_final_nonthreatened$BirdTree.Tip
pruned.trees_nonthreatened <- lapply(trees,keep.tip,tip=tips_nonthreatened)
class(pruned.trees_nonthreatened) <- "multiPhylo" #Need to change back to multiPhylo because lapply makes it a list
mcc_pruned_nonthreatened <- maxCladeCred(pruned.trees_nonthreatened)
mcc_pruned_nonthreatened <- ladderize(mcc_pruned_nonthreatened)
mcc_multi_nonthreatened <- mcc_pruned_nonthreatened$tip.label
mcc_ult_nonthreatened <- compute.brlen(mcc_pruned_nonthreatened, power = 1) 
bird_phylo_cor_nonthreatened <- vcv(mcc_ult_nonthreatened, cor=T) 

overall_phylo_sp_manage_NONTHREAT <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_nonthreatened), mods = ~relevel(factor(Wild.Managed.Captive),ref="Wild"), test="t", method="REML", data=ies.da_nonthreatened)
summary(overall_phylo_sp_manage_NONTHREAT) #Management Level is significant
predict(overall_phylo_sp_manage_NONTHREAT,newmods=rbind(c(0,0), c(1,0), c(0,1)),transf=transf.ipft.hm, targ=list(ni=ies.da_nonthreatened$Total))

#### SENSITIVITY ANALYSES ####

## Excluding Briskie & Mackintosh (2004) due to uncertainty around management details (see Appendix S4 for further details)
data_final_noBriskie <- data_final[data_final$Study.Name!="Briskie & Mackintosh 2004",]
ies.da_no.Briskie=escalc(xi=Total-Cases, ni=Total, data=data_final_noBriskie, measure="PFT", add=0) 

final_species_no.Briskie <-sort(unique(as.character(data_final_noBriskie$BirdTree.Tip))) #219

tips_no.Briskie <- data_final_noBriskie$BirdTree.Tip
pruned.trees_no.Briskie <- lapply(trees,keep.tip,tip=tips_no.Briskie)
class(pruned.trees_no.Briskie) <- "multiPhylo"
mcc_pruned_no.Briskie <- maxCladeCred(pruned.trees_no.Briskie)
mcc_pruned_no.Briskie <- ladderize(mcc_pruned_no.Briskie)
mcc_pruned_no.Briskie$tip.label #219 species
mcc_multi_no.Briskie <- mcc_pruned_no.Briskie$tip.label
mcc_ult_no.Briskie <- compute.brlen(mcc_pruned_no.Briskie, power = 1)
bird_phylo_cor_no.Briskie <- vcv(mcc_ult_no.Briskie, cor=T) 

overall_phylo_sp_no.Briskie <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor_no.Briskie),test="t",method="REML",data=ies.da_no.Briskie)
summary(overall_phylo_sp_no.Briskie)
predict(overall_phylo_sp_no.Briskie, transf=transf.ipft.hm, targ=list(ni=ies.da_no.Briskie$Total))

W <- diag(1/ies.da_no.Briskie$vi)
X <- model.matrix(overall_phylo_sp_no.Briskie)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

100 * sum(overall_phylo_sp_no.Briskie$sigma2) / (sum(overall_phylo_sp_no.Briskie$sigma2) + (overall_phylo_sp_no.Briskie$k-overall_phylo_sp_no.Briskie$p)/sum(diag(P))) #Total I^2 = 99.17601
100 * overall_phylo_sp_no.Briskie$sigma2 / (sum(overall_phylo_sp_no.Briskie$sigma2) + (overall_phylo_sp_no.Briskie$k-overall_phylo_sp_no.Briskie$p)/sum(diag(P))) #Within-study I^2 = 25.86437; Between-study I^2 = 19.87115; Between-species I^2 = 18.48729; Phylogenetic history = 34.95320 - Table S9

## Excluding potential outliers identified through externally studentized residuals and Cook's Distance
data_final_no.outliers<-data_final[-c(13,25,76,187,197,201,207,214,258,259,264,291,306,314,366,401,402,404,471,472),]
ies.da_no.outliers=escalc(xi=Total-Cases, ni=Total, data=data_final_no.outliers, measure="PFT", add=0)
dim(data_final_no.outliers)
final_species_no.outliers <- sort(unique(as.character(data_final_no.outliers$BirdTree.Tip))) #228

tips_no.outliers <- data_final_no.outliers$BirdTree.Tip
pruned.trees_no.outliers <- lapply(trees,keep.tip,tip=tips_no.outliers)
class(pruned.trees_no.outliers)<-"multiPhylo"
mcc_pruned_no.outliers <- maxCladeCred(pruned.trees_no.outliers)
mcc_pruned_no.outliers <- ladderize(mcc_pruned_no.outliers)
mcc_pruned_no.outliers$tip.label #228 species
mcc_multi_no.outliers <- mcc_pruned_no.outliers$tip.label
mcc_ult_no.outliers <- compute.brlen(mcc_pruned_no.outliers, power = 1) 
bird_phylo_cor_no.outliers <- vcv(mcc_ult_no.outliers, cor=T)

overall_phylo_sp_no.outliers <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor_no.outliers),test="t",method="REML",data=ies.da_no.outliers)
summary(overall_phylo_sp_no.outliers)
predict(overall_phylo_sp_no.outliers, transf=transf.ipft.hm, targ=list(ni=ies.da_no.outliers$Total))

W <- diag(1/ies.da_no.outliers$vi)
X <- model.matrix(overall_phylo_sp_no.outliers)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

100 * sum(overall_phylo_sp_no.outliers$sigma2) / (sum(overall_phylo_sp_no.outliers$sigma2) + (overall_phylo_sp_no.outliers$k-overall_phylo_sp_no.outliers$p)/sum(diag(P))) #Total I^2 = 98.85184
100 * overall_phylo_sp_no.outliers$sigma2 / (sum(overall_phylo_sp_no.outliers$sigma2) + (overall_phylo_sp_no.outliers$k-overall_phylo_sp_no.outliers$p)/sum(diag(P))) #Within-study I^2 = 33.06409; Between-study I^2 = 21.39931; Between-species I^2 = 25.57844; Phylogenetic history = 18.80999 - Table S9

## Excluding records with sample sizes estimated based on mean clutch size from Cooney et al. (2020)
data_final_no.Cooney <- data_final[data_final$TotalEggs_Summary!="Estimated from Cooney",]
ies.da_no.Cooney=escalc(xi=Total-Cases, ni=Total, data=data_final_no.Cooney, measure="PFT", add=0) 

final_species_no.Cooney <-sort(unique(as.character(data_final_no.Cooney$BirdTree.Tip))) #210

tips_no.Cooney <- data_final_no.Cooney$BirdTree.Tip
pruned.trees_no.Cooney <- lapply(trees,keep.tip,tip=tips_no.Cooney)
class(pruned.trees_no.Cooney) <- "multiPhylo"
mcc_pruned_no.Cooney <- maxCladeCred(pruned.trees_no.Cooney)
mcc_pruned_no.Cooney <- ladderize(mcc_pruned_no.Cooney)
mcc_pruned_no.Cooney$tip.label #210 species
mcc_multi_no.Cooney <- mcc_pruned_no.Cooney$tip.label
mcc_ult_no.Cooney <- compute.brlen(mcc_pruned_no.Cooney, power = 1)
bird_phylo_cor_no.Cooney <- vcv(mcc_ult_no.Cooney, cor=T) 

overall_phylo_sp_no.Cooney <- rma.mv(yi, vi,random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo),R = list(phylo = bird_phylo_cor_no.Cooney),test="t",method="REML",data=ies.da_no.Cooney)
summary(overall_phylo_sp_no.Cooney)
predict(overall_phylo_sp_no.Cooney, transf=transf.ipft.hm, targ=list(ni=ies.da_no.Cooney$Total))

W <- diag(1/ies.da_no.Cooney$vi)
X <- model.matrix(overall_phylo_sp_no.Cooney)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

100 * sum(overall_phylo_sp_no.Cooney$sigma2) / (sum(overall_phylo_sp_no.Cooney$sigma2) + (overall_phylo_sp_no.Cooney$k-overall_phylo_sp_no.Cooney$p)/sum(diag(P))) #Total I^2 = 99.28562
100 * overall_phylo_sp_no.Cooney$sigma2 / (sum(overall_phylo_sp_no.Cooney$sigma2) + (overall_phylo_sp_no.Cooney$k-overall_phylo_sp_no.Cooney$p)/sum(diag(P))) #Within-study I^2 = 23.26942; Between-study I^2 = 15.38820; Between-species I^2 = 18.99683; Phylogenetic history = 41.63116 - Table S9

## NOTE: Re-run all mixed effects multilevel meta-analysis models for the _no.outliers, _no.Briskie, and _no.Cooney datasets for Tables S10, S25 and S26

#### PUBLICATION BIAS ####

## Taxonomy - Table S21
## Comparing final dataset to global distribution - By number of records
data_final %>%
  count(Order)
dataset_taxa_records = c(3,24,2,1,0,6,58,3,0,3,3,2,1,14,25,0,13,0,0,0,0,2,245,12,0,0,4,1,10,12,6,11,11,5,1,0)
global_taxa <- c(251/10994,170/10994,72/10994,597/10994,2/10994,7/10994,377/10994,20/10994,6/10994,352/10994,187/10994,150/10994,2/10994,64/10994,306/10994,5/10994,169/10994,1/10994,3/10994,24/10994,1/10994,26/10994,6602/10994,110/10994,3/10994,6/10994,484/10994,20/10994,145/10994,403/10994,16/10994,18/10994,238/10994,61/10994,53/10994,43/10994)
taxa_records_result <- chisq.test(x=dataset_taxa_records,p=global_taxa)
taxa_records_result #Significant p value = reject H0 = final data doesn't follow global distribution
taxa_records_result$expected
taxa_bias_records <- matrix(data=c(dataset_taxa_records,taxa_records_result$expected),nrow=2,ncol=36,byrow=T)
dimnames(taxa_bias_records) = list(Dataset=c('Final_data','World'),Order=c('ACCIPITRIFORMES','ANSERIFORMES','BUCEROTIFORMES','CAPRIMULGIFORMES','CARIAMIFORMES','CATHARTIFORMES','CHARADRIIFORMES','CICONIIFORMES','COLIIFORMES','COLUMBIFORMES','CORACIIFORMES','CUCULIFORMES','EURYPYGIFORMES','FALCONIFORMES','GALLIFORMES','GAVIIFORMES','GRUIFORMES','LEPTOSOMIFORMES','MESITORNITHIFORMES','MUSOPHAGIFORMES','OPISTHOCOMIFORMES','OTIDIFORMES','PASSERIFORMES','PELECANIFORMES','PHAETHONTIFORMES','PHOENICOPTERIFORMES','PICIFORMES','PODICIPEDIFORMES','PROCELLARIIFORMES','PSITTACIFORMES','PTEROCLIFORMES','SPHENISCIFORMES','STRIGIFORMES','STRUTHIONIFORMES','SULIFORMES','TROGONIFORMES'))
taxa_bias_records
## Comparing final dataset to global distribution - By number of species
data_final[data_final$Order=="Accipitriformes",] %>% 
  count(BirdTree.Tip) 
dataset_taxa_species = c(2,11,1,1,0,1,36,1,0,3,3,2,1,5,14,0,7,0,0,0,0,1,106,8,0,0,4,1,6,10,1,6,5,4,1,0)
global_taxa <- c(251/10994,170/10994,72/10994,597/10994,2/10994,7/10994,377/10994,20/10994,6/10994,352/10994,187/10994,150/10994,2/10994,64/10994,306/10994,5/10994,169/10994,1/10994,3/10994,24/10994,1/10994,26/10994,6602/10994,110/10994,3/10994,6/10994,484/10994,20/10994,145/10994,403/10994,16/10994,18/10994,238/10994,61/10994,53/10994,43/10994)
taxa_species_result <- chisq.test(x=dataset_taxa_species,p=global_taxa)
taxa_species_result #Significant p value = reject H0 = final data doesn't follow global distribution
taxa_species_result$expected
taxa_bias_species <- matrix(data=c(dataset_taxa_species,taxa_species_result$expected),nrow=2,ncol=36,byrow=T)
dimnames(taxa_bias_species) = list(Dataset=c('Final_data','World'),Order=c('ACCIPITRIFORMES','ANSERIFORMES','BUCEROTIFORMES','CAPRIMULGIFORMES','CARIAMIFORMES','CATHARTIFORMES','CHARADRIIFORMES','CICONIIFORMES','COLIIFORMES','COLUMBIFORMES','CORACIIFORMES','CUCULIFORMES','EURYPYGIFORMES','FALCONIFORMES','GALLIFORMES','GAVIIFORMES','GRUIFORMES','LEPTOSOMIFORMES','MESITORNITHIFORMES','MUSOPHAGIFORMES','OPISTHOCOMIFORMES','OTIDIFORMES','PASSERIFORMES','PELECANIFORMES','PHAETHONTIFORMES','PHOENICOPTERIFORMES','PICIFORMES','PODICIPEDIFORMES','PROCELLARIIFORMES','PSITTACIFORMES','PTEROCLIFORMES','SPHENISCIFORMES','STRIGIFORMES','STRUTHIONIFORMES','SULIFORMES','TROGONIFORMES'))
taxa_bias_species

## Threat Status - Table S22
## Comparing final dataset to global distribution - By number of records
data_final %>%
  count(IUCN.Status.Current_SUB)
dataset_threat_records <- c(347,35,37,39,20)
global_threat <- c(8460/10942,1001/10942,798/10942,460/10942,223/10942)
threat_records_result <- chisq.test(x=dataset_threat_records,p=global_threat)
threat_records_result #Significant p value = reject H0 = final data doesn't follow global distribution
threat_records_result$expected #Expected number in each category: LC = 370, NT = 44, VU = 35, EN = 20, CR = 10
threat_bias_records <- matrix(data=c(dataset_threat_records,threat_records_result$expected),nrow=2,ncol=5,byrow=T)
dimnames(threat_bias_records) = list(Dataset=c('Final_data','World'),RedList=c('LC','NT','VU','EN','CR'))
threat_bias_records
## Comparing final dataset to global distribution - By number of species
data_final[data_final$IUCN.Status.Current_SUB=="CR",] %>% 
  count(BirdTree.Tip) 
dataset_threat_species <- c(158,23,25,25,13)
global_threat <- c(8460/10942,1001/10942,798/10942,460/10942,223/10942)
threat_species_result <- chisq.test(x=dataset_threat_species,p=global_threat)
threat_species_result #Significant p value = reject H0 = final data doesn't follow global distribution
threat_species_result$expected #Expected number in each category: LC = 193, NT = 23, VU = 18, EN = 10, CR = 5
threat_bias_species <- matrix(data=c(dataset_threat_species,threat_species_result$expected),nrow=2,ncol=5,byrow=T)
dimnames(threat_bias_species) = list(Dataset=c('Final_data','World'),RedList=c('LC','NT','VU','EN','CR'))
threat_bias_species

## Geographical Distribution - Figure 2
world <- map_data("world")
includedMAP <- data_final  %>% 
  filter(Longitude != "NA")
long <- factor(includedMAP$Longitude)
lat <- factor(includedMAP$Latitude)
as.numeric(as.character(long))
as.numeric(as.character(lat))
glimpse(includedMAP)
ggplot() +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey", alpha=0.3) +
  geom_point(data=includedMAP, aes(x=as.numeric(as.character(long)), y=as.numeric(as.character(lat)), color=IUCN.Status.Current_SUB)) +
  scale_colour_manual(name="IUCN Red List Status", breaks = c("CR","EN","VU","NT","LC"),
                      values=c("#C00000","#ED7D31","#FFC000","#92D050","#00B050")) +
  theme_void() +
  ggtitle("Study Locations of Final Dataset Population Records")

#### ADDITIONAL FIGURES ####

## Figure 3 - Forest Plot
forest_data <- tibble(mean=c(NA,21.02,15.24,NA,27.95,23.41,16.60,17.57,15.13,NA,38.07,20.05,13.78,NA,37.38,14.38,NA,33.47,15.37,NA,21.34,15.89,NA,16.79),
                      lower = c(NA,11.62,7.42,NA,15.26,12.83,7.61,8.22,7.16,NA,26.60,12.02,7.25,NA,26.04,7.37,NA,23.03,8.57,NA,11.34,7.57,NA,8.28),
                      upper = c(NA,32.22,25.07,NA,42.67,35.93,28.02,29.35,25.24,NA,50.23,29.46,21.87,NA,49.45,23.09,NA,44.77,23.62,NA,33.37,26.38,NA,27.40),
                      category = c(NA,"Threatened","Non-threatened",NA,"Critically Endangered","Endangered","Vulnerable","Near Threatened","Least Concern",NA,"Captive","Wild Managed","Wild",NA,"Artificial","Natural",NA,"Fed","Not Fed",NA,"Provided","Not Provided",NA,"All"),
                      k = c(NA,"96","382",NA,"20","39","37","35","347",NA,"41","172","265",NA,"53","425",NA,"51","427",NA,"161","317",NA,"478"))
header <- tibble(category=c("","Category"), k = c("","k"), summary=TRUE)
empty_row <- tibble(mean=NA_real_)
final_forest <- bind_rows(header, forest_data, empty_row)

final_forest %>% 
  forestplot(labeltext = c(category,k), title = "Forest Plot of Mean Effect Size (Hatching Failure %)", xlab = "Mean Hatching Failure (%)",
             clip = c(0, 55), zero = 0, boxsize = 0.5,vertices = TRUE,
             xticks = c(0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0,47.5,50.0,52.5,55.0),
             txt_gp = fpTxtGp(
               label = list(
                 gpar(fontface="plain"),gpar(fontface="plain")),
               ticks = gpar(cex=0.6),
               xlab = gpar(cex=1)),
             fn.ci_norm = c(fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawDiamondCI,fpDrawNormalCI,fpDrawDiamondCI),
             grid = structure(c(5,10,15,20,25,30,35,40,45,50,55),
                              gp=gpar(lty=1,col="lightgray")),
             shapes_gp = fpShapesGp(
               box = list(
                 gpar(fill=NA),gpar(fill=NA),gpar(fill=NA),gpar(fill="#f1eef6"),gpar(fill="#f1eef6"),gpar(fill=NA),gpar(fill="#d0d1e6"),gpar(fill="#d0d1e6"),gpar(fill="#d0d1e6"),gpar(fill="#d0d1e6"),gpar(fill="#d0d1e6"),gpar(fill=NA),gpar(fill="#a6bddb"),gpar(fill="#a6bddb"),gpar(fill="#a6bddb"),gpar(fill=NA),gpar(fill="#74a9cf"),gpar(fill="#74a9cf"),gpar(fill=NA),gpar(fill="#2b8cbe"),gpar(fill="#2b8cbe"),gpar(fill=NA),gpar(fill="#045a8d"),gpar(fill="#045a8d"),gpar(fill=NA),gpar(fill="#c51b8a"),gpar(fill=NA)),
               lines = gpar(col="black"),
               zero = gpar(col="black")))

## Figure S3 - Bubble Plot (publication year)
ies.da$size <- (1/sqrt(ies.da$vi))/500
ies.da$HF<-((data_final$Total-data_final$Cases)/(data_final$Total))*100
plot(NA,NA,xlim=c(2004,2020),ylim=c(0,80), xlab="Publication Year", ylab="Hatching Failure (%)", las=1,bty="l", xaxt="none")
axis(1,seq(2004,2020,2))
symbols(ies.da$Publication.Year, ies.da$HF, circles=ies.da$size, inches=FALSE, add=TRUE, fg= rgb(0,0,0, max=255, alpha=150), pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150))
lines(2004:2020, (preds$pred)*100, lty="dashed")
lines(2004:2020, (preds$ci.lb)*100, lty="dotted")
lines(2004:2020, (preds$ci.ub)*100, lty="dotted")

#slope estimate of the fitted line
x <- 2004:2020
y <- (preds$pred)*100
plot(x,y)
abline(mod <- lm(y~x))
coef(mod) #slope estimate is -0.00199844 

## Figure S4 - Bubble Plot - (original source publication year)
data_final_originalyear <- data_final[-which(is.na(data_final$Original.Source_Publication.Year)),]

data_final_originalyear %>%
  count(Original.Source_Publication.Year)
ies.da_originalyear=escalc(xi=Total-Cases, ni=Total, data=data_final_originalyear, measure="PFT", add=0)

tips_originalyear <- data_final_originalyear$BirdTree.Tip
pruned.trees_originalyear <- lapply(trees,keep.tip,tip=tips_originalyear)
class(pruned.trees_originalyear)<-"multiPhylo"
mcc_pruned_originalyear <- maxCladeCred(pruned.trees_originalyear)
mcc_pruned_originalyear <- ladderize(mcc_pruned_originalyear)
mcc_pruned_originalyear$tip.label #219 species
mcc_multi_originalyear <- mcc_pruned_originalyear$tip.label
mcc_ult_originalyear <- compute.brlen(mcc_pruned_originalyear, power = 1)
bird_phylo_cor_originalyear <- vcv(mcc_ult_originalyear, cor=T) 

overall_phylo_sp_pubyear_original <- rma.mv(yi, vi, random=list(~1|effect_size_id,~1|Study.Name,~1|BirdTree.Tip,~1|phylo), R = list(phylo = bird_phylo_cor_originalyear), mods = ~Original.Source_Publication.Year, test="t", method="REML", data=ies.da_originalyear)
summary(overall_phylo_sp_pubyear_original)
preds_original=predict(overall_phylo_sp_pubyear_original,newmods=c(1940,1949,1959,1961,1963,1969,1971,1973,1981,1982,1983,1984,1986,1987,1988,1990,1992,1993,1994,1995,1997,1999:2020),transf=transf.ipft.hm, targ=list(ni=ies.da_originalyear$Total))

ies.da_originalyear$size <- (1/sqrt(ies.da_originalyear$vi))/500
ies.da_originalyear$HF<-((data_final_originalyear$Total-data_final_originalyear$Cases)/(data_final_originalyear$Total))*100
plot(NA,NA,xlim=c(1940,2020),ylim=c(0,80), xlab="Publication Year", ylab="Hatching Failure (%)", las=1,bty="l", xaxt="none")
axis(1,seq(1940,2020,2))
symbols(ies.da_originalyear$Original.Source_Publication.Year, ies.da_originalyear$HF, circles=ies.da_originalyear$size, inches=FALSE, add=TRUE, fg= rgb(0,0,0, max=255, alpha=150), pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150))
lines(c(1940,1949,1959,1961,1963,1969,1971,1973,1981,1982,1983,1984,1986,1987,1988,1990,1992,1993,1994,1995,1997,1999:2020), (preds_original$pred)*100, lty="dashed")
lines(c(1940,1949,1959,1961,1963,1969,1971,1973,1981,1982,1983,1984,1986,1987,1988,1990,1992,1993,1994,1995,1997,1999:2020), (preds_original$ci.lb)*100, lty="dotted")
lines(c(1940,1949,1959,1961,1963,1969,1971,1973,1981,1982,1983,1984,1986,1987,1988,1990,1992,1993,1994,1995,1997,1999:2020), (preds_original$ci.ub)*100, lty="dotted")

#slope estimate of the fitted line
x <- c(1940,1949,1959,1961,1963,1969,1971,1973,1981,1982,1983,1984,1986,1987,1988,1990,1992,1993,1994,1995,1997,1999:2020)
y <- (preds_original$pred)*100
plot(x,y)
abline(mod <- lm(y~x))
coef(mod) #slope estimate is 0.176156

## Figure S5
## Final Dataset
my.labels_tax_final <- c("Accipitriformes (n=3)","Anseriformes (n=24)","Bucerotiformes (n=2)","Caprimulgiformes (n=1)","Cathartiformes (n=6)","Charadriiformes (n=58)","Ciconiiformes (n=3)","Columbiformes (n=3)","Coraciiformes (n=3)","Cuculiformes (n=2)","Eurypygiformes (n=1)","Falconiformes (n=14)","Galliformes (n=25)","Gruiformes (n=13)","Otidiformes (n=2)","Passeriformes (n=245)","Pelecaniformes (n=12)","Piciformes (n=4)","Podicipediformes (n=1)","Procellariiformes (n=10)","Psittaciformes (n=12)","Pterocliformes (n=6)","Sphenisciformes (n=11)","Strigiformes (n=11)","Struthioniformes (n=5)","Suliformes (n=1)")
taxonomy_plot_final <- data_final %>%
  mutate(Order = fct_relevel(Order,"Accipitriformes","Anseriformes","Bucerotiformes","Caprimulgiformes","Cathartiformes","Charadriiformes","Ciconiiformes","Columbiformes","Coraciiformes","Cuculiformes","Eurypygiformes","Falconiformes","Galliformes","Gruiformes","Otidiformes","Passeriformes","Pelecaniformes","Piciformes","Podicipediformes","Procellariiformes","Psittaciformes","Pterocliformes","Sphenisciformes","Strigiformes","Struthioniformes","Suliformes")) %>%
  ggplot(aes(x=Order, y=Mean.Hatching.Failure.std, fill=Order)) + 
  geom_boxplot() + 
  labs(title="Plot of Mean Hatching Failure % according to Taxonomic Order\nFinal Dataset (n = 478)", x="Order", y="Mean Hatching Failure (%)")
taxonomy_final <- taxonomy_plot_final + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=90,hjust=0.95)) + 
  scale_x_discrete(label=my.labels_tax_final) +  
  scale_fill_discrete(name="Order",labels=my.labels_tax_final)
taxonomy_final

#Final Dataset - limited to orders with at least 10 effect sizes
data_final_10 <- data_final[data_final$Order!="Accipitriformes" & data_final$Order!="Bucerotiformes" & data_final$Order!="Caprimulgiformes" & data_final$Order!="Cathartiformes" & data_final$Order!="Ciconiiformes" & data_final$Order!="Columbiformes" & data_final$Order!="Coraciiformes" & data_final$Order!="Cuculiformes" & data_final$Order!="Eurypygiformes" & data_final$Order!="Otidiformes" & data_final$Order!="Piciformes" & data_final$Order!="Podicipediformes" & data_final$Order!="Pterocliformes" & data_final$Order!="Struthioniformes" & data_final$Order!="Suliformes", ]
my.labels_tax_final_10 <- c("Anseriformes (n=24)","Charadriiformes (n=58)","Falconiformes (n=14)","Galliformes (n=25)","Gruiformes (n=13)","Passeriformes (n=245)","Pelecaniformes (n=12)","Procellariiformes (n=10)","Psittaciformes (n=12)","Sphenisciformes (n=11)","Strigiformes (n=11)")
taxonomy_plot_final_10 <- data_final_10 %>%
  mutate(Order = fct_relevel(Order,"Anseriformes","Charadriiformes","Falconiformes","Galliformes","Gruiformes","Passeriformes","Pelecaniformes","Procellariiformes","Psittaciformes","Sphenisciformes","Strigiformes")) %>%
  ggplot(aes(x=Order, y=Mean.Hatching.Failure.std, fill=Order)) + 
  geom_boxplot() + 
  labs(title="Plot of Mean Hatching Failure % according to Taxonomic Order\nFinal Dataset - at least 10 records per order", x="Order", y="Mean Hatching Failure (%)")
taxonomy_final_10 <- taxonomy_plot_final_10 + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=90,hjust=0.95)) +  
  scale_x_discrete(label=my.labels_tax_final_10) +  
  scale_fill_discrete(name="Order",labels=my.labels_tax_final_10) + 
  scale_fill_manual(values=c("#EE8043","#A8A401","#3FC090","#3FC2AB","#3DBFC4","#37B5ED","#32ACFC","#B684FF","#D575FE","#F562DD","#F561C6"),labels=my.labels_tax_final_10)
taxonomy_final_10

## Figure S6
## Wild Only
my.labels_tax_final_wild <- c("Accipitriformes (n=3)","Anseriformes (n=10)","Bucerotiformes (n=1)","Caprimulgiformes (n=1)","Charadriiformes (n=49)","Ciconiiformes (n=3)","Columbiformes (n=1)","Coraciiformes (n=2)","Cuculiformes (n=2)","Eurypygiformes (n=1)","Falconiformes (n=4)","Galliformes (n=15)","Gruiformes (n=11)","Passeriformes (n=120)","Pelecaniformes (n=8)","Piciformes (n=4)","Podicipediformes (n=1)","Procellariiformes (n=9)","Psittaciformes (n=5)","Pterocliformes (n=1)","Sphenisciformes (n=10)","Struthioniformes (n=4)")
taxonomy_plot_final_wild <- data_final[data_final$Wild.Managed.Captive=="Wild",]  %>%
  mutate(Order = fct_relevel(Order,"Accipitriformes","Anseriformes","Bucerotiformes","Caprimulgiformes","Charadriiformes","Ciconiiformes","Columbiformes","Coraciiformes","Cuculiformes","Eurypygiformes","Falconiformes","Galliformes","Gruiformes","Passeriformes","Pelecaniformes","Piciformes","Podicipediformes","Procellariiformes","Psittaciformes","Pterocliformes","Sphenisciformes","Struthioniformes")) %>%
  ggplot(aes(x=Order, y=Mean.Hatching.Failure.std, fill=Order)) + 
  geom_boxplot() + 
  labs(title="Plot of Mean Hatching Failure % according to Taxonomic Order\nFinal Dataset (Wild Only | n = 265)", x="Order", y="Mean Hatching Failure (%)")
taxonomy_final_wild <- taxonomy_plot_final_wild + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=90,hjust=0.95)) + 
  scale_x_discrete(label=my.labels_tax_final_wild) + 
  scale_fill_discrete(name="Order",labels=my.labels_tax_final_wild) + 
  scale_fill_manual(values=c("#F7766D","#EE8043","#E18A00","#D19300","#A8A401","#8CAB02","#68B100","#3EB700","#3FBB49","#3FBE70","#3FC090","#3FC2AB","#3DBFC4","#37B5ED","#32ACFC","#42A0FF","#8B93FF","#B684FF","#D575FE","#EB69F0","#F562DD","#F665AC"),labels=my.labels_tax_final_wild)
taxonomy_final_wild

## Wild Managed Only
my.labels_tax_final_wildmanaged <- c("Anseriformes (n=10)","Bucerotiformes (n=1)","Cathartiformes (n=3)","Charadriiformes (n=9)","Columbiformes (n=2)","Coraciiformes (n=1)","Falconiformes (n=9)","Galliformes (n=2)","Gruiformes (n=1)","Passeriformes (n=111)","Pelecaniformes (n=3)","Procellariiformes (n=1)","Psittaciformes (n=5)","Pterocliformes (n=1)","Strigiformes (n=11)","Struthioniformes (n=1)","Suliformes (n=1)")
taxonomy_plot_final_wildmanaged <- data_final[data_final$Wild.Managed.Captive=="ManagedWild",] %>%
  mutate(Order = fct_relevel(Order,"Anseriformes","Bucerotiformes","Cathartiformes","Charadriiformes","Columbiformes","Coraciiformes","Falconiformes","Galliformes","Gruiformes","Passeriformes","Pelecaniformes","Procellariiformes","Psittaciformes","Pterocliformes","Strigiformes","Struthioniformes","Suliformes")) %>%
  ggplot(aes(x=Order, y=Mean.Hatching.Failure.std, fill=Order)) + 
  geom_boxplot() +  
  labs(title="Plot of Mean Failure % according to Taxonomic Order\nFinal Dataset (Wild Managed Only | n = 172)", x="Order", y="Mean Hatching Failure (%)")
taxonomy_final_wildmanaged <- taxonomy_plot_final_wildmanaged + 
  theme_classic() +  
  theme(axis.text.x=element_text(angle=90,hjust=0.95)) + 
  scale_x_discrete(label=my.labels_tax_final_wildmanaged) + 
  scale_fill_discrete(name="Order",labels=my.labels_tax_final_wildmanaged) +  
  scale_fill_manual(values=c("#EE8043","#E18A00","#BE9C00","#A8A401","#68B100","#3EB700","#3FC090","#3FC2AB","#3DBFC4","#37B5ED","#32ACFC","#B684FF","#D575FE","#EB69F0","#F561C6","#F665AC","#F66C91"),labels=my.labels_tax_final_wildmanaged)
taxonomy_final_wildmanaged

data_final[data_final$Wild.Managed.Captive=="Captive",]  %>%
  count(Order)

## Captive Only
my.labels_tax_final_captive <- c("Anseriformes (n=4)","Cathartiformes (n=3)","Falconiformes (n=1)","Galliformes (n=8)","Gruiformes (n=1)","Otidiformes (n=2)","Passeriformes (n=14)","Pelecaniformes (n=1)","Psittaciformes (n=2)","Pterocliformes (n=4)","Sphenisciformes (n=1)")
taxonomy_plot_final_captive <- data_final[data_final$Wild.Managed.Captive=="Captive",] %>%
  mutate(Order = fct_relevel(Order,"Anseriformes","Cathartiformes","Falconiformes","Galliformes","Gruiformes","Otidiformes","Passeriformes","Pelecaniformes","Psittaciformes","Pterocliformes","Sphenisciformes")) %>%
  ggplot(aes(x=Order, y=Mean.Hatching.Failure.std, fill=Order)) + 
  geom_boxplot() +  
  labs(title="Plot of Mean Failure % according to Taxonomic Order\nFinal Dataset (Captive Only | n = 41)", x="Order", y="Mean Hatching Failure (%)")
taxonomy_final_captive <- taxonomy_plot_final_captive + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=90,hjust=0.95)) + 
  scale_x_discrete(label=my.labels_tax_final_captive) + 
  scale_fill_discrete(name="Order",labels=my.labels_tax_final_captive) +  
  scale_fill_manual(values=c("#EE8043","#BE9C00","#3FC090","#3FC2AB","#3DBFC4","#3ABBDA","#37B5ED","#32ACFC","#D575FE","#EB69F0","#F562DD"),labels=my.labels_tax_final_captive) 
taxonomy_final_captive

