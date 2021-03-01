rm(list=ls())

#### Load libraries and data ####
library(tidyverse)
library(ggplot2)


load("ResultsNicheSize.RData")
load("ResultsNicheFill.RData")

#### Merge niche size and niche fill data ####

# create columns to save niche fill data
dati$fill.mean<-NA
dati$fill.sd<-NA

# calculate mean and sd for niche fill
for (i in dati$name) {
  dati$fill.mean[dati$name==i] <- (mean(niche.fill[[i]][["share.included"]]))*100
  dati$fill.sd[dati$name==i] <- (sd(niche.fill[[i]][["share.included"]]))*100
}

head(dati)

# PLOTS ####
# Plot niche size data #
ggplot(dati, aes(dispersal, volume, fill=dispersal)) +
  geom_boxplot() + theme_bw() +
  # geom_dotplot(binaxis='y', 
  #              stackdir='center', 
  #              dotsize = .5) +
  scale_fill_viridis_d(
    name="Dispersal type",
    breaks=c("ANE", "END", "EPI", "HYD"),
    labels=c("Anemochorous", "Endozoochorous", "Epizoochorous", "Hydrochorous")) +
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12)) +
  ylab("Mean hypervolume size") + theme(axis.text=element_text(size=12)) +  xlab( "Dispersal type") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) 

# Plot niche fill data #
ggplot(dati, aes(dispersal, fill.mean, fill=dispersal)) +
  geom_boxplot() + theme_bw() +
  # geom_dotplot(binaxis='y', 
  #              stackdir='center', 
  #              dotsize = .5) +
  scale_fill_viridis_d(
    name="Dispersal type",
    breaks=c("ANE", "END", "EPI", "HYD"),
    labels=c("Anemochorous", "Endozoochorous", "Epizoochorous", "Hydrochorous")) +
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size = 12)) +
  ylab("Niche fill percentage") + theme(axis.text=element_text(size=12)) +  xlab( "Dispersal type") +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank()) 

#### GLM ####
# Niche size #

# GLM
glm.volume <-glm(volume~dispersal,family="gaussian",data=dati)
summary(glm.volume)

# AOV
aov.volume<-aov(volume~dispersal,data=dati)
summary(aov.volume)

# Post Hoc analysis
library(agricolae)
res2<-HSD.test(aov.volume, "dispersal", group=TRUE)

# Niche fill #

# GLM
glm.fill <- glm(fill.mean~dispersal, family=gaussian, data=dati)
summary(glm.fill)

# AOV 
aov.fill<-aov(fill.mean~dispersal,data=dati)
summary(aov.fill)

# Post Hoc analysis
library(agricolae)
res<-HSD.test(aov.fill, "dispersal", group=TRUE)

# count sp by dispersal syndrome
dati %>% 
  group_by(dispersal) %>%
  summarise(no_rows = length(dispersal))

#### PGLS ####

# load libraries
library(ape)
library(geiger)
library(nlme)
library(phytools)

# load tree
tree <- ape::read.tree("/tree.nex")

# in case names in the tree have "_" and in the database have " "
dati$name <- gsub(" ", "_", dati$name, fixed=TRUE)

# check whether the names in the database and the tree are the same
name.check(phy = tree, data.names = dati$name)

# correct names of sp in dati 3 so they match the tree, repeat
# until $data_not_tree is character(0)

dati$name[dati$name=="Corema_album"] <- "Corema_album_ssp_azoricum"
dati$name[dati$name=="Ilex_perado"] <- "Ilex_perado_ssp_azorica"
dati$name[dati$name=="Lysimachia_tenella"] <- "Anagallis_tenella"
dati$name[dati$name=="Morella_faya"] <- "Myrica_faya"
dati$name[dati$name=="Platanthera_pollostantha"] <- "Platanthera_azorica"

a <- name.check(phy = tree, data.names = dati$name)

# now we have to prune the tree because pgls function needs the same number of 
# observations in the db and tips in the tree

prune.sp <- a$tree_not_data
pruned.tree <- drop.tip(phy = tree, tip = tree$tip.label[match(prune.sp, tree$tip.label)])

# check if the tree is rooted
is.rooted(pruned.tree)

# next is to set rownames in data frame so they match tree tip names
rownames(dati) <- dati$name

# PGLS calculated with caper package #
library(caper)

# compute lambda for discrete data
dispersal <- dati$dispersal
names(dispersal) <- dati$name

lambda.est <- fitDiscrete(pruned.tree, dispersal)
lambda.est

# pgls requieres to combine our tree with the dataset with comparative.data()
comp.data<-comparative.data(pruned.tree, dati, names.col="name", vcv.dim=3, warn.dropped=TRUE)

# pgls for niche size
gls.volume<-pgls(volume~dispersal, data=comp.data, lambda = "ML") # we can do the same without lambda
summary(pgls.volume)

# pgls for niche fill
pgls.fill<-pgls(fill.mean~dispersal, data=comp.data, lambda = "ML") # we can do the same without lambda
summary(pgls.fill)

# optional: plot likelihood profile for lambda in niche size model
lm.lk.volume<-pgls.profile(pgls.volume, which="lambda")
plot(lm.lk.volume) # 0 is the most likely value for lambda

# optional plot likelihood profile for lambda in niche fill model
lm.lk.fill<-pgls.profile(pgls.fill, which="lambda")
plot(lm.lk.fill)  # 0 is the most likely value for lambda
