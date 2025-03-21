# rm(list=ls())
library(corHMM)
library(OUwie)
library(parallel)
# setwd("/Users/tvasc/Desktop/bee_sociality_dist")

source("00_utility_functions.R")
#--------------------------------------
# First organizing dataset:
# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv", full.names = T)

# Let's take bio1 and bio12 (temperature and precipitation) and seasonalities (bio4 and bio15)
all_climatic_vars <- all_climatic_vars[grep(paste(c("bio_1_","bio_12_","bio_4_",
                                                    "bio_15_","ai"),collapse="|"), all_climatic_vars)]
climatic_list <- lapply(all_climatic_vars, read.csv)

# Now merge everything in one table
merged_climatic_vars <- climatic_list[[1]] 
for(i in 2:length(climatic_list)) {
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(merged_climatic_vars, one_climatic_var, by="species") 
}

# Select only mean columns
merged_climatic_vars <- merged_climatic_vars[,c(1, grep("mean", colnames(merged_climatic_vars)))]

# And finally merge to the trait data
merged_traits <- merge(traits, merged_climatic_vars, by.x="tips",by.y="species")

# merged_traits0 <- merged_traits[1,]
# while(length(unique(merged_traits0$sociality))<3) {
#   merged_traits0 <- merged_traits[sample(1:nrow(merged_traits), 100),]
# }
# merged_traits <- merged_traits0

# log all continuous variables:
merged_traits$mean_bio_1 <- log((merged_traits$mean_bio_1)+273) # transform celcius to kelvin for temperature
merged_traits$mean_bio_12 <- log(merged_traits$mean_bio_12)
merged_traits$mean_bio_15 <- log(merged_traits$mean_bio_15)
merged_traits$mean_bio_4 <- log(merged_traits$mean_bio_4)
merged_traits$mean_awi_pm_sr_yr <- log(merged_traits$mean_awi_pm_sr_yr)

merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_1))
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_4))
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_awi_pm_sr_yr))

phy <- keep.tip(phy, which(phy$tip.label %in% merged_traits$tips))

#------------------------------------
dat <- merged_traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]
dat <- dat[,c("tips","sociality_binary","nest_binary")]

corhmm_fits <- corHMM:::fitCorrelationTest(phy, dat) 
save(corhmm_fits, file = "corhmm_fits_cortest.Rsave")
#load("corhmm_fits_cortest.Rsave")
corhmm_tbl <- corHMM:::getModelTable(corhmm_fits)

# LIKELIHOOD RATIO TEST
teststat <- -2 * (corhmm_tbl$lnLik[2] - corhmm_tbl$lnLik[4])
p.val <- pchisq(teststat, df = 8, lower.tail = FALSE)

# determining corHMM models with corhmm dredge
dredge_sociality <- corHMM:::corHMMDredge(phy, dat, max.rate.cat=2)
?corHMM
#dredge_sociality[[8]]
save(dredge_sociality, file="corhmm_dredge_binary.Rsave")
load("corhmm_dredge_binary.Rsave")
corhmm_tbl_sociality <- corHMM:::getModelTable(dredge_sociality)
write.csv(corhmm_tbl_sociality, file="corhmm_tbl_dredge.csv")
# 
# dredge_nesting <- corHMM:::corHMMDredge(phy, merged_traits[,c("tips","nest")],max.rate.cat=3)
# save(dredge_nesting, file="corhmm_dredge_nesting.Rsave")
# corhmm_tbl_nesting <- corHMM:::getModelTable(dredge_nesting)
# write.csv(corhmm_tbl_nesting, file="corhmm_tbl_nesting.csv")


# Ancestral state plot

# this script will reconstruc the ancestral state based on our model averaged corHMM models

## functions 
WtQ <- function(Q, Weights){
  Weights <- Weights[!is.na(Q)]/sum(Weights[!is.na(Q)])
  AvgQ <- sum(Q[!is.na(Q)] * Weights)
  return(AvgQ)
}


# model average rates from corHMM
getModelAvgRate <- function(file){
  load(file)
  rate.mat <- obj$res_ARD.ARD.2$index.mat
  AICc <- unlist(lapply(obj, function(x) x$AICc))
  AICwt <- exp(-0.5 * AICc - min(AICc))/sum(exp(-0.5 * AICc - min(AICc)))
  # Solutions <- lapply(obj, function(x) x$solution)
  Solutions <- lapply(obj, function(x) c(na.omit(c(x$solution))))
  Solutions[[1]] <-  c(Solutions[[1]][1], NA, Solutions[[1]][2], NA,
                       NA, Solutions[[1]][1], NA, Solutions[[1]][2])
  Solutions[[2]] <-  c(Solutions[[2]][1], NA, Solutions[[2]][2], NA,
                       NA, Solutions[[2]][1], NA, Solutions[[2]][2])
  Rates <- do.call(rbind, Solutions)
  p.wt <- apply(Rates, 2, function(x) WtQ(x, AICwt))
  rate.mat[!is.na(rate.mat)] <- p.wt 
  return(rate.mat)
}


getTipRecon <- function(file){
  load(file)
  phy <- obj$res_ER$phy
  data <- obj$res_ER$data
  root.p <- obj$res_ARD$root.p
  index.mat <- obj$res_ARD.ARD.2$index.mat
  p <- getModelAvgRate(file)[sapply(1:max(index.mat, na.rm = TRUE), function(x) match(x, index.mat))]
  res <- corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = index.mat, node.states = "marginal", p = p, root.p = root.p, get.tip.states = TRUE)
  return(res)
}

getTipRecon()

getModelAvgRate(corhmm_fits$hidden_Markov_correlated_model_fit)

anc_recon <- corhmm_fits$hidden_Markov_correlated_model_fit$states

boxplot(corhmm_fits$hidden_Markov_correlated_model_fit$states)

dev.off()

plotRECON <- function(phy, likelihoods, piecolors=NULL, cex=0.5, pie.cex=0.25, file=NULL, height=11, width=8.5, show.tip.label=TRUE, title=NULL, ...){
  if(is.null(piecolors)){
    piecolors=c("pink","black","red","yellow","forestgreen","blue","coral","aquamarine","darkorchid","gold","grey","yellow","#3288BD","#E31A1C")
  }
  if(!is.null(file)){
    pdf(file, height=height, width=width,useDingbats=FALSE)
  }
  plot(phy, cex=cex, show.tip.label=show.tip.label, ...)
  
  if(!is.null(title)){
    title(main=title)
  }
  nodelabels(pie=likelihoods,piecol=piecolors, cex=pie.cex)
  states <- colnames(likelihoods)
  #legend(x="topleft", states, cex=0.8, pt.bg=piecolors,col="black",pch=21);
  
  if(!is.null(file)){
    dev.off()
  }
}

pal1 <- hcl.colors(8, palette = "Viridis", alpha = 0.7)
# pal2 <- hcl.colors(4, palette = "Inferno", alpha = 0.7)
# 
# custom_colors <- c(pal1, pal2)
# names(custom_colors) <- colnames(anc_recon)

pdf("corhmm_recon_test_no_dredge.pdf", height=45, width=10)
plotRECON(
  phy = phy,
  likelihoods = anc_recon,
  pie.cex = 0.3,  # Size of pie charts
  show.tip.label = T,
  cex=0.1
)
axisPhylo()
dev.off()


?plotRECON
