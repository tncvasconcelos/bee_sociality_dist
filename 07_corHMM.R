# ==============================================================================
# 06.1 Ancestral state reconstruction of sociality and nesting strategy
# ==============================================================================
# Uses corHMM to model the correlated evolution of sociality and nesting strategy 
# in bees under hidden rate models, selects the best-fitting model via AIC, and 
# reconstructs ancestral states across the phylogeny. Visualizes the inferred 
# transition matrix and ancestral state reconstructions. 
# ==============================================================================

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
rm(list=ls())
wd <- "/Users/lenarh/Desktop/bee_sociality_dist"
setwd(wd)
results_wd <- file.path(wd, "results/corHMM/corHMMdredge_results")

library(corHMM)
library(phytools)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(scales)
#devtools::install_github("thej022214/corHMM")

# Load utility functions
source("00_utility_functions.R")

# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")

#-------------------------------------------------------------------------------
# Organizing dataset
#-------------------------------------------------------------------------------
# Drop tips in the tree that don’t match any in the trait data
phy <- keep.tip(phy, which(phy$tip.label %in% traits$tips))

# Find species shared by both datasets
dat <- traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

# Reorder the data and tree to align species
dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]

# Keep relevant trait columns
dat <- dat[,c("tips","sociality_binary","nest_binary")]

#-------------------------------------------------------------------------------
# Model selection (run or load corHMM dredge) with root fixed as solitary/ground
#-------------------------------------------------------------------------------
# Run model selection with corHMMDredge (if not already done)
#dredge_sociality <- corHMM::corHMMDredge(phy, dat, max.rate.cat=2, root.p=c(0,0,0,1))
#save(dredge_sociality, file="corHMMdredge_results/corhmm_dredge_binary_Oct5.Rsave")

# Load saved dredge results
load("results/corHMM/corHMMdredge_results/corhmm_dredge_binary.Rsave") # old
load("corHMMdredge_results/corhmm_dredge_binary_Oct5.Rsave") # new

# Save summary table to compare models with AIC values
corhmm_tbl_sociality <- corHMM:::getModelTable(dredge_sociality)
#write.csv(corhmm_tbl_sociality, file="corHMMdredge_results/corhmm_tbl_dredge.csv")

#-------------------------------------------------------------------------------
# Visualize transition matrix
#-------------------------------------------------------------------------------
# Define states
states <- c(
  "R1 social|aboveground", "R1 solitary|aboveground", "R1 social|ground", "R1 solitary|ground",
  "R2 social|aboveground", "R2 solitary|aboveground", "R2 social|ground", "R2 solitary|ground"
)

# Extract the transition rate matrix from the 10th (best) model
rates_mat <- dredge_sociality[[21]]$solution

# Manually create transition matrix
# Initialize an 8x8 empty matrix for states
rates_mat <- matrix(NA, nrow = 8, ncol = 8, dimnames = list(states, states)) 

# Convert matrix to edge list
edges <- as.data.frame(as.table(rates_mat)) %>% 
  filter(!is.na(Freq)) %>%
  rename(from = Var1, to = Var2, rate = Freq)

# Define state colors
state_colors <- c( 
  "social|aboveground" = "darkblue",
  "solitary|aboveground" = "darkred",
  "social|ground" = "darkgreen",
  "solitary|ground" = "purple"
)

# Assign colors to nodes depending on rate class (R1 = opaque, R2 = transparent)
assign_color <- function(state) {
  part <- sub("R[12] ", "", state)
  rate_class <- substr(state, 2, 2)  # R1 or R2
  
  base_color <- state_colors[[part]]
  
  if (rate_class == "1") {
    return(base_color)
  } else {
    return(alpha(base_color, 0.4))  # lighter for R2
  }
}

nodes <- tibble(name = states) %>% # Creates a node table with colors
  mutate(color = sapply(name, assign_color))

# Create graph and plot
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
tg <- as_tbl_graph(g)

pdf(file.path(results_wd, "transition_matrix_corhmm.pdf"))

ggraph(tg, layout = 'circle') +
  geom_edge_link(aes(width = rate),
                 arrow = arrow(length = unit(3, 'mm')),
                 end_cap = circle(3, 'mm'),
                 edge_alpha = 0.8,
                 color = "gray40") +
  geom_node_point(aes(color = name), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_edge_width(range = c(0.5, 5)) +
  scale_color_manual(values = nodes$color) +
  theme_void() +
  ggtitle("Transition Diagram for Sociality (Model 10)\nColors by State, Lighter for R2")

dev.off()

#===============================================================================
# Ancestral State Reconstruction (ASR) & Plotting
#===============================================================================
# Reconstruct the ancestral state based on corHMM models
#-------------------------------------------------------------------------------
# Functions 
#-------------------------------------------------------------------------------
WtQ <- function(Q, Weights){ # Weighted average of transition rates using AIC weights
  Weights <- Weights[!is.na(Q)]/sum(Weights[!is.na(Q)])
  AvgQ <- sum(Q[!is.na(Q)] * Weights)
  return(AvgQ)
}

getModelAvgRate <- function(file){ # Loads model set, computes AIC weights, and creates a model-averaged rate matrix
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

getTipRecon <- function(file){ # Reconstructs states using model-averaged transition matrix, returns marginal likelihoods
  load(file)
  phy <- obj$res_ER$phy
  data <- obj$res_ER$data
  root.p <- obj$res_ARD$root.p
  index.mat <- obj$res_ARD.ARD.2$index.mat
  p <- getModelAvgRate(file)[sapply(1:max(index.mat, na.rm = TRUE), function(x) match(x, index.mat))]
  res <- corHMM(phy = phy, data = data, rate.cat = 2, rate.mat = index.mat, node.states = "marginal", p = p, root.p = root.p, get.tip.states = TRUE)
  return(res)
}

# Plots ancestral state reconstructions with colored pie charts at nodes
plotRECON <- function(phy, likelihoods, piecolors=NULL, cex=0.5, pie.cex=0.25, file=NULL, height=11, width=8.5, show.tip.label=TRUE, title=NULL, ...){
  if(is.null(piecolors)){
    #piecolors=c("pink","black","red","yellow","forestgreen","blue","coral","aquamarine")
    # piecolors=c("#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF", 
    #             "#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF")
    piecolors=rev(c("#00204DFF", "#575C6DFF", "#A69D75FF", "#FFEA46FF",
                    "#00204DFF", "#575C6DFF", "#A69D75FF", "#FFEA46FF"))
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
  legend(x="topleft", states, cex=0.8, pt.bg=piecolors,col="black",pch=21);
  
  if(!is.null(file)){
    dev.off()
  }
}

#-------------------------------------------------------------------------------
# Plot ancestral state reconstruction
#-------------------------------------------------------------------------------
# Extract likelihoods for internal nodes from best model (10th)
anc_recon <- dredge_sociality[[10]]$states
#anc_recon[1,] # checking root state

pdf(file.path(results_wd, "corhmm_dredge_recon_final.pdf"), height = 45, width = 10)

plotRECON(
  phy = phy,
  likelihoods = anc_recon,
  pie.cex = 0.3,
  show.tip.label = T,
  cex = 0.1
)

axisPhylo()
dev.off()

#-------------------------------------------------------------------------------
# Estimate number and timing of trait transitions via SIMMAP
# (using corHMM dev version: thej022214/corHMM)
#-------------------------------------------------------------------------------
