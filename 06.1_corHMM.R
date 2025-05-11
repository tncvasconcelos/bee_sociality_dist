# rm(list=ls())
# setwd("/Users/tvasc/Desktop/bee_sociality_dist")
library(corHMM)

source("00_utility_functions.R")
#--------------------------------------
# First organizing dataset:
# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")
phy <- keep.tip(phy, which(phy$tip.label %in% traits$tips))

#------------------------------------
dat <- traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]
dat <- dat[,c("tips","sociality_binary","nest_binary")]

#-------------------------------------
# determining best corHMM models with corhmm dredge
# Note: we are fixing the root to be state 4 ("solitary/ground") with root.p=c(0,0,0,1)
# dredge_sociality <- corHMM:::corHMMDredge(phy, dat, max.rate.cat=2, root.p=c(0,0,0,1))
#save(dredge_sociality, file="corHMMdredge_results/corhmm_dredge_binary.Rsave")
load("corHMMdredge_results/corhmm_dredge_binary.Rsave")

corhmm_tbl_sociality <- corHMM:::getModelTable(dredge_sociality)
write.csv(corhmm_tbl_sociality, file="corHMMdredge_results/corhmm_tbl_dredge.csv")

#---------------------
# Rate matrix figure

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(scales)

# Define states
states <- c(
  "R1 social|aboveground", "R1 solitary|aboveground", "R1 social|ground", "R1 solitary|ground",
  "R2 social|aboveground", "R2 solitary|aboveground", "R2 social|ground", "R2 solitary|ground"
)

rates_mat <- dredge_sociality[[10]]$solution

# Rebuild the transition matrix manually (you pasted it but it's not in R object form)
# So let's manually create it in a simple way for now:

# Create a named 8x8 matrix
rates_mat <- matrix(NA, nrow = 8, ncol = 8, dimnames = list(states, states))

# Convert matrix to edge list
edges <- as.data.frame(as.table(rates_mat)) %>%
  filter(!is.na(Freq)) %>%
  rename(from = Var1, to = Var2, rate = Freq)

# Define biological colors
state_colors <- c(
  "social|aboveground" = "darkblue",
  "solitary|aboveground" = "darkred",
  "social|ground" = "darkgreen",
  "solitary|ground" = "purple"
)

# Assign colors to nodes
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

nodes <- tibble(name = states) %>%
  mutate(color = sapply(name, assign_color))

# Build graph
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
tg <- as_tbl_graph(g)

# Plot
pdf("plots/transition_matrix_corhmm.pdf")
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

#---------------------

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


anc_recon <- dredge_sociality[[10]]$states
#anc_recon[1,] # checking root state


#boxplot(dredge_sociality[[10]]$states)

dev.off()

plotRECON <- function(phy, likelihoods, piecolors=NULL, cex=0.5, pie.cex=0.25, file=NULL, height=11, width=8.5, show.tip.label=TRUE, title=NULL, ...){
  if(is.null(piecolors)){
    #piecolors=c("pink","black","red","yellow","forestgreen","blue","coral","aquamarine")
    piecolors=c("#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF", 
                "#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF")
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

#pal1 <- hcl.colors(8, palette = "Viridis", alpha = 0.7)
#pal1 <- hcl.colors(4, palette = "Inferno", alpha = 1)
#pal2 <- hcl.colors(4, palette = "Inferno", alpha = 0.25)

# 
# custom_colors <- c(pal1, pal2)
# names(custom_colors) <- colnames(anc_recon)

pdf("corHMMdredge_results/corhmm_dredge_recon_final.pdf", height=45, width=10)

plotRECON(
  phy = phy,
  likelihoods = anc_recon,
  pie.cex = 0.3,  # Size of pie charts
  show.tip.label = T,
  cex=0.1
)

axisPhylo()
dev.off()

