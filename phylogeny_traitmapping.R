# ==============================================================================
# Mapping Sociality and Nesting on a Circular Tree
# ==============================================================================
# Plots a circular phylogeny with outer rings showing sociality and 
# nesting traits of each species tip.
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: clear environment, set working directory, load libraries
# ------------------------------------------------------------------------------

rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(ape)
library(phytools)

# If getting the getYmult() error:
# # Install devtools if not already installed
# install.packages("devtools")
# 
# # Reinstall phytools from GitHub
# devtools::install_github("liamrevell/phytools")


# ------------------------------------------------------------------------------
# Load trait data and phylogenetic tree
# ------------------------------------------------------------------------------

traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")


# ------------------------------------------------------------------------------
# Define coordinate conversion functions (polar <-> cartesian)
# ------------------------------------------------------------------------------

toCart <- function(r, th, deg = FALSE) {
  if (deg) th = th * pi / 180
  x = r * cos(th)
  y = r * sin(th)
  return(list(x = x, y = y))
}

toPolar <- function(x, y) {
  r = sqrt(x^2 + y^2)
  th = atan2(y, x)
  return(list(r = r, th = th))
}


# ------------------------------------------------------------------------------
# Prepare trait vectors and align with tree tip labels
# ------------------------------------------------------------------------------

soci_x = traits$sociality_binary
nest_x = traits$nest_binary

names(soci_x) <- traits$tips
names(nest_x) <- traits$tips

# Sanity check: ensure trait names match tree tip labels
identical(sort(names(soci_x)), sort(tree$tip.label))  # all the same
identical(sort(names(nest_x)), sort(tree$tip.label))  # all the same


# ------------------------------------------------------------------------------
# Set color palettes for each trait
# ------------------------------------------------------------------------------

soci_cols = c("slategray1", "steelblue")
names(soci_cols) <- c("solitary", "social")

nest_cols = c("thistle1", "violet")
names(nest_cols) <- c("ground", "aboveground")

N = length(tree$tip.label) # Number of terminal tips


# ------------------------------------------------------------------------------
# Begin plotting
# ------------------------------------------------------------------------------

pdf(file = "plots/phylogeny_traitmapping.pdf", width = 7, height = 7)
plotTree(tree, ftype = "off", lwd = 0.3, type = "fan") # Draw the circular tree without tip labels


# ------------------------------------------------------------------------------
# Extract tip coordinates from tree
# ------------------------------------------------------------------------------

obj <- get("last_plot.phylo",envir=.PlotPhyloEnv)  
  # note that this will extract coordinates for whatever tree you just plotted (in our case, the tree plotted on line 137)
xx = obj$xx # includes all nodes
tips_xx = xx[1:N] # subsets to only terminal nodes
yy = obj$yy
tips_yy = yy[1:N] 


# ------------------------------------------------------------------------------
# Define parameters for ring layout
# ri = radial increment for each ring
# len = number of rings per trait
# space = spacing between trait rings
# ------------------------------------------------------------------------------

ri = 0.2
len = 10
space = 4

# Save original polar coordinates for reuse in rings and tip labels
p_original <- toPolar(tips_xx, tips_yy)


# ------------------------------------------------------------------------------
# Draw trait rings around the phylogeny
# First empty space, then colored ring for sociality
# Then more space, then colored ring for nesting
# ------------------------------------------------------------------------------

for (i in 0:(len*2 + (space*2))) {
  
  p = toPolar(tips_xx, tips_yy)
  p$r = p$r + ri*i # Expand radius outward for each ring
  
  c = toCart(p$r, p$th, deg = F) # Convert back to Cartesian
  
  if (i < space) {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.1) # inner blank space
  } else if (i >= space & i < space+len) {
    points(c$x, c$y, col = soci_cols[soci_x[tree$tip.label]], pch = 16, cex = 0.1) # sociality ring
  } else if (i >= space+len & i < len+space*2) {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.1) # middle blank space
  } else {
    points(c$x, c$y, col = nest_cols[nest_x[tree$tip.label]], pch = 16, cex = 0.1) # nesting ring
  }
}


# ------------------------------------------------------------------------------
# Optional: Add species names outside outer ring
# Useful for manually checking if the traits mapped correctly
# ------------------------------------------------------------------------------

# label_radius <- max(p_original$r) + ri * (len*2 + space*2 + 1)
# label_coords <- toCart(label_radius, p_original$th)

# for (i in seq(1, N, by = 20)) {
#   angle_deg <- p_original$th[i] * 180 / pi
#   align <- ifelse(angle_deg > 90 & angle_deg < 270, 1, 0)
#   rot_angle <- ifelse(align == 1, angle_deg + 180, angle_deg)
#   text(label_coords$x[i], label_coords$y[i], labels = tree$tip.label[i],
#        cex = 0.3, srt = rot_angle, adj = align)
# }


# ------------------------------------------------------------------------------
# Add legends for both traits
# ------------------------------------------------------------------------------

legend("topleft", legend = c("Solitary", "Social"), col = soci_cols, bty = "n", pch = 16, cex = 1.2)
legend("topright", legend = c("Ground", "Above-ground"), col = nest_cols, bty = "n", pch = 16, cex = 1.2)

dev.off() # Close PDF output
