# Setup
#rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(phytools)
library(ape)

# Loading traits and tree
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")

# Loading functions
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

# Traits subset
soci_x = traits$sociality_binary
nest_x = traits$nest_binary

names(soci_x) <- traits$tips
names(nest_x) <- traits$tips

identical(sort(names(soci_x)), sort(tree$tip.label))  # all the same
identical(sort(names(nest_x)), sort(tree$tip.label))  # all the same

# plotting with phytools

# set color palettes
soci_cols = c("slategray1", "steelblue")
names(soci_cols) <- c("solitary", "social")

nest_cols = c("thistle1", "violet")
names(nest_cols) <- c("ground", "aboveground")

N = length(tree$tip.label)

# START HERE FOR PLOTTING

pdf(file = "traitmapped_phylogeny.pdf", width = 7, height = 7)

plotTree(tree, ftype = "off", lwd = 0.3, type = "fan")

# extracting location of tips
obj <- get("last_plot.phylo",envir=.PlotPhyloEnv)  # please note this will extract coordinates for whatever tree you just plotted (in our case, the tree plotted on line 137)
xx = obj$xx   # includes all nodes
tips_xx = xx[1:N]   # subsets to only terminal nodes
yy = obj$yy
tips_yy = yy[1:N] 

# extension

ri = 0.2
len = 10
space = 4

# Save original polar coordinates once
p_original <- toPolar(tips_xx, tips_yy)


for (i in 0:(len*2 + (space*2))) {
  
  p = toPolar(tips_xx, tips_yy)
  p$r = p$r + ri*i    # increments small value to make dot go out a little bit
  
  c = toCart(p$r, p$th, deg = F)   # conversion back into cartesian space
  
  if (i < space) {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.1)
  } else if (i >= space & i < space+len) {
    points(c$x, c$y, col = soci_cols[soci_x[tree$tip.label]], pch = 16, cex = 0.1)
  } else if (i >= space+len & i < len+space*2) {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.1)
  } else {
    points(c$x, c$y, col = nest_cols[nest_x[tree$tip.label]], pch = 16, cex = 0.1)
  }
}

##############################################################################
# Add species labels outside the last ring
#label_radius = max(p$r) + 0.3  # adjust distance as needed
#label_coords = toCart(label_radius, p$th)

label_radius = max(p_original$r) + ri * (len*2 + space*2 + 1)  # far enough out
label_coords = toCart(label_radius, p_original$th)


# Add tip labels
for (i in seq(1, N, by = 20)) {
  angle_deg <- p_original$th[i] * 180 / pi
  align <- ifelse(angle_deg > 90 & angle_deg < 270, 1, 0)
  rot_angle <- ifelse(align == 1, angle_deg + 180, angle_deg)
  
  text(label_coords$x[i], label_coords$y[i], labels = tree$tip.label[i],
       cex = 0.3, srt = rot_angle, adj = align)
}


##############################################################################

legend("topleft", legend = c("Solitary", "Social"), col = soci_cols, bty = "n", pch = 16, cex = 1.2)
legend("topright", legend = c("Ground", "Above-ground"), col = nest_cols, bty = "n", pch = 16, cex = 1.2)

dev.off()
