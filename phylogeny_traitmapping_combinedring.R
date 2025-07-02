# ==============================================================================
# Mapping Combined Sociality and Nesting on a Circular Tree
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------
rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(ape)
library(phytools)
library(viridis)
library(dplyr)

# ------------------------------------------------------------------------------
# Load tree and traits
# ------------------------------------------------------------------------------
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")

# ------------------------------------------------------------------------------
# Create combined trait column and clean labels
# ------------------------------------------------------------------------------
traits$trait_combo <- paste(traits$sociality_binary, traits$nest_binary, sep = "/")
traits$trait_combo <- recode(traits$trait_combo,
                             "solitary/ground" = "Solitary/Ground",
                             "solitary/aboveground" = "Solitary/Above-ground",
                             "social/ground" = "Social/Ground",
                             "social/aboveground" = "Social/Above-ground")

# Reorder trait dataframe to match tree tip label order
traits <- traits[match(tree$tip.label, traits$tips), ]
stopifnot(identical(traits$tips, tree$tip.label))  # Sanity check

# Create named vector for trait combinations
trait_combo_vector <- setNames(traits$trait_combo, traits$tips)

# ------------------------------------------------------------------------------
# Color palette for combined trait combo using viridis::cividis
# ------------------------------------------------------------------------------
combo_levels <- c("Solitary/Ground", "Solitary/Above-ground", "Social/Ground", "Social/Above-ground")
combo_colors <- viridis::viridis(n = 4, option = "cividis")
names(combo_colors) <- combo_levels

# ------------------------------------------------------------------------------
# Coordinate conversion functions
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
# Plot tree
# ------------------------------------------------------------------------------
N <- length(tree$tip.label)
png(filename = "plots/phylogeny_traitmapping_combinedring.png",
    width = 1800, height = 1800, res = 300)
plotTree(tree, ftype = "off", lwd = 0.3, type = "fan")

# Extract tip coordinates
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips_xx <- obj$xx[1:N]
tips_yy <- obj$yy[1:N]
p_original <- toPolar(tips_xx, tips_yy)

# ------------------------------------------------------------------------------
# Define ring layout parameters
# ------------------------------------------------------------------------------
ri <- 0.2    # radial increment
len <- 15    # ring thickness in "pixels"
space <- 4   # blank space before and after

# ------------------------------------------------------------------------------
# Sanity check: ensure all trait combos have a color
# ------------------------------------------------------------------------------
stopifnot(all(trait_combo_vector %in% names(combo_colors)))

# ------------------------------------------------------------------------------
# Draw the single trait ring
# ------------------------------------------------------------------------------
for (i in 0:(len + 2 * space)) {
  p <- toPolar(tips_xx, tips_yy)
  p$r <- p$r + ri * i
  c <- toCart(p$r, p$th)
  
  if (i < space) {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.15)
  } else if (i >= space & i < space + len) {
    trait_colors <- combo_colors[trait_combo_vector[tree$tip.label]]
    points(c$x, c$y, col = trait_colors, pch = 16, cex = 0.25)
  } else {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.15)
  }
}

# ------------------------------------------------------------------------------
# Add legend
# ------------------------------------------------------------------------------
legend("bottomright",
       legend = names(combo_colors),
       col = combo_colors,
       pch = 16,
       pt.cex = 0.9, # dot size
       title = "Trait Combination",
       bty = "n",
       cex = 0.6, # text size
       text.font = 1,
       x.intersp = 1.3, # spacing
       inset = c(0, 0.025),  # move legend outward
       xpd = TRUE)            # allow plotting outside plot region


dev.off()  # Close the PDF device
