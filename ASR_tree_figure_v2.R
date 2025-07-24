# # ============================================================================
# # Circular Phylogeny: ASR Branch Coloring + Trait Ring + Family Labels
# # ============================================================================
# # Circular phylogeny with branches colored by ancestral state reconstruction (ASR), 
# # a trait ring at the tips, and labeled families with MRCA node dots.
# # ==============================================================================
# ============================================================================
# Circular Phylogeny: ASR Branch Coloring + Trait Ring + Time Scale + Family Labels
# ============================================================================

# ------------------------------------------------------------------------------
# Load libraries and data
# ------------------------------------------------------------------------------
library(ape)
library(phytools)
library(viridis)
library(dplyr)

setwd("/Users/lenarh/Desktop/bee_sociality_dist")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
traits <- read.csv("curated_data/bees_traits.csv")

traits <- traits[match(tree$tip.label, traits$tips), ]
traits$trait_combo <- paste(traits$sociality_binary, traits$nest_binary, sep = "/")
traits$trait_combo <- recode(traits$trait_combo,
                             "solitary/ground" = "Solitary/Ground",
                             "solitary/aboveground" = "Solitary/Above-ground",
                             "social/ground" = "Social/Ground",
                             "social/aboveground" = "Social/Above-ground")
trait_combo_vector <- setNames(traits$trait_combo, traits$tips)

combo_colors <- viridis::viridis(n = 4, option = "cividis")
names(combo_colors) <- c("Solitary/Ground", "Solitary/Above-ground",
                         "Social/Ground", "Social/Above-ground")

# ------------------------------------------------------------------------------
# Load ASR from corHMM
# ------------------------------------------------------------------------------
load("corHMMdredge_results/corhmm_dredge_binary.Rsave")
anc_recon <- dredge_sociality[[10]]$states

asr_state_key <- c(
  "R1 social|aboveground"     = "Social/Above-ground",
  "R1 solitary|aboveground"   = "Solitary/Above-ground",
  "R1 social|ground"          = "Social/Ground",
  "R1 solitary|ground"        = "Solitary/Ground",
  "R2 social|aboveground"     = "Social/Above-ground",
  "R2 solitary|aboveground"   = "Solitary/Above-ground",
  "R2 social|ground"          = "Social/Ground",
  "R2 solitary|ground"        = "Solitary/Ground"
)

n_tips <- length(tree$tip.label)
branch_colors <- rep("gray80", nrow(tree$edge))
node_states <- apply(anc_recon, 1, function(x) {
  if (max(x) > 0.5) names(which.max(x)) else NA
})

for (i in 1:nrow(tree$edge)) {
  parent_node <- tree$edge[i, 1]
  if (parent_node > n_tips) {
    state <- node_states[parent_node - n_tips]
    if (!is.na(state)) {
      trait_label <- asr_state_key[state]
      branch_colors[i] <- combo_colors[trait_label]
    }
  }
}

# ------------------------------------------------------------------------------
# Functions: Coordinate conversion and axis drawing
# ------------------------------------------------------------------------------
toCart <- function(r, th, deg = FALSE) {
  if (deg) th <- th * pi / 180
  list(x = r * cos(th), y = r * sin(th))
}
toPolar <- function(x, y) {
  list(r = sqrt(x^2 + y^2), th = atan2(y, x))
}

add_time_axis_radial <- function(tree, ring_times, angle_deg = 0, tick_length = 0.015, label_cex = 0.4) {
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  n_tips <- length(tree$tip.label)
  root_radius <- max(sqrt(obj$xx^2 + obj$yy^2)[(n_tips + 1):length(obj$xx)])
  tree_height <- max(nodeHeights(tree))
  reversed_radii <- (tree_height - ring_times) / tree_height * root_radius
  theta <- angle_deg * pi / 180
  dx <- cos(theta)
  dy <- sin(theta)
  segments(0, 0, root_radius * dx, root_radius * dy, col = "black", lwd = 1)
  
  for (i in seq_along(reversed_radii)) {
    r <- reversed_radii[i]
    x_tick <- r * dx
    y_tick <- r * dy
    perp_dx <- -dy
    perp_dy <- dx
    
    # Draw tick mark
    tick_x0 <- x_tick - tick_length * root_radius * perp_dx
    tick_y0 <- y_tick - tick_length * root_radius * perp_dy
    tick_x1 <- x_tick + tick_length * root_radius * perp_dx
    tick_y1 <- y_tick + tick_length * root_radius * perp_dy
    segments(tick_x0, tick_y0, tick_x1, tick_y1, col = "black", lwd = 1)
    
    # Skip 100 label
    if (ring_times[i] != 100) {
      # Position label above tick
      label_offset <- 2 * tick_length * root_radius
      label_x <- x_tick + label_offset * perp_dx
      label_y <- y_tick + label_offset * perp_dy
      
      label_text <- if (ring_times[i] == 0) "0 Ma" else as.character(ring_times[i])
      text(label_x, label_y, labels = label_text, cex = label_cex, adj = c(0.5, 0))
    }
  }
}
  

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------
png("plots/ASR_phylogeny.png", width = 1800, height = 1800, res = 300)
par(mar = c(2, 2, 2, 2), xpd = TRUE)

# Step 1: Initialize plotting space and coordinate system
plot(tree, type = "fan", edge.color = "transparent", edge.width = 0.5,
     show.tip.label = FALSE, no.margin = FALSE, open.angle = 15, rotate.tree = 7.5)

# Step 2: Re-plot tree with edge colors on top
plot(tree, type = "fan", edge.color = branch_colors, edge.width = 0.5,
     show.tip.label = FALSE, no.margin = FALSE, open.angle = 15, rotate.tree = 7.5, add = TRUE)

# Step 3: Add time axis
ring_times <- c(0, 20, 40, 60, 80, 100)
add_time_axis_radial(tree, ring_times, angle_deg = 0)

# ------------------------------------------------------------------------------
# Trait Ring at Tips
# ------------------------------------------------------------------------------
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips_xx <- obj$xx[1:n_tips]
tips_yy <- obj$yy[1:n_tips]
ri <- 0.2
len <- 15
space <- 4
for (i in 0:(len + 2 * space)) {
  p <- toPolar(tips_xx, tips_yy)
  p$r <- p$r + ri * i
  c <- toCart(p$r, p$th)
  if (i < space || i >= space + len) {
    points(c$x, c$y, col = "white", pch = 16, cex = 0.15)
  } else {
    trait_colors <- combo_colors[trait_combo_vector[tree$tip.label]]
    points(c$x, c$y, col = trait_colors, pch = 16, cex = 0.25)
  }
}

# ------------------------------------------------------------------------------
# Family Labels and MRCAs
# ------------------------------------------------------------------------------
family_list <- unique(traits$family)

family_nodes <- sapply(family_list, function(fam) {
  species <- traits %>% filter(family == fam) %>% pull(tips)
  tryCatch(findMRCA(tree, species, type = "node"), error = function(e) NA)
})

family_nodes <- na.omit(family_nodes)

r_ring <- max(sqrt(tips_xx^2 + tips_yy^2)) + ri * (space + len)

r_label <- r_ring + 0.02 * r_ring

for (i in seq_along(family_nodes)) {
  fam <- names(family_nodes)[i]
  fam_tips <- which(traits$family == fam)
  fam_tip_labels <- traits$tips[fam_tips]
  fam_tip_indices <- match(fam_tip_labels, tree$tip.label)
  x_vals <- obj$xx[fam_tip_indices]
  y_vals <- obj$yy[fam_tip_indices]
  thetas <- atan2(y_vals, x_vals)
  mean_theta <- atan2(mean(sin(thetas)), mean(cos(thetas)))
  x_out <- r_label * cos(mean_theta)
  y_out <- r_label * sin(mean_theta)
  hjust <- ifelse(cos(mean_theta) >= 0, 0, 1)
  text(x_out, y_out, labels = fam, cex = 0.6, srt = 0, adj = c(hjust, 0.5))
}

mrca_color <- "#9467bd"
for (node in family_nodes) {
  points(obj$xx[node], obj$yy[node], pch = 16, cex = 1, col = mrca_color)
}

# ------------------------------------------------------------------------------
# Legend
# ------------------------------------------------------------------------------
legend("topleft",
       legend = names(combo_colors),
       col = combo_colors,
       pch = 16,
       pt.cex = 0.9,
       title = expression(bold("Trait Combination")),
       bty = "n",
       cex = 0.6,
       x.intersp = 1.3,
       inset = c(-0.08, -0.06))

dev.off()

#length(tree$tip.label) # Number of species
# 
# # Which model was used (can also describe in text)
# dredge_sociality[[10]]$model
# 
# # Confirm model details
# str(dredge_sociality[[10]])
# length(family_nodes)  # How many MRCAs were successfully plotted


