# # ============================================================================
# # Circular Phylogeny: ASR Branch Coloring + Trait Ring + Family Labels
# # ============================================================================
# # Circular phylogeny with branches colored by ancestral state reconstruction (ASR), 
# # a trait ring at the tips, and labeled families with MRCA node dots.
# # ==============================================================================

# ------------------------------------------------------------------------------
# Setup
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
                             "social/ground"   = "Social/Ground",
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
# Family MRCA nodes
# ------------------------------------------------------------------------------
family_list <- unique(traits$family)
family_nodes <- sapply(family_list, function(fam) {
  species <- traits %>% filter(family == fam) %>% pull(tips)
  tryCatch(findMRCA(tree, species, type = "node"), error = function(e) NA)
})
family_nodes <- na.omit(family_nodes)

# ------------------------------------------------------------------------------
# Coordinate conversions
# ------------------------------------------------------------------------------
toCart <- function(r, th, deg = FALSE) {
  if (deg) th <- th * pi / 180
  list(x = r * cos(th), y = r * sin(th))
}
toPolar <- function(x, y) {
  list(r = sqrt(x^2 + y^2), th = atan2(y, x))
}

# # ----------------------------------------------------------------------------
# # Pre-plot: draw filled concentric gray rings for periods
# # ----------------------------------------------------------------------------
# tree_height <- max(nodeHeights(tree))
# periods <- data.frame(
#   name = c("Cretaceous", "Paleogene", "Neogene"),
#   start = c(145, 66, 23),
#   end   = c(66, 23, 0)
# )
# convert_to_radius <- function(t) max(sqrt(tips_xx^2 + tips_yy^2)) * ((tree_height - t) / tree_height)
# 
# # Pre-compute tip radius for circle scaling
# fake_plot <- plot(tree, type = "fan", plot = FALSE)
# tips_xx <- fake_plot$xx[1:length(tree$tip.label)]
# tips_yy <- fake_plot$yy[1:length(tree$tip.label)]
# ring_radii <- convert_to_radius(periods$start)
# 
# # Load plotrix and draw rings
# library(plotrix)
# ring_colors <- c("#F0F0F0", "#D9D9D9", "#BDBDBD")
# for (i in seq_along(ring_radii)) {
#   r_inner <- ifelse(i == length(ring_radii), 0, ring_radii[i + 1])
#   r_outer <- ring_radii[i]
#   draw.circle(0, 0, radius = r_outer, nv = 200, border = NA, col = ring_colors[i])
#   if (r_inner > 0) {
#     draw.circle(0, 0, radius = r_inner, nv = 200, border = "white", col = "white")
#   }
# }

# ------------------------------------------------------------------------------
# Plot
# # ----------------------------------------------------------------------------
png("plots/ASR_phylogeny.png", width = 1800, height = 1800, res = 300)
par(mar = c(2, 2, 2, 2), xpd = TRUE)

plot(tree, type = "fan", edge.color = branch_colors, edge.width = 0.5,
     show.tip.label = FALSE, no.margin = FALSE)

obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips_xx <- obj$xx[1:n_tips]
tips_yy <- obj$yy[1:n_tips]
p_original <- toPolar(tips_xx, tips_yy)

# Trait ring
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

# Family labels
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

# Family MRCAs
mrca_color <- "#9467bd"
for (node in family_nodes) {
  points(obj$xx[node], obj$yy[node], pch = 16, cex = 1, col = mrca_color)
}

# Legend
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

# ------------------------------------------------------------------------------
# Time axis
# ------------------------------------------------------------------------------
# theta_axis <- 0
# n_ticks <- 5
# time_values <- seq(tree_height, 0, length.out = n_ticks)
# radii <- max(sqrt(tips_xx^2 + tips_yy^2)) * ((tree_height - time_values) / tree_height)
# 
# line_coords <- toCart(radii, rep(theta_axis, length(radii)))
# segments(0, 0, line_coords$x[1], line_coords$y[1], col = "black", lwd = 1)
# 
# tick_length <- 0.02 * max(radii)
# for (i in seq_along(radii)) {
#   tick_r <- radii[i]
#   tick_coords <- toCart(tick_r, theta_axis)
#   tick_perp <- toCart(tick_length, theta_axis + pi/2)
#   segments(tick_coords$x - tick_perp$x,
#            tick_coords$y - tick_perp$y,
#            tick_coords$x + tick_perp$x,
#            tick_coords$y + tick_perp$y,
#            col = "black", lwd = 1)
#   text(tick_coords$x - 1.2 * tick_perp$x,
#        tick_coords$y - 1.2 * tick_perp$y,
#        labels = paste0(round(time_values[i], 1), " Ma"),
#        cex = 0.5, pos = 1)
# }

# # Add period labels above the axis
# for (i in 1:nrow(periods)) {
#   mid_time <- (periods$start[i] + periods$end[i]) / 2
#   mid_r <- convert_to_radius(mid_time)
#   label_coords <- toCart(mid_r, theta_axis)
#   text(label_coords$x, label_coords$y + 0.08 * max(radii),
#        labels = periods$name[i], cex = 0.5, font = 3, pos = 3)
# }

dev.off()



length(tree$tip.label)              # Number of species

# Which model was used (can also describe in text)
dredge_sociality[[10]]$model

# Confirm model details
str(dredge_sociality[[10]])
length(family_nodes)  # How many MRCAs were successfully plotted

