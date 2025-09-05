# ==============================================================================
# Data exploration
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: clear environment, set working directory, load data
# ------------------------------------------------------------------------------
rm(list = ls())

wd <- "/Users/lenarh/Desktop/bee_sociality_dist"
setwd(wd)
plots_wd <- file.path(wd, "plots")

library(ggplot2)
library(dplyr)
library(patchwork)

traits <- read.csv("curated_data/bees_traits.csv")


# ------------------------------------------------------------------------------
# Create combined trait column and clean labels
# ------------------------------------------------------------------------------
traits$trait_combo <- paste(traits$sociality_binary, traits$nest_binary, sep = "/")

traits$trait_combo <- recode(traits$trait_combo,
                             "solitary/ground" = "Solitary/Ground",
                             "solitary/aboveground" = "Solitary/Above-ground",
                             "social/ground" = "Social/Ground",
                             "social/aboveground" = "Social/Above-ground"
)


# ------------------------------------------------------------------------------
# Create summary table for total species per trait syndrome
# ------------------------------------------------------------------------------
trait_summary_df <- as.data.frame(table(traits$trait_combo))
colnames(trait_summary_df) <- c("Trait_Combination", "Species_Count")

# Set factor level order for consistent plotting
trait_levels <- c("Solitary/Ground", "Solitary/Above-ground", "Social/Ground", "Social/Above-ground")
trait_summary_df$Trait_Combination <- factor(trait_summary_df$Trait_Combination, levels = trait_levels)


# ------------------------------------------------------------------------------
# Plot 1: Barplot of species per trait syndrome
# ------------------------------------------------------------------------------
plot1 <- ggplot(trait_summary_df, aes(x = Trait_Combination, y = Species_Count, fill = Trait_Combination)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Species_Count), vjust = -0.5, size = 5) +  # Removed family = "Times"
  scale_fill_viridis_d(option = "cividis", name = "Trait Combination") +
  labs(x = "Trait Combination",
       y = "Species Count") +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),  # Inside top-right
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

quartz()
print(plot1)

ggsave(
  filename = file.path(plots_wd, "species_per_trait.png"),
  plot = plot1,
  width = 8,
  height = 7,
  dpi = 900
)


# ------------------------------------------------------------------------------
# Plot 2: Stacked barplot of species per trait syndrome by family (sorted)
# ------------------------------------------------------------------------------
# Recode trait_combo factor levels (if not already done)
traits$trait_combo <- factor(traits$trait_combo, levels = trait_levels)

# Count number of species per trait_combo within each family
family_combo_counts <- traits %>%
  group_by(family, trait_combo) %>%
  summarise(Species_Count = n(), .groups = "drop")

# Total number of species per family
family_species_counts <- family_combo_counts %>%
  group_by(family) %>%
  summarise(Total_Species = sum(Species_Count), .groups = "drop")

# Sort family factor levels by descending species count
ordered_families <- family_species_counts %>%
  arrange(desc(Total_Species)) %>%
  pull(family)

# Apply ordering
family_combo_counts$family <- factor(family_combo_counts$family, levels = ordered_families)
family_species_counts$family <- factor(family_species_counts$family, levels = ordered_families)

# Plot
plot2 <- ggplot(family_combo_counts, aes(x = family, y = Species_Count, fill = trait_combo)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(data = family_species_counts,
            aes(x = family, y = Total_Species, label = Total_Species),
            vjust = -0.5, size = 4, inherit.aes = FALSE) +
  scale_fill_viridis_d(option = "cividis", name = "Trait Combination") +
  labs(x = "Family", y = "Species Count") +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# Show plot
quartz()
print(plot2)

# Save plot
ggsave(
  filename = file.path(plots_wd, "family_trait_breakdown.png"),
  plot = plot2,
  width = 8,
  height = 5,
  dpi = 900
)

# ------------------------------------------------------------------------------
# Combined plot
# ------------------------------------------------------------------------------
plot1 <- plot1 +
  theme(legend.position = "none") + # remove legend from plot 1
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 25, face = "bold"))

plot2 <- plot2 +
  theme(axis.title.y = element_blank()) + # remove y-axis label from plot 2
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 25, face = "bold"))

combined_plot <- plot1 + plot2

quartz()
combined_plot 

ggsave(
  filename = file.path(plots_wd, "combined_species_trait_plot.png"),
  plot = combined_plot,
  width = 14,
  height = 9,
  dpi = 900
)

ggsave(
  filename = file.path(plots_wd, "combined_species_trait_plot.pdf"),
  plot = combined_plot,
  width = 14,
  height = 8,
  dpi = 900
)
