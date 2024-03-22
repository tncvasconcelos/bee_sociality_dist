library(corHMM)
library(ggplot2)
# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

# Some exploratory analyses
# Reloading traits, tree and climatic data:
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv")

traits <- traits[,c(3:5)]

# traits <- subset(traits, grepl("Xylocopa",traits$tips))
# tree <- keep.tip(tree, which(tree$tip.label%in%traits$tips))

# Correlation between nesting habit and sociality
corhmm_fits <- corHMM:::fitCorrelationTest(tree, traits) 
save(corhmm_fits, file = "results/corhmm_fits.Rsave")
load("results/corhmm_fits.Rsave")
corhmm_tbl <- corHMM:::getModelTable(corhmm_fits)

# conducting a lrt
teststat <- -2 * (corhmm_tbl$lnLik[1] - corhmm_tbl$lnLik[2])
p.val50 <- pchisq(teststat, df = 1, lower.tail = FALSE)
print(p.val50)

#---------------------
# plot
#pal <- hcl.colors(5, palette = "Viridis", alpha = 0.7)

subset_traits <- subset(traits, traits$sociality!="parasite")
subset_traits <- subset(subset_traits, subset_traits$nest!="parasite")

subset_traits$combined <- paste(subset_traits$sociality,subset_traits$nest,sep="_")

table(subset_traits$combined)

t0 <- ggplot(subset_traits, aes(x=sociality, fill=nest)) + 
  geom_bar(position = "fill", alpha=0.8)+
  theme_bw(base_size = 8) +
  theme(legend.position = "none") + 
  scale_fill_manual(values = c("#66a182","#edae49")) +
  #annotate(geom="text", x=0.5, y=max(corolla_diam), size=4, hjust=0, label=paste0("p=",phylanova_results$Pf)) + 
  xlab("sociality") +
  ylab("Proportion in ground") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 
