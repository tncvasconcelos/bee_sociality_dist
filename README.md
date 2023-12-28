# bee_sociality_dist


Note: packages rgeos, rgdal, and maptools were unfortunately retired as of October 2023. Much of the code here depends on functions from these three packages and, although the functionalities of these packages have been migrated to other packages, I still have to make changes to the code so that it works in computers that dont already have rgeos, rgdal, and maptools installed. 

----
Description of folders: 
 
- **curated_data/** 

Data that has been processed by one of the pipelines below, including filtered and thinned occurrence points, pruned phylogenetic tree, simplified sociality categories and summary statistics for environmental variables.

- **original_data/** 

Data that comes straight from source and has not been curated yet. It also include original data produced by this work (e.g. sociality categories, exotic species to remove, etc)

- **environmental_layers/** 

All environmental layers that are used to calculate summary statistics for each species (proxy for environmental preferences for each species)

- **shapefile_bees/**

Resulting shapefiles of geographic ranges for each sampled bee species based on binarized versions of species distribution modeling. 

- **shapefile_angiosperms/**

Resulting shapefiles of geographic ranges for each sampled angiosperm species based on binarized versions of species distribution modeling. Note: right now it only includes about 14,500 species for which I had data from another project. I'm using those just to test some of the pipeline, but we will need to download data again to get a more accurate picture of angiosperm diversity in the Americas if this is our focal area (which will probably sum up to closer to 120,000 species)
 
- **plots/** 

A folder to organize plots for visual inspection of data and results. 



----
Scripts:

> 00_utility_functions.R

Functions that can be sourced to perform tasks in other scripts.

> 01_data_curation.R

Script to curate data. It simplifies sociality categories and prunes the bee supermatrix phylogeny of Henríquez-Piskulich et al. (2024) to include only species with data in our dataset. It will also create a tree plot highlighting sociality categories at the tips.

> 02_getting_climate_dist.R

This script will: (1) load occurrence points from Dorey, et al. (2023); (2) filter those to keep only species sampled in the bee phylogeny of Henríquez-Piskulich et al. (2024); (3) thin occurrence points to include only one occurrence per grid cell per species (to decrease bias that comes from overcollecting in some areas); (4) remove species with "suspicious" distributions (e.g. domesticated, invasive, etc) based on table created by Aline after visual inspection of points; (5) overlay occurrence points with environmental layers from folder environmental_layers/ to calculate summary statistics (e.g. mean, standard deviation) for each species; (6) save tables with summary statistics in folder curated_data/. Note: table with occurrence points is too large to be on github. It has to be downloaded from the original publication: see Dorey et al. (2023) supplementary information, "OutputData/05_cleaned_database.csv".

> 03_phylANOVA.R

Runs phylANOVAS between sociality categories and environmental variables based on curated data, using the pruned phylogeny. It will also plot simple boxplots showing the distribution of environmental variables by socliaty category and perform a postdhoc test to test significancy of differences between categories.  

> 04_modeling_bee_dist.R

This script will use thinned occurrence points to perfom species distribution modeling analyses for indiviual species using the maxent algorithm from R package dismo. The function "GetRanges" wraps several other functions from the 00_utility_functions.R to automatize the many steps required to perform SDM analyses, including: (1) it will download climatic layers from worldclim based on preferred resolution and crop those to limit background point to area where points are; (2) it will check for colinearity among the 19 climatic layers and keep only those that are not strongly correlated based on the preferred threshold; (3) it will create background points and randomly separate train and test sets before analyses; (4) it will report back AUC values and it will also binarize models to calculate likely range size based on preferred threshold. Finally, the function will also automatically save individual shape files in folder shapefile_bees/.

> 05_ploting_sp_rich.R

Loads binarized shapefiles from shapefile_bees/ and shapefile_angiosperms/ and overlays them to create species richness maps.


> 06_mapping_points.R

Map filtered and thinned occurrence points for visual inspection.

----
  
**Climate layers from:**
  
Karger, et al. (2017). Climatologies at high resolution for the earth’s land surface areas. Scientific data, 4(1), 1-20.  
  
Trabucco, A., & Zomer, R. (2018). Global aridity index and potential evapo- transpiration (ET0) climate database v2. CGIAR Consortium for Spatial Information (CGIAR-CSI). Published online, available from the CGIAR-CSI GeoPortal at https://cgiarcsi.community
  
  
**Filtered GBIF points from:**
  
Dorey, et al. (2023). A globally synthesised and flagged bee occurrence dataset and cleaning workflow. Scientific Data, 10(1), 747.


**Bee phylogeny from:**

Henríquez-Piskulich, et al. (2024). A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila). Molecular Phylogenetics and Evolution, 190, 107963.
