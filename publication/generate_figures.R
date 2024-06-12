########################
#### generate figures ##
########################

## setup
setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

## generate figures
source("./generate_figures_GIM.R")
source("./generate_figures_simulations.R")
source("./generate_figures_TCGA.R")