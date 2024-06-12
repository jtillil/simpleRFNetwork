######################################################################################################################
#### run this script to reproduce the results of Tillil, Hu: Group-based random forest for disease module detection ##
######################################################################################################################

## setup
setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

########################
#### simulation study ##
########################

## generate 100 classfication datasets in 7 scenarios
source("./gen_classification_data.R")

## group Lasso
source("./overlapgglasso.R")
source("./group_lasso.R")

## group-based RF
# !!!! this will take months if not heavily parallelized !!!!
# !!!! only execute on a dedicated server !!!!
source("./boruta_LDA_Classification.R")
source("./boruta_logridge1_Classification.R")
source("./boruta_PCA_Classification.R")

## collect results
source("./analyze_group_lasso.R")
source("./analyze_results.R")

##########################
#### TCGA dataset study ##
##########################

## load data
source("./load_TCGA_data.R")

## group-based RF
source("./tcga_LDA.R")
source("./tcga_Ridge.R")
source("./tcga_PCA.R")
# source("./tcga_LDA_combined.R")
# source("./tcga_Ridge_combined.R")
# source("./tcga_PCA_combined.R")

## standard RF
source("./tcga_RF.R")

## collect results
source("./analyze_TCGA.R")

########################
#### generate figures ##
########################

source("./generate_figures_GIM.R")
source("./generate_figures_simulations.R")
source("./generate_figures_TCGA.R")

