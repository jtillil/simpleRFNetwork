setwd(getSrcDirectory(function(){})[1])

library(SeqNet)
library(parallel)
library(pracma)
library(tictoc)
library(Rfast)
library(matrixcalc)

library(simpleRFNetwork)

# source("../R/utility.R")
# source("../R/Data.R")
# source("../R/splitAlgorithms.R")
# source("../R/TreeVarClusters.R")
# source("../R/TreeVarClustersClassification.R")
# source("../R/Forest.R")
# source("../R/ForestClassification.R")
# source("../R/simpleRFNetwork.R")
# 
# source("../R/genGeneNetworkData.R")
# source("../R/genGeneNetworkDataClassification.R")
# source("../R/genGeneNetworkDataRegression.R")

source("./boruta.R")