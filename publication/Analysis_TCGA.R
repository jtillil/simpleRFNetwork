library(curatedTCGAData)
library(TCGAutils)
library(MultiAssayExperiment)

## download multiassayexperiments
RNAseq_data = curatedTCGAData(diseaseCode = "BRCA", assays = "RNASeqGene", version = "2.1.0", dry.run=FALSE)
mRNA_data = curatedTCGAData(diseaseCode = "BRCA", assays = "mRNAArray", version = "2.1.0", dry.run=FALSE)

## extract data from multiassayexperiments


## run boruta
