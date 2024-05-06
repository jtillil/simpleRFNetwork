setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

datroot = paste0(
  # "./resclassif",
  "./data/ndclassif",
  "_nn", 100,
  "_ng", 1000,
  "_ns", 1000,
  "_ndm", 0,
  "_mdg", 0,
  "_pdg", 0.5,
  "_ab", 0,
  ".Rdata"
)
load(datroot)

num_threads = 60

implist = list()
for (i in 1:100) {
  print(i)
  rfdat = dat[[i]]$data[1:500,]
  
  rf = simpleRFNetwork(pheno ~ .,
                       data = rfdat,
                       num_trees=500,
                       num_threads=num_threads,
                       splitobject="module",
                       splitmethod="logridge1",
                       varselection="none",
                       mtry="root",
                       varclusters = dat[[i]]$modules,
                       seed = 1L
  )
  
  # run var importance
  implist[[i]] = rf$variableImportance(type = "permutation", num_threads = num_threads)
}

save(implist, file = "./results/logridge_GIM_dist.Rdata")