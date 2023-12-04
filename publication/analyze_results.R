# load borutares
setwd(getSrcDirectory(function(){})[1])
saveroot = paste0(
  "./results/resclassif",
  "_", "LDA",
  "_", "Gini",
  "_ni", 20,
  "_nn", 100,
  "_ng", 1000,
  "_ns", 1000,
  "_ndm", 1,
  "_mdg", 0,
  "_ab", 1,
  ".Rdata"
)
load(saveroot)

# show classification
for (res in borutares) {
  print(sum(res$causalmodules %in% which(res$classification == 1)))
  # print(sum(res$classification))
}

for (i in 1:20) {
  print(sum(borutares[[2*i]]))
}


