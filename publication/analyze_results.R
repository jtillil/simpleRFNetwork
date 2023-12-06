# load borutares
setwd(getSrcDirectory(function(){})[1])
saveroot = paste0(
  "./serverresults/resclassif",
  "_", "SVM",
  "_", "permutation",
  "_ni", 20,
  "_nn", 100,
  "_ng", 1000,
  "_ns", 1000,
  "_ndm", 1,
  "_mdg", 0,
  "_ab", 0.5,
  ".Rdata"
)
load(saveroot)
datroot = paste0("./data/ndclassif_", substr(saveroot, gregexpr("nn", saveroot)[[1]][1], nchar(saveroot)))
load(datroot)

# summarize classification
grandsum = 0
for (res in borutares) {
  # grandsum = grandsum + sum(res$causalmodules %in% which(res$classification == 1))
  # print(sum(res$causalmodules %in% which(res$classification == 1)))
  print(res$binomresults[res$causalmodules])
  # print(sum(res$classification))
}
# print(grandsum)

# analyze overlap between selected modules and causal modules








