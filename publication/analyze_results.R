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
num_classified_modules = c()
num_selected = c()

for (res in borutares) {
  # grandsum = grandsum + sum(res$causalmodules %in% which(res$classification == 1))
  # print(sum(res$causalmodules %in% which(res$classification == 1)))
  # print(res$second_binomresults[res$causalmodules])
  # print(res$second_classification[res$causalmodules])
  # num_classified_modules = c(num_classified_modules, match(1, res$second_classification))
  # num_classified_modules = c(num_classified_modules, which(res$second_classification %in% c(1)))
  # num_selected = rbind(num_selected, res$second_classification[res$causalmodules])
  # print(which(res$second_classification %in% c(1)))
  # print(res$causalmodules)
  # print(res$second_classification[res$causalmodules])
  # print(sum(res$classification))
  # print(c(length(res$first_classification), length(which(res$first_classification %in% c(-1))), length(which(res$second_classification %in% c(1)))))
  print(c(length(res$first_classification), length(which(res$first_classification %in% c(0, 1))), length(which(res$second_classification %in% c(0, 1)))))
  
  
  ## classify modules
  # alt_classification = rep(0, length(res$first_binomresults))
  # alt_classification[(res$first_classification != -1)] = NA
  # upcutoff = qbinom(0.99999, 20, 0.5)
  # 
  # for (i in 1:length(res$updated_modules) ) {
  #   if (res$second_binomresults[(res$first_classification != -1)][i] > upcutoff) {
  #     alt_classification[(res$first_classification != -1)][i] = 1
  #   }
  # }
  # print(which(alt_classification %in% c(1)))
  
}
print(colSums(num_selected, TRUE))

hist(num_classified_modules, 18)

# analyze overlap between selected modules and causal modules


#### individual result ####

res = borutares[[100]]



