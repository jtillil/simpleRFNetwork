# load borutares
setwd(getSrcDirectory(function(){})[1])
# saveroot = paste0(
#   # "./resclassif",
#   "./serverresults/serverres_24_03_09/resclassif",
#   "_", "LDA",
#   "_", "permutation",
#   "_ni", 20,
#   "_nn", 100,
#   "_ng", 1000,
#   "_ns", 1000,
#   "_ndm", 2,
#   "_mdg", 0,
#   "_pdg", 0.5,
#   "_ab", 0.5,
#   ".Rdata"
# )
# load(saveroot)
# print("load res done")
# datroot = paste0("./data/ndclassif_", substr(saveroot, gregexpr("nn", saveroot)[[1]][1], nchar(saveroot)))
# load(datroot)
# print("load dat done")

# summarize classification
# num_classified_modules = c()
# num_selected = c()
# found_causal_modules = c()

# for (res in borutares) {
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
  # print(c(length(res$first_classification), length(which(res$first_classification %in% c(1))), length(which(res$second_classification %in% c(1)))))
  # found_causal_modules = c(found_causal_modules, sum(which(res$aggregated_classifications %in% c(1)) %in% res$causalmodules))
  
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
  
# }
# print(colSums(num_selected, TRUE))

# hist(num_classified_modules, 18)

# analyze overlap between selected modules and causal modules
# print(sum(found_causal_modules))
# hist(found_causal_modules)

#### Null case ####
# library(latex2exp)
# library(ggplot2)
# library(dplyr)
# 
# # modulecountdat <- data.frame(
# #   Method = c(),
# #   Size = c(),
# #   Count = c()
# # )
# 
# GIMdistributiondat <- data.frame(
#   Method = character(),
#   Size = double(),
#   GIM = double()
# )
# 
# for (method in c("LDA", "logridge1", "PCA")) {
#   saveroot = paste0(
#     # "./resclassif",
#     "./serverresults/serverres_24_04_18/resclassif",
#     "_", method,
#     "_", "permutation",
#     "_ni", 20,
#     "_nn", 100,
#     "_ng", 1000,
#     "_ns", 1000,
#     "_ndm", 0,
#     "_mdg", 0,
#     "_pdg", 0.5,
#     "_ab", 0,
#     ".Rdata"
#   )
#   datroot = paste0(
#     # "./resclassif",
#     "./data/ndclassif",
#     "_nn", 100,
#     "_ng", 1000,
#     "_ns", 1000,
#     "_ndm", 0,
#     "_mdg", 0,
#     "_pdg", 0.5,
#     "_ab", 0,
#     ".Rdata"
#   )
#   if (method != "logridge1") {
#     if (file.exists(saveroot)) {
#       print(c(method))
#       load(saveroot)
#       load(datroot)
#       
#       for (i in 1:length(borutares)) {
#         print(i)
#         res = borutares[[i]]
#         if (!is.null(res)) {
#           network = dat[[i]]
#           for (j in 1:length(network$modules)) {
#             for (k in 1:1) {
#               GIMdistributiondat[nrow(GIMdistributiondat) + 1,] = c(method, as.numeric(lengths(network$modules)[j]), as.numeric(res$first_vim[k, j]))
#             }
#           }
#         }
#       }
#     }
#   } else {
#     for (i in 1:length(borutares)) {
#       print(i)
#       imp = implist[[i]]
#       if (!is.null(res)) {
#         network = dat[[i]]
#         for (j in 1:length(network$modules)) {
#           for (k in 1:1) {
#             GIMdistributiondat[nrow(GIMdistributiondat) + 1,] = c(method, as.numeric(lengths(network$modules)[j]), as.numeric(imp[j]))
#           }
#         }
#       }
#     }
#   }
# }
# GIMdistributiondat$Method[GIMdistributiondat$Method == "logridge1"] = "Ridge"
# GIMdistributiondat$Method = factor(GIMdistributiondat$Method, levels = c("LDA", "Ridge", "PCA"))
# GIMdistributiondat$Size = as.numeric(GIMdistributiondat$Size)
# GIMdistributiondat$GIM = as.numeric(GIMdistributiondat$GIM)
# 
# # LDAdat = GIMdistributiondat[GIMdistributiondat$Method == "LDA", ]
# 
# p = ggplot(GIMdistributiondat, aes(Size, GIM)) +
#   geom_point() +
#   scale_x_continuous(breaks = 25*(0:4)) +
#   theme_bw() +
#   xlab("Module size") +
#   ylab("Group permutation importance") +
#   facet_grid(rows = vars(Method), scales = "free")
# plot(p)
# 
# ggsave("scatter_GIM_Size_Null.pdf", width = 7, height = 5)

#### individual result ####

# res = borutares[[100]]

#### ggplot ####
library(latex2exp)
library(ggplot2)
library(dplyr)

detectiondat <- data.frame(
  Method = rep(rep(c("Lasso", "LDA", "logridge1", "PCA"), each = 3), 6),
  Number = rep(rep(c("0", "1", "2"), times = 4), 6),
  Count = rep(0, 72),
  ndm = c(rep(1, 36), rep(2, 36)),
  ndm_plot = c(rep("1 disease module per network", 36), rep("2 disease modules per network", 36)),
  ab = rep(c(rep(0.5, 12), rep(1, 12), rep(2, 12)), times = 2),
  ab_plot = rep(c(rep("beta = 0.5", 12), rep("beta = 1", 12), rep("beta = 2", 12)), times = 2)
)

for (method in c("LDA", "PCA", "logridge1")) {
  for (ndm in c(1, 2)) {
    for (ab in c(0.5, 1, 2)) {
      # setwd(getSrcDirectory(function(){})[1])
      saveroot = paste0(
        # "./resclassif",
        "./serverresults/serverres_24_04_18/resclassif",
        "_", method,
        "_", "permutation",
        "_ni", 20,
        "_nn", 100,
        "_ng", 1000,
        "_ns", 1000,
        "_ndm", ndm,
        "_mdg", 0,
        "_pdg", 0.5,
        "_ab", ab,
        ".Rdata"
      )
      if (file.exists(saveroot)) {
        print(c(method, ndm, ab))
        load(saveroot)

        for (res in borutares) {
          found_causal_modules = sum(which(res$aggregated_classifications %in% c(1)) %in% res$causalmodules)
          idx = (detectiondat$Method == method & detectiondat$ndm == ndm & detectiondat$ab == ab & detectiondat$Number == found_causal_modules)
          detectiondat$Count[idx] = detectiondat$Count[idx] + 1
        }
      }
    }
  }
}

for (ndm in c(1, 2)) {
  for (ab in c(0.5, 1, 2)) {
    # setwd(getSrcDirectory(function(){})[1])
    saveroot = paste0(
      # "./resclassif",
      "./results/grplasso",
      "_nn", 100,
      "_ng", 1000,
      "_ns", 1000,
      "_ndm", ndm,
      "_mdg", 0,
      "_pdg", 0.5,
      "_ab", ab,
      ".Rdata"
    )
    if (file.exists(saveroot)) {
      print(c("Lasso", ndm, ab))
      load(saveroot)
      
      for (res in grplassores) {
        found_causal_modules = res$number_correctly_selected
        idx = (detectiondat$Method == "Lasso" & detectiondat$ndm == ndm & detectiondat$ab == ab & detectiondat$Number == found_causal_modules)
        detectiondat$Count[idx] = detectiondat$Count[idx] + 1
      }
    }
  }
}

# detectiondat <- detectiondat %>%
#   group_by(Method, Number) %>%
#   summarise(TotalCount = sum(Count), .groups = "drop_last")

detectiondat = detectiondat[!(detectiondat$ndm == 1 & detectiondat$Number == 2),]

detectiondat <- detectiondat %>%
  group_by(Method, ndm, ab) %>%
  mutate(PercentWithinGroup = Count / sum(Count) * 100)

detectiondat$ab[detectiondat$ab == 0.5] = "\beta = 0.5"
detectiondat$ab[detectiondat$ab == 1] = "\beta = 1"
detectiondat$ab[detectiondat$ab == 2] = "\beta = 2"

ggplot(detectiondat, aes(x = Number, y = PercentWithinGroup, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  labs(
    # title = "Number of detected disease modules per network",
       x = "Number of disease modules detected",
       y = "Percentage of replications") +
  theme_bw() +
  scale_fill_discrete(name = "Method", labels = c("Group Lasso", "Group RF: LDA", "Group RF: Ridge", "Group RF: PCA")) +
  # legend(c("LDA", "Ridge", "PCA")) +
  facet_grid(ndm_plot ~ ab_plot, scales = "free")#, labeller = label_parsed)
             # labeller = label_parsed(ab = c("beta=0.5", "beta=1", "beta=2"), ndm = c("disease modules = 1", "disease modules = 2")))

ggsave("bar_number_detected_modules.pdf", width = 7, height = 5)
