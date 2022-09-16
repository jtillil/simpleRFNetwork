% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genGeneNetworkData.R
\name{genGeneNetworkData}
\alias{genGeneNetworkData}
\title{simpleRFNetwork}
\usage{
genGeneNetworkData(
  num_networks,
  num_genes,
  num_modules,
  num_observations,
  num_causal_modules,
  num_causal_genes,
  effect_measure
)
}
\arguments{
\item{num_networks}{Integer, number of networks to generate.}

\item{num_genes}{Integer, number of genes per network.}

\item{num_modules}{Integer, number of modules per network.}

\item{num_observations}{Integer, number of expression data observations per network.}

\item{num_causal_modules}{Integer, how many modules should be causal for the phenotype.}

\item{num_causal_genes}{Integer, how many genes in each causal module should be causal for the phenotype. Can also be "all".}

\item{effect_measure}{Float, standardized effect measure of causal genes in the network.}
}
\description{
Generates random networks of genes, their associations, expression data and
binary labels where the effect measures of individual genes can be set by
the user.
}
\examples{
\donttest{
library(simpleRFNetwork)
library(SeqNet) 


}

}
\references{
Breiman, L. (2001). Random forests. Mach Learn, 45(1), 5-32. \cr
}
\author{
Johannes Tillil
}