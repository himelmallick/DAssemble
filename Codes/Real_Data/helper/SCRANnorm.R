#### Apply SCRAN Normalization To A Dataset ####

SCRANnorm = function(features) {

  #### Extract Features ####

  features <- as.matrix(features)

  #### SCRAN Normalizing the Data ####

  sce <- t(features)
  sizeFactors <- scran::calculateSumFactors(sce)

  #### Reset library size ####

  libSize <- sizeFactors
  return(libSize)
}
