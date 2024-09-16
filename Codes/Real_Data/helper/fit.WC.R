#####################
# Wilcoxon Wrappper #
#####################

fit.Wilcoxon  = function(features, metadata, libSize, ID, MultTestCorrection = "fdr") {

  ###################
  # Subset by group #
  ###################

  group0indexes <- which(metadata == "astrocytes")
  group1indexes <- which(metadata == "oligodendrocytes")

  #################
  # Main function #
  #################

  spe <- function(x){
    tryCatch({
      fit1 <- wilcox.test(x[group0indexes], x[group1indexes])
    }, error=function(err){
      fit1 <- try({wilcox.test(x[group0indexes], x[group1indexes])})
      return(fit1)
    })}

  ###################
  # Collect results #
  ###################

  spes <- apply(features, 2, spe)
  output_df <- data.frame(pval = sapply(spes, function(x) x$p.value))
  output_df$coef <- sapply(spes, function(x) x$statistic)
  output_df$metadata<-colnames(metadata)
  output_df$feature<-colnames(features)
  output_df$qval<-as.numeric(p.adjust(output_df$pval, method = MultTestCorrection))
  output_df<-output_df[order(output_df$qval, decreasing=FALSE),]
  output_df<-dplyr::select(output_df, c('feature', 'metadata'), everything())
  rownames(output_df)<-NULL
  return(output_df)
}
