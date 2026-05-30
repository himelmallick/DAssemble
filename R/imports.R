#' @importFrom dplyr everything
#' @importFrom stats as.formula binomial coef glm model.matrix na.exclude na.pass p.adjust pcauchy plogis sd update
#' @importFrom utils capture.output combn write.table
NULL

utils::globalVariables(c(
  ".",
  "bin_cat",
  "Pr(>Chisq)",
  "ci.hi",
  "ci.lo",
  "component",
  "contrast",
  "data",
  "name",
  "primerid",
  "tweedie.index",
  "xnames"
))
