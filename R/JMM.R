#' @useDynLib JMM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats  as.formula  pnorm  pchisq complete.cases
#' @importFrom statmod  gauss.quad
#' @importFrom utils  read.table
#' @importFrom parallel mclapply parLapply makeCluster stopCluster
NULL