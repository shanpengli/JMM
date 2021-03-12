##' @title Get fitted values for multiple markers over time
##' @param object a object inheriting from class JMM.
##' @param long.covar a matrix of longitudinal covariates for biomarkers. 
##' Each row of observations stand for a particuler biomarker. 
##' If NULL the mean values of covariates in the longitudianl data will be adopted. Default is NULL.
##' @export
##'

fittedval <- function(object, long.covar = NULL) {
  if (!inherits(object, "JMM"))
    stop("Use only with 'JMM' objects.\n")
  
  ydata <- object$ydata
  
  if (is.null(long.covar)) {
    longvar <- matrix(0, nrow = object$j_max, ncol = max(object$p01, object$p02))
    long.covar <- colMeans(ydata[, c(3:(2+object$p01), (4+object$p01):(3+object$p01+object$p02))])
    for (i in 1:object$p01) longvar[1, i] <- long.covar[i]
    for (i in 1:object$p02) longvar[2, i] <- long.covar[object$p01 + i]
    longvar <- as.data.frame(longvar)
  } else if (is.data.frame(long.covar) && ncol(long.covar) == max(object$p01, object$p02) &&
             nrow(long.covar) == object$j_max) {
    longvar <- long.covar
  } else {
    stop("The input matrix are incorrect. Make sure if the dimension mataches the object.")
  }
  longvarnew=tempfile(pattern = "", fileext = ".txt")
  writenh(longvar,longvarnew)
  
  tL <- object$tL
  tU <- object$tU
  nbreak <- object$nbreak
  k_max <- object$k_max
  j_max <- object$j_max
  p01 <- object$p01
  p02 <- object$p02
  sigmau_inv <- object$sigmau_inv
  sigmau_invnew=tempfile(pattern = "", fileext = ".txt")
  writenh(sigmau_inv,sigmau_invnew)
  
  theta <- as.data.frame(object$theta_estimate)
  thetanew=tempfile(pattern = "", fileext = ".txt")
  writenh(theta,thetanew)
  
  btheta <- as.data.frame(object$btheta_matrix)
  bthetanew=tempfile(pattern = "", fileext = ".txt")
  writenh(btheta,bthetanew)
  
  beta0 <- as.data.frame(object$beta0_matrix)
  beta0new=tempfile(pattern = "", fileext = ".txt")
  writenh(beta0,beta0new)
  beta1 <- as.data.frame(object$beta1_estimate)
  beta1new=tempfile(pattern = "", fileext = ".txt")
  writenh(beta1,beta1new)
  
  myresult <- getfitted_main(tL, tU, nbreak, k_max, j_max, p01, p02, 
                             sigmau_invnew, thetanew, bthetanew, beta0new, beta1new, longvarnew)
  
  return(myresult)
  
  
}