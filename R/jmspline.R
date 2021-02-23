##' Joint modeling of longitudinal continuous and survival data
##' @title Joint Modelling for Continuous outcomes
##' @param ydata a multibiomarker longitudinal data frame in long format.
##' @param cdata a survival data frame with single failure.
##' Each subject has one data entry.
##' @param mdata a data vector of number of repeated measurements.
##' @param sigmau_inv sigmau_inv.
##' @param tbtheta a data frame of the normalizing constant for spline basis.
##' @param tL the lower limit of time interval for spline basis.
##' @param tU the upper limit of time interval for spline basis.
##' @param nbreak number of knots to be specified.
##' @param p01 number of covariates in the first biomarker.
##' @param p02 number of covariates in the second biomarker.
##' @param j_max Total number of biomarkers in the longitudinal data.
##' @param k_max Number of random effects to be specified. Default is 2.
##' @param quadpoint Number of quadrature points.
##' @param maxiter The maximum number of iterations for EM process.
##' @param do.trace A logic value to control the iterated results to be printed. Default is FALSE.
##' @param beta0init Initial values of beta0. 
##' @param beta1init Initial values of beta1. 
##' @param sigmainit Initial values of sigma^2. 
##' @param thetainit Initial values of theta. 
##' @param sigmadinit Initial values of sigmad. 
##' @param gammainit Initial values of gamma.
##' @param survVar logical; TRUE if the survival sub-model include covariates. Default is TRUE.  
##' @export
##'

jmspline <- function(ydata, cdata, mdata, sigmau_inv, tbtheta, 
                     tL, tU, nbreak, p01, p02, j_max, k_max = 2,
                     quadpoint = 10, maxiter = 2500, do.trace = FALSE,
                     beta0init = NULL, beta1init = NULL, sigmainit = NULL,
                     thetainit = NULL, sigmadinit = NULL, gammainit = NULL, survVar = TRUE) {
  
  if (do.trace) {
    trace=1;
  } else {
    trace=0;
  }
  
  if (survVar) {
    survvar=1;
  } else {
    survvar=0;
  }
  
  #Gaussian-Hermite quadrature nodes and weights
  #The dimension of xs/ws is half of the quadpoint value since they are symmetric
  
  if (quadpoint %% 2 == 1)
  {
    stop("Number of quadrature quadpoints can only be even!")
  }
  
  gq_vals <- statmod::gauss.quad(n = quadpoint, kind = "hermite")
  
  xs <- gq_vals$nodes[(quadpoint / 2 + 1) : quadpoint]
  
  ws <- gq_vals$weights[(quadpoint / 2 + 1) : quadpoint]
  
  ydim = dim(ydata)
  cdim = dim(cdata)
  n1 = ydim[1]
  
  ##q_b: dimension of cubic B-spline basis functions
  q_b = 4 + nbreak - 2
  
  ## sample size
  n = cdim[1]
  
  ## number of repeated measurements
  n_total = ydim[1]
  
  ## max # of observations per subject
  t_max = max(mdata)
  
  ## number of survival covariates
  q_eta = cdim[2] - 2
  
  ##pass the initials into C++
  if (nrow(beta0init) != j_max || ncol(beta0init) != max(p01, p02)) {
    stop("The dimension of beta0 initials is incorrect!")
  }
  
  if (nrow(beta1init) != j_max) {
    stop("The dimension of beta1 initials is incorrect!")
  }
  
  if (nrow(sigmainit) != j_max) {
    stop("The dimension of sigma^2 initials is incorrect!")
  }
  
  if (nrow(thetainit) != q_b) {
    stop("The dimension of theta initials is incorrect!")
  }
  
  if (nrow(sigmadinit) != k_max) {
    stop("The dimension of sigmad initials is incorrect!")
  }
  
  if (!is.numeric(gammainit)) {
    stop("Gamma initials can only be numeric.")
  }
  
  
  ##pass data to jmspline
  ydatanew=tempfile(pattern = "", fileext = ".txt")
  writenh(ydata,ydatanew)
  mdatanew=tempfile(pattern = "", fileext = ".txt")
  writenh(mdata,mdatanew)
  cdatanew=tempfile(pattern = "", fileext = ".txt")
  writenh(cdata,cdatanew)
  sigmau_invnew=tempfile(pattern = "", fileext = ".txt")
  writenh(sigmau_inv,sigmau_invnew)
  tbthetanew=tempfile(pattern = "", fileext = ".txt")
  writenh(tbtheta,tbthetanew)
  
  
  #pass initals to jmspline
  beta0initnew=tempfile(pattern = "", fileext = ".txt")
  writenh(beta0init,beta0initnew)
  beta1initnew=tempfile(pattern = "", fileext = ".txt")
  writenh(beta1init,beta1initnew)
  sigmainitnew=tempfile(pattern = "", fileext = ".txt")
  writenh(sigmainit,sigmainitnew)
  thetainitnew=tempfile(pattern = "", fileext = ".txt")
  writenh(thetainit,thetainitnew)
  sigmadinitnew=tempfile(pattern = "", fileext = ".txt")
  writenh(sigmadinit,sigmadinitnew)
  
  myresult = jmspline_main(n, n_total, tL, tU, p01, p02, q_b, q_eta, j_max, 
                           t_max, nbreak, k_max, quadpoint, maxiter, trace, 
                           ydatanew, mdatanew, cdatanew, sigmau_invnew, 
                           tbthetanew, xs, ws, beta0initnew, beta1initnew,
                           sigmainitnew, thetainitnew, sigmadinitnew, gammainit,
                           survvar)
  myresult$type="jmspline"
  myresult$quadpoint <- quadpoint
  myresult$n = n
  myresult$n1 = n_total
  myresult$sigmau_inv = sigmau_inv
  myresult$tbtheta = tbtheta
  myresult$tL = tL
  myresult$tU = tU
  myresult$p01 = p01
  myresult$p02 = p02
  myresult$nbreak = nbreak
  myresult$j_max = j_max
  myresult$k_max = k_max
  
  
  
  class(myresult) <- "JMM"
  
  return (myresult)
  
  
  
}