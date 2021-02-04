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
##' @export
##'

jmspline <- function(ydata, cdata, mdata, sigmau_inv, tbtheta, 
                     tL, tU, nbreak, p01, p02, j_max, k_max = 2,
                     quadpoint = 10, maxiter = 2500, do.trace = FALSE) {
  
  if (do.trace) {
    trace=1;
  } else {
    trace=0;
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
  
  ##pass data datas to jmspline
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
  
  myresult = jmspline_main(n, n_total, tL, tU, p01, p02, q_b, q_eta, j_max, 
                           t_max, nbreak, k_max, quadpoint, maxiter, trace, 
                           ydatanew, mdatanew, cdatanew, sigmau_invnew, 
                           tbthetanew, xs, ws)
  myresult$type="jmspline"
  myresult$quadpoint <- quadpoint
  myresult$n = n
  myresult$n1 = n_total
  class(myresult) <- "JMM"
  
  return (myresult)
  
  
  
}