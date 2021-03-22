##' A Flexible Joint Model for Multiple Longitudinal Biomarkers and A Time-to-Event Outcome
##' @title A Flexible Joint Model for Multiple Longitudinal Biomarkers and A Time-to-Event Outcome
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
##' @param conversigmad logical; TRUE if sigmad is considered into convergence criteria. Default is FALSE.
#' @return Object of class \code{JMM} with elements
##'   \tabular{ll}{
##'       \code{ydata}    \tab  a data frame of inputted longitudinal data with multiple biomarkers. \cr
##'       \code{cdata} \tab a data frame of inputted survival data. \cr
##'       \code{N1} \tab a vector containing the number of repeated measurements for all subjects. \cr
##'       \code{beta0_matrix} \tab  a matrix of point estimates of \eqn{\beta_0}. \cr
##'       \code{beta1_estimate}   \tab  a vector of point estimates of \eqn{\beta_1}. \cr
##'       \code{sigma2_estimate} \tab a vector of point estimates of \eqn{\sigma^2} for each biomarker. \cr
##'       \code{theta_estimate}    \tab a vector of point estimates of \eqn{\theta}. \cr
##'       \code{sigmad_estimate}     \tab  a vector of point estimates of \eqn{\sigma_d^2}, i.e., the diagonal elements of the variance-covariance matrix of random effects.\cr
##'       \code{eta_estimate}     \tab  a vector of The point estimates of survival fixed effects \eqn{\eta}.\cr
##'       \code{btheta_matrix}     \tab a matrix of point estimates of \eqn{\Theta}.\cr
##'       \code{gamma}     \tab the association parameter \eqn{\gamma}. \cr
##'       \code{BaselineHazard}     \tab a matrix of baseline hazard. The first row denotes a distinct survival time, 
##'       the second row denotes the number of events, and the third row denotes the hazard rate on the corresponding survival time.\cr
##'       \code{FUNA_matrix}     \tab a matrix of expected value of random effects \eqn{\alpha}. \cr
##'       \code{iter} \tab 	the number of iterations in the EM algorithm. \cr
##'       \code{loglike} \tab 	the log-likelihood value. \cr
##'       \code{AIC} \tab 	the AIC value. \cr
##'       \code{TotalPara} \tab 	the total number of parameters of the model fit. \cr
##'       \code{SurvivalVar} \tab 	an identifier of whether survival covariates are considered in the joint model. 
##'       1 corresponds to yes whereas 0 to no. \cr
##'       \code{type} \tab 	object identifier. \cr
##'       \code{quadpoint} \tab 	the number of standard Gaussian Hermite quadrature points used for the E step. \cr
##'       \code{n} \tab 	total number of subjects. \cr
##'       \code{n1} \tab 	total number of observations. \cr
##'       \code{sigma_inv} \tab 	a matrix of \eqn{\Sigma^{-1}}. \cr
##'       \code{tbtheta} \tab 	a matrix of inputted \eqn{\Theta}. \cr
##'       \code{tL} \tab 	lower limit of time for spline basis. \cr
##'       \code{tU} \tab 	upper limit of time for spline basis. \cr
##'       \code{p01} \tab the number of longitudinal fixed effects for the first biomarker. \cr
##'       \code{p02} \tab the number of longitudinal fixed effects for the second biomarker. \cr
##'       \code{nbreak} \tab 	the number of knots used for spline basis. \cr
##'       \code{j_max} \tab 	total number of biomarkers in the longitudinal data. \cr
##'       \code{k_max} \tab 	total number of random effects. \cr
##'       
##'   }
##' @examples
##'
##' require(JMM)
##' ## load supporting data for data generation and model fit 
##' data(simdata)
##' tbeta0 <- as.data.frame(matrix(c(-0.4, -0.3), nrow = 2))
##' tbeta1 <- as.data.frame(c(1, 0.8))
##' tsigma <- as.data.frame(c(1, 1))
##' tsigmad <- as.data.frame(c(2, 1))
##' ttheta <- as.data.frame(c(5, 0.148183, -1.606910, 1.083067, -1.788186, -1.205329,
##'                           0.245283, -1.433481, -0.603650, -0.134991))
##' teta <- as.data.frame(c(0))
##' tgamma <- 0.5
##' tL <- 0
##' tU <- 9
##' ## generate a simulated data
##' Data <- SimData(sim = 1, n = 215, tL = tL, tU = tU, nbreak = 8, p01 = 1, p02 = 1,
##' q_eta = 1, t_max = 19, k_max = 2, distr = "Gaussian", m_age = 0,
##' std_age = 2, sigmau_inv = sigmau_inv, tbtheta = tbtheta, 
##' tbeta0 = tbeta0, tbeta1 = tbeta1, tsigma = tsigma, tsigmad = tsigmad,
##' ttheta = ttheta, teta = teta, tgamma = tgamma)
##' ydata <- as.data.frame(Data$ydata_0)
##' cdata <- as.data.frame(Data$cdata_0)
##' mdata <- as.data.frame(Data$mdata_0)
##' \donttest{
##' ## fit a joint model
##' fit <- jmspline(ydata = ydata, cdata = cdata, mdata = mdata, sigmau_inv = sigmau_inv, 
##'                 tbtheta = btheta, tL = tL, tU = tU, nbreak = 8, p01 = 1, p02 = 1, 
##'                 j_max = 2, beta0init = tbeta0, beta1init = tbeta1, sigmainit = tsigma,
##'                 thetainit = ttheta, sigmadinit = tsigmad, gammainit = tgamma,
##'                 survVar = FALSE, conversigmad = TRUE)
##' }
##' @references
##' \itemize{
##' \item Ning Li, Yi Liu, Shanpeng Li, Robert M. Elashoff, and Gang Li. "A Flexible Joint Model for Multiple Longitudinal Biomarkers and A Time-to-Event Outcome: 
##' With Applications to Dynamic Prediction Using Highly Correlated Biomarkers." Biometrical Journal (Accepted).
##' }
##' @seealso \code{\link{SimData}}
##' @export
##'

jmspline <- function(ydata, cdata, mdata, sigmau_inv, tbtheta, 
                     tL, tU, nbreak, p01, p02, j_max, k_max = 2,
                     quadpoint = 10, maxiter = 2500, do.trace = FALSE,
                     beta0init = NULL, beta1init = NULL, sigmainit = NULL,
                     thetainit = NULL, sigmadinit = NULL, gammainit = NULL, survVar = TRUE,
                     conversigmad = FALSE) {
  
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
  
  if (conversigmad) {
    conversigmad=1;
  } else {
    conversigmad=0;
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
                           survvar, conversigmad)
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