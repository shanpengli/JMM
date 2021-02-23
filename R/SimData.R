##' @title Data generation for simulation
##' @param sim number of simulated datasets to be generated.
##' @param n sample size for each dataset.
##' @param tL lower limit of time point for spline basis.
##' @param tU upper limit of time point for spline basis.
##' @param nbreak number of knots for spline basis.
##' @param p01 number of covariates in the first biomarker.
##' @param p02 number of covariates in the second biomarker.
##' @param q_eta number of survival covariates.
##' @param t_max maximum number of repeated measurements allowed for each subject.
##' @param k_max dimension of random effects.
##' @param distr distributional assumption of random effects.
##' @param m_age true mean value of age variable in the longitudinal sub-model.
##' @param std_age true standard deviation of age variable in the longitudinal sub-model.
##' @param sigma_inv sigma_inv.
##' @param tbtheta true btheta matrix.
##' @param tbeta0 true value of beta0 matrix.
##' @param tbeta1 true value of beta1 vector.
##' @param tsigma true value of sigma square vector.
##' @param tsigmad true value of variance of random effects.
##' @param ttheta true value of theta vector.
##' @param teta true value of survival fixed effects.
##' @param tgamma true value of latent association parameter.
##' @param lambda0 true baseline hazard. 
##' @export
##'
##'

SimData <- function(sim, n, tL, tU, nbreak = 8, p01 = 1, p02 = 1, q_eta = 1,
                    t_max = 19, k_max = 2, distr = "Gaussian", m_age = 0, std_age = 2, 
                    sigmau_inv, tbtheta, tbeta0, tbeta1, tsigma, tsigmad, ttheta, 
                    teta, tgamma, lambda0 = 0.08) {
  
  j_max = 2;
  
  if (distr == "Gaussian") {
    distr <- 1
  } else if (distr == "Gamma") {
    distr <- 0
  } else {
    stop("The available assumptions of random effects are either Gaussian or Gamma distribution. 
         Please choose one of the aformentioned distributions for simulation.")
  }
  
  ##pass data to Simdata
  sigmau_invnew=tempfile(pattern = "", fileext = ".txt")
  writenh(sigmau_inv,sigmau_invnew)
  tbthetanew=tempfile(pattern = "", fileext = ".txt")
  writenh(tbtheta,tbthetanew)
  tbeta0new=tempfile(pattern = "", fileext = ".txt")
  writenh(tbeta0,tbeta0new)
  tbeta1new=tempfile(pattern = "", fileext = ".txt")
  writenh(tbeta1,tbeta1new)
  tsigmanew=tempfile(pattern = "", fileext = ".txt")
  writenh(tsigma,tsigmanew)
  tsigmadnew=tempfile(pattern = "", fileext = ".txt")
  writenh(tsigmad,tsigmadnew)
  tthetanew=tempfile(pattern = "", fileext = ".txt")
  writenh(ttheta,tthetanew)
  tetanew=tempfile(pattern = "", fileext = ".txt")
  writenh(teta,tetanew)
  
  myresult <- Simdata_main(n, sim, nbreak, tL, tU, q_eta, 
                           j_max, p01, p02, t_max, distr, m_age,
                           std_age, k_max,
                           tbthetanew, sigmau_invnew,
                           tbeta0new, tbeta1new,
                           tsigmanew, tthetanew,
                           tsigmadnew, tetanew,
                           tgamma, lambda0)
  
  
  return(myresult)
  
  
  
}