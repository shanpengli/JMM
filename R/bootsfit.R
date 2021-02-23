##' @export
##'
bootsfit <- function(i, Data = Data, nycol = nycol, nccol = nccol, 
                     sigmau_inv = sigmau_inv, tbtheta = tbtheta, tL = tL, tU = tU,
                     nbreak = nbreak, p01 = p01, p02 = p02, j_max = j_max,
                     k_max = k_max, quadpoint = quadpoint, maxiter = maxiter,
                     beta0init = tbeta0, beta1init = tbeta1,
                     sigmainit = tsigma2, thetainit = ttheta, 
                     sigmadinit = tsigmad,
                     gammainit = tgamma,
                     survVar = survVar) {
  coef <- vector()
  bootsydata <- matrix(unlist(Data[1+(i-1)*3]), ncol = nycol)
  bootsydata <- as.data.frame(bootsydata)
  bootscdata <- matrix(unlist(Data[2+(i-1)*3]), ncol = nccol)
  bootscdata <- as.data.frame(bootscdata)
  bootsmdata <- matrix(unlist(Data[3*i]), ncol = 1)
  bootsmdata <- as.data.frame(bootsmdata)
  
  p_max <- max(p01, p02)
  q_b <- 4 + nbreak - 2
  cdim <- dim(bootscdata)
  
  q_eta = cdim[2] - 2
  
  fit <- jmspline(ydata = bootsydata, cdata = bootscdata, mdata = bootsmdata,
                      sigmau_inv = sigmau_inv, tbtheta = tbtheta, tL = tL, tU = tU,
                      nbreak = nbreak, p01 = p01, p02 = p02, j_max = j_max,
                      k_max = k_max, quadpoint = quadpoint, maxiter = maxiter,
                      do.trace = FALSE, beta0init = beta0init, beta1init = beta1init,
                      sigmainit = sigmainit, thetainit = thetainit, sigmadinit = sigmadinit,
                      gammainit = gammainit, survVar = survVar)
  
  if (fit$iter == maxiter) {
    totalp <- fit$TotalPara
    coef <- rep(NA, totalp+1)
    return(coef)
  } else {
    beta0 <- fit$beta0_matrix
    beta1 <- fit$beta1_estimate
    sigma2 <- fit$sigma2_estimate
    theta <- fit$theta_estimate
    sigmad <- fit$sigmad_estimate
    eta <- fit$eta_estimate
    btheta <- fit$btheta_matrix
    gamma <- fit$gamma
    
    pp = 1
    for (t in 1:j_max) {
      for (u in 1:p_max) {
        coef[pp] = beta0[t, u]
        pp = pp+1
      }
    }
    
    for (t in 1:j_max) {
      coef[pp] = beta1[t]
      pp = pp+1
    }
    
    for (t in 1:j_max) {
      coef[pp] = sigma2[t]
      pp = pp+1
    }
    
    for (t in 1:q_b) {
      coef[pp] = theta[t]
      pp = pp+1
    }
    
    for (t in 1:k_max) {
      coef[pp] = sigmad[t]
      pp = pp+1
    }
    
    for (t in 1:q_eta) {
      coef[pp] = eta[t]
      pp = pp+1
    }
    
    for (t in 1:k_max) {
      for (u in 1:q_b) {
        coef[pp] = btheta[u, t]
        pp = pp+1
      }
    }
    
    coef[pp] = gamma
    
    return(coef)
    
  }
}