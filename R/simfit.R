##' @title Simulation results 
##' @param sim number of simulated datasets to be generated.
##' @param n sample size for each dataset.
##' @param tL lower limit of time point for spline basis.
##' @param tU upper limit of time point for spline basis.
##' @param nbreak number of knots for spline basis for a model fit.
##' @param p01 number of covariates in the first biomarker.
##' @param p02 number of covariates in the second biomarker.
##' @param j_max Total number of biomarkers in the longitudinal data.
##' @param quadpoint Number of quadrature points.
##' @param maxiter The maximum number of iterations for EM process.
##' @param distr distributional assumption of random effects.
##' @param m_age true mean value of age variable in the longitudinal sub-model.
##' @param std_age true standard deviation of age variable in the longitudinal sub-model. 
##' @param k_max dimension of random effects. 
##' @param q_eta number of survival covariates.
##' @param t_max maximum number of repeated measurements allowed for each subject.
##' @param tsigma_inv true sigma_inv matrix.
##' @param tnbreak number of knots for spline basis for data generation.
##' @param tbtheta true btheta matrix.
##' @param tbeta0 true value of beta0 matrix.
##' @param tbeta1 true value of beta1 vector.
##' @param tsigma true value of sigma square vector.
##' @param tsigmad true value of variance of random effects.
##' @param ttheta true value of theta vector.
##' @param teta true value of survival fixed effects.
##' @param tgamma true value of latent association parameter. 
##' @param sigma_invinit Initial values of sigma_inv matrix.
##' @param bthetainit a data frame of the normalizing constant for spline basis.
##' @param beta0init Initial values of beta0. 
##' @param beta1init Initial values of beta1. 
##' @param sigmainit Initial values of sigma^2. 
##' @param thetainit Initial values of theta. 
##' @param sigmadinit Initial values of sigmad. 
##' @param gammainit Initial values of gamma.
##' @param survVar logical; TRUE if the survival sub-model include covariates. Default is FALSE. 
##' @param lambda0 true baseline hazard.
##' @param out.bs logical; TRUE if output baseline hazard. Default is FALSE.
##' @param conversigmad logical; TRUE if sigmad is considered into convergence criteria. Default is FALSE.
##' @param ncores number of cores to proceed parallel computation.
##' @seealso \code{\link{SimData}, \link{jmspline}}
##' @export
##'

Simfit <- function(sim = 100, n = 215, tL = 0, tU = 9,
                   nbreak = 8, p01 = 1, p02 = 1, j_max = 2, quadpoint = 10,
                   maxiter = 4000, distr = "Gaussian",
                   m_age = 0, std_age = 2, k_max = 2, q_eta = 1, t_max = 19,
                   tsigmau_inv = NULL, tnbreak = 8,
                   tbtheta = NULL, tbeta0 = NULL,
                   tbeta1 = NULL, tsigma = NULL, tsigmad = NULL,
                   ttheta = NULL, teta = NULL, tgamma = NULL, sigmau_invinit = NULL,
                   bthetainit = NULL, beta0init = NULL, 
                   beta1init = NULL, sigmainit = NULL, thetainit = NULL,
                   sigmadinit = NULL, gammainit = NULL, survVar = FALSE,
                   lambda0 = NULL, out.bs = FALSE, conversigmad = FALSE, ncores = 2) {
  
  Fsim <- sim*2
  
  a <- SimData(sim = Fsim, n = n, tL = tL, tU = tU, sigmau_inv = tsigmau_inv,
               tbtheta = tbtheta, tbeta0 = tbeta0, tbeta1 = tbeta1,
               tsigma = tsigma, tsigmad = tsigmad, ttheta = ttheta,
               teta = teta, tgamma = tgamma, nbreak = tnbreak, distr = distr, 
               m_age = m_age, std_age = std_age, k_max = k_max, q_eta = q_eta, 
               t_max = t_max, lambda0 = lambda0)
  
  if (ncores<1 || ncores > sim) {
    stop("The specification of number of cores for parallel computation is not appropriate. 
         Please try another number.")
  }
  
  k_cubic = 4
  p_max <- max(p01, p02)
  q_b <- k_cubic + nbreak - 2
  ncolM <- j_max*p_max + 2*j_max + q_b + k_max + q_eta + k_max*q_b + 1
  
  nycol <- dim(a$ydata_0)[2]
  nccol <- dim(a$cdata_0)[2]
  
  BHMatrix <- NULL
  if (ncores == 1) {
    ## Allocate a matrix for simulated estimates
    ParaMatrix <- matrix(NA, nrow = Fsim, ncol = ncolM)
    ParaMatrix <- as.data.frame(ParaMatrix)
    Realsim <- 0
    for (i in 1:Fsim) {
      writeLines(paste0("Try ", i, " th sample now!"))
      simydata <- matrix(unlist(a[1+(i-1)*3]), ncol = nycol)
      simydata <- as.data.frame(simydata)
      simcdata <- matrix(unlist(a[2+(i-1)*3]), ncol = nccol)
      simcdata <- as.data.frame(simcdata)
      simmdata <- matrix(unlist(a[3*i]), ncol = 1)
      simmdata <- as.data.frame(simmdata)
      
      fit <- jmspline(ydata = simydata, cdata = simcdata, mdata = simmdata,
                      sigmau_inv = sigmau_invinit, tbtheta = bthetainit, tL = tL, tU = tU,
                      nbreak = nbreak, p01 = p01, p02 = p02, j_max = j_max,
                      k_max = k_max, quadpoint = quadpoint, maxiter = maxiter,
                      do.trace = FALSE, beta0init = beta0init, beta1init = beta1init,
                      sigmainit = sigmainit, thetainit = thetainit, sigmadinit = sigmadinit,
                      gammainit = gammainit, survVar = survVar, conversigmad = conversigmad)
      
      if (fit$iter == maxiter) {
        ParaMatrix[i, ] <- NA
      } else if (fit$loglike == -Inf) {
        ParaMatrix[i, ] <- NA
      } else {
        Realsim <- Realsim + 1
        writeLines(paste0(Realsim, " th sample's parameter estimates is collected!"))
        
        beta0 <- fit$beta0_matrix
        beta1 <- fit$beta1_estimate
        sigma2 <- fit$sigma2_estimate
        theta <- fit$theta_estimate
        sigmad <- fit$sigmad_estimate
        eta <- fit$eta_estimate
        btheta <- fit$btheta_matrix
        gamma <- fit$gamma
        
        ##ncolM <- j_max*p_max + 2*j_max + q_b + k_max + q_eta + k_max*q_b + 1
        pp = 1
        for (t in 1:j_max) {
          for (u in 1:p_max) {
            ParaMatrix[i, pp] = beta0[t, u]
            colnames(ParaMatrix)[pp] <- paste0("beta0_", t, u)
            pp = pp+1
          }
        }
        
        for (t in 1:j_max) {
          ParaMatrix[i, pp] = beta1[t]
          colnames(ParaMatrix)[pp] <- paste0("beta1_", t)
          pp = pp+1
        }
        
        for (t in 1:j_max) {
          ParaMatrix[i, pp] = sigma2[t]
          colnames(ParaMatrix)[pp] <- paste0("sigma2_", t)
          pp = pp+1
        }
        
        for (t in 1:q_b) {
          ParaMatrix[i, pp] = theta[t]
          colnames(ParaMatrix)[pp] <- paste0("theta_", t)
          pp = pp+1
        }
        
        for (t in 1:k_max) {
          ParaMatrix[i, pp] = sigmad[t]
          colnames(ParaMatrix)[pp] <- paste0("sigmad_", t)
          pp = pp+1
        }
        
        for (t in 1:q_eta) {
          ParaMatrix[i, pp] = eta[t]
          colnames(ParaMatrix)[pp] <- paste0("eta_", t)
          pp = pp+1
        }
        
        for (t in 1:k_max) {
          for (u in 1:q_b) {
            ParaMatrix[i, pp] = btheta[u, t]
            colnames(ParaMatrix)[pp] <- paste0("btheta_", u, t)
            pp = pp+1
          }
        }
        
        ParaMatrix[i, pp] = gamma
        colnames(ParaMatrix)[pp] <- "gamma"
        
        
        ##output baselinehazard
        if (out.bs == TRUE) {
          BS <- t(fit$BaselineHazard)
          BS <- cbind(BS, Realsim)
          BHMatrix <- rbind(BHMatrix, BS) 
        } 
        
      }
      
      if (Realsim == sim) break
      
    }
    ParaMatrix <- ParaMatrix[complete.cases(ParaMatrix), ]
    if (out.bs == TRUE) {
      BHMatrix <- BHMatrix[, -2]
      BHMatrix <- as.data.frame(BHMatrix)
    }
    ParaMatrix <- as.data.frame(ParaMatrix)
    a <- list(ParaMatrix, BHMatrix)
    names(a) <- c("ParaMatrix", "BHMatrix")
    return(a)
    
  } else {

    if (out.bs == TRUE) {
      # ParaMatrixRaw <- parallel::mclapply(1:sim, bootsfit, Data = a, 
      #                                     nycol = nycol, nccol = nccol, 
      #                                     sigmau_inv = sigmau_invinit, 
      #                                     tbtheta = bthetainit, tL = tL, tU = tU,
      #                                     nbreak = nbreak, p01 = p01, p02 = p02, 
      #                                     j_max = j_max,
      #                                     k_max = k_max, quadpoint = quadpoint, 
      #                                     maxiter = maxiter,
      #                                     beta0init = beta0init, beta1init = beta1init,
      #                                     sigmainit = sigmainit, thetainit = thetainit, 
      #                                     sigmadinit = sigmadinit,
      #                                     gammainit = gammainit,
      #                                     survVar = survVar,
      #                                     conversigmad = conversigmad,
      #                                     out.bs = out.bs,
      #                                     mc.cores = ncores)
      
      cl <- parallel::makeCluster(ncores)
      ParaMatrixRaw <- parallel::parLapply(cl, 1:sim, bootsfit, 
                                           Data = a, 
                                           nycol = nycol, nccol = nccol, 
                                           sigmau_inv = sigmau_invinit, 
                                           tbtheta = bthetainit, tL = tL, tU = tU,
                                           nbreak = nbreak, p01 = p01, p02 = p02, 
                                           j_max = j_max,
                                           k_max = k_max, quadpoint = quadpoint, 
                                           maxiter = maxiter,
                                           beta0init = beta0init, beta1init = beta1init,
                                           sigmainit = sigmainit, thetainit = thetainit, 
                                           sigmadinit = sigmadinit,
                                           gammainit = gammainit,
                                           survVar = survVar,
                                           conversigmad = conversigmad,
                                           out.bs = out.bs)
      parallel::stopCluster(cl)
      
      count <- 0
      ParaMatrix <- NULL
      BHMatrix <- NULL
      for (i in 1:length(ParaMatrixRaw)) {
        if (!is.na(ParaMatrixRaw[i][[1]]$coef1[1])) {
          count <- count + 1
          subParaMatrix <- ParaMatrixRaw[i][[1]]$coef1
          ParaMatrix <- rbind(ParaMatrix, subParaMatrix)
          subBHMatrix <- cbind(ParaMatrixRaw[i][[1]]$coef2, count)
          BHMatrix <- rbind(BHMatrix, subBHMatrix)
        }
      }
      
      FrowPara <- nrow(ParaMatrix)
      u <- FrowPara
      t <- sim
      
      while (u < sim && t < Fsim) {
        nncores <- min((sim - u), ncores)
        # ParaMatrixRaw <- parallel::mclapply((t+1):(t + sim - u), bootsfit, 
        #                                     Data = a, 
        #                                     nycol = nycol, nccol = nccol, 
        #                                     sigmau_inv = sigmau_invinit, 
        #                                     tbtheta = bthetainit, tL = tL, tU = tU,
        #                                     nbreak = nbreak, p01 = p01, p02 = p02, 
        #                                     j_max = j_max,
        #                                     k_max = k_max, quadpoint = quadpoint, 
        #                                     maxiter = maxiter,
        #                                     beta0init = beta0init, beta1init = beta1init,
        #                                     sigmainit = sigmainit, thetainit = thetainit, 
        #                                     sigmadinit = sigmadinit,
        #                                     gammainit = gammainit,
        #                                     survVar = survVar,
        #                                     conversigmad = conversigmad,
        #                                     out.bs = out.bs,
        #                                     mc.cores = nncores)
        # 
        cl <- parallel::makeCluster(nncores)
        ParaMatrixRaw <- parallel::parLapply(cl, (t+1):(t + sim - u), bootsfit, 
                                             Data = a, 
                                             nycol = nycol, nccol = nccol, 
                                             sigmau_inv = sigmau_invinit, 
                                             tbtheta = bthetainit, tL = tL, tU = tU,
                                             nbreak = nbreak, p01 = p01, p02 = p02, 
                                             j_max = j_max,
                                             k_max = k_max, quadpoint = quadpoint, 
                                             maxiter = maxiter,
                                             beta0init = beta0init, beta1init = beta1init,
                                             sigmainit = sigmainit, thetainit = thetainit, 
                                             sigmadinit = sigmadinit,
                                             gammainit = gammainit,
                                             survVar = survVar,
                                             conversigmad = conversigmad,
                                             out.bs = out.bs)
        parallel::stopCluster(cl)
        
        subcount <- 0
        
        for (i in 1:length(ParaMatrixRaw)) {
          if (!is.na(ParaMatrixRaw[i][[1]]$coef1[1])) {
            count <- count + 1
            subcount <- subcount + 1
            subParaMatrix <- ParaMatrixRaw[i][[1]]$coef1
            ParaMatrix <- rbind(ParaMatrix, subParaMatrix)
            subBHMatrix <- cbind(ParaMatrixRaw[i][[1]]$coef2, count)
            BHMatrix <- rbind(BHMatrix, subBHMatrix)
          }
        }
        
        if (subcount == 0) {
          FrowPara <- 0
          u <- u + FrowPara
          t <- t + sim - u
        } else {
          FrowPara <- subcount
          u <- u + FrowPara
          t <- t + sim - u
        }
      }
      ParaMatrix <- as.data.frame(ParaMatrix)
      BHMatrix <- as.data.frame(BHMatrix)
      BHMatrix <- BHMatrix[, -2]
      a <- list(ParaMatrix, BHMatrix)
      names(a) <- c("ParaMatrix", "BHMatrix")
      return(a)
    } else {
      # ParaMatrixRaw <- parallel::mclapply(1:sim, bootsfit, Data = a, 
      #                                     nycol = nycol, nccol = nccol, 
      #                                     sigmau_inv = sigmau_invinit, 
      #                                     tbtheta = bthetainit, tL = tL, tU = tU,
      #                                     nbreak = nbreak, p01 = p01, p02 = p02, 
      #                                     j_max = j_max,
      #                                     k_max = k_max, quadpoint = quadpoint, 
      #                                     maxiter = maxiter,
      #                                     beta0init = beta0init, beta1init = beta1init,
      #                                     sigmainit = sigmainit, thetainit = thetainit, 
      #                                     sigmadinit = sigmadinit,
      #                                     gammainit = gammainit,
      #                                     survVar = survVar,
      #                                     conversigmad = conversigmad,
      #                                     mc.cores = ncores)
      
      cl <- parallel::makeCluster(ncores)
      ParaMatrixRaw <- parallel::parLapply(cl, 1:sim, bootsfit, 
                                           Data = a, 
                                           nycol = nycol, nccol = nccol, 
                                           sigmau_inv = sigmau_invinit, 
                                           tbtheta = bthetainit, tL = tL, tU = tU,
                                           nbreak = nbreak, p01 = p01, p02 = p02, 
                                           j_max = j_max,
                                           k_max = k_max, quadpoint = quadpoint, 
                                           maxiter = maxiter,
                                           beta0init = beta0init, beta1init = beta1init,
                                           sigmainit = sigmainit, thetainit = thetainit, 
                                           sigmadinit = sigmadinit,
                                           gammainit = gammainit,
                                           survVar = survVar,
                                           conversigmad = conversigmad)
      parallel::stopCluster(cl)
      
      
      ParaMatrix <- t(matrix(unlist(ParaMatrixRaw), nrow = ncolM))
      ParaMatrix <- ParaMatrix[complete.cases(ParaMatrix), ]
      
      FrowPara <- nrow(ParaMatrix)
      u <- FrowPara
      t <- sim
      
      while (u < sim && t < Fsim) {
        nncores <- min((sim - u), ncores)
        # ParaMatrixRaw <- parallel::mclapply((t+1):(t + sim - u), bootsfit, 
        #                                     Data = a, 
        #                                     nycol = nycol, nccol = nccol, 
        #                                     sigmau_inv = sigmau_invinit, 
        #                                     tbtheta = bthetainit, tL = tL, tU = tU,
        #                                     nbreak = nbreak, p01 = p01, p02 = p02, 
        #                                     j_max = j_max,
        #                                     k_max = k_max, quadpoint = quadpoint, 
        #                                     maxiter = maxiter,
        #                                     beta0init = beta0init, beta1init = beta1init,
        #                                     sigmainit = sigmainit, thetainit = thetainit, 
        #                                     sigmadinit = sigmadinit,
        #                                     gammainit = gammainit,
        #                                     survVar = survVar,
        #                                     conversigmad = conversigmad,
        #                                     mc.cores = nncores)
        
        cl <- parallel::makeCluster(nncores)
        ParaMatrixRaw <- parallel::parLapply(cl, (t+1):(t + sim - u), bootsfit, 
                                             Data = a, 
                                             nycol = nycol, nccol = nccol, 
                                             sigmau_inv = sigmau_invinit, 
                                             tbtheta = bthetainit, tL = tL, tU = tU,
                                             nbreak = nbreak, p01 = p01, p02 = p02, 
                                             j_max = j_max,
                                             k_max = k_max, quadpoint = quadpoint, 
                                             maxiter = maxiter,
                                             beta0init = beta0init, beta1init = beta1init,
                                             sigmainit = sigmainit, thetainit = thetainit, 
                                             sigmadinit = sigmadinit,
                                             gammainit = gammainit,
                                             survVar = survVar,
                                             conversigmad = conversigmad)
        parallel::stopCluster(cl)
        
        SubParaMatrix <- t(matrix(unlist(ParaMatrixRaw), nrow = ncolM))
        
        if (sum(is.na(SubParaMatrix)) == nrow(SubParaMatrix)*ncolM) {
          FrowPara <- 0
          u <- u + FrowPara
          t <- t + sim - u
        } else if (nrow(SubParaMatrix) == 1) {
          FrowPara <- 1
          u <- u + FrowPara
          t <- t + sim - u
          ParaMatrix <- rbind(ParaMatrix, SubParaMatrix)
        } else {
          SubParaMatrix <- SubParaMatrix[complete.cases(SubParaMatrix), ]
          FrowPara <- nrow(SubParaMatrix)
          u <- u + FrowPara
          t <- t + sim - u
          ParaMatrix <- rbind(ParaMatrix, SubParaMatrix)
        }
      }
      ParaMatrix <- as.data.frame(ParaMatrix)
      
      ## name the parameters
      pp = 1
      for (t in 1:j_max) {
        for (u in 1:p_max) {
          colnames(ParaMatrix)[pp] <- paste0("beta0_", t, u)
          pp = pp+1
        }
      }
      
      for (t in 1:j_max) {
        colnames(ParaMatrix)[pp] <- paste0("beta1_", t)
        pp = pp+1
      }
      
      for (t in 1:j_max) {
        colnames(ParaMatrix)[pp] <- paste0("sigma2_", t)
        pp = pp+1
      }
      
      for (t in 1:q_b) {
        colnames(ParaMatrix)[pp] <- paste0("theta_", t)
        pp = pp+1
      }
      
      for (t in 1:k_max) {
        colnames(ParaMatrix)[pp] <- paste0("sigmad_", t)
        pp = pp+1
      }
      
      for (t in 1:q_eta) {
        colnames(ParaMatrix)[pp] <- paste0("eta_", t)
        pp = pp+1
      }
      
      for (t in 1:k_max) {
        for (u in 1:q_b) {
          colnames(ParaMatrix)[pp] <- paste0("btheta_", u, t)
          pp = pp+1
        }
      }
      colnames(ParaMatrix)[pp] <- "gamma"
      a <- list(ParaMatrix, BHMatrix)
      names(a) <- c("ParaMatrix", "BHMatrix")
      return(a)
      }
    
  }
  
}