##' @title Standard error estimates by bootstrap sampling
##' @param object an object inheriting from class JMM.
##' @param nboots number of bootstrap samples for standard error estimation. 
##' Default is 100.
##' @param print.bootspara logical; if TRUE a matrix of parameter estimates from 
##' all bootstrap samples will be returned.
##' @param maxiter number of iterations for EM algorithm. Default is 1000.
##' @param ncores number of cores to proceed parallel computation.
##' @param quadpoint number of Gauss-Hermite quadrature points used for numerical integration. 
##' @export
##'

bootsSE <- function(object, nboots = 100, print.bootspara = TRUE, maxiter = 1000,
                      ncores = 2, quadpoint = 10) {
  if (!inherits(object, "JMM"))
    stop("Use only with 'JMM' objects.\n")
  
  ydata <- as.data.frame(object$ydata)
  mdata <- as.data.frame(object$N1)
  cdata <- as.data.frame(object$cdata)
  
  ##pass data datas to jmspline
  ydatanew=tempfile(pattern = "", fileext = ".txt")
  writenh(ydata,ydatanew)
  mdatanew=tempfile(pattern = "", fileext = ".txt")
  writenh(mdata,mdatanew)
  cdatanew=tempfile(pattern = "", fileext = ".txt")
  writenh(cdata,cdatanew)
  
  sigmau_inv <- object$sigmau_inv
  tbtheta <- object$tbtheta
  
  tL <- object$tL
  tU <- object$tU
  
  nbreak <- object$nbreak
  q_b <- 4 + nbreak - 2
  p01 <- object$p01
  p02 <- object$p02
  
  p_max <- max(p01, p02)
  
  j_max <- object$j_max
  k_max <- object$k_max
  
  n <- object$n
  n1 <- object$n1
  
  cdim <- dim(cdata)
  
  q_eta = cdim[2] - 2
  
  t_max = max(mdata[, 1])
  
  Fnboots <- nboots*2
  
  Data = bootsdata_main(n, n1, tL, tU, q_eta, j_max, p01, p02, t_max, 
                            ydatanew, cdatanew, mdatanew, Fnboots)
  
  nycol <- dim(object$ydata)[2]
  nccol <- dim(object$cdata)[2]
  
  ## Initialize the parameter estimates
  tbeta0 <- as.data.frame(matrix(c(-0.01, -0.02), nrow = 2, ncol = 1))
  tbeta1 <- as.data.frame(c(1, 1.26))
  tsigma2 <- as.data.frame(c(3.2, 4.2))
  ttheta <- as.data.frame(c(15, -0.3, -7, 0.2, -4.7, 0.1, 1.4, -0.2, -0.2, 0))
  tsigmad <- as.data.frame(c(10.8, 6.4))
  tgamma <- 0.26
  
  if (ncores<1 || ncores > nboots) {
    stop("The specification of number of cores for parallel computation is not appropriate. 
         Please try another number.")
  }
  
  if (ncores == 1)
  {
    ## Allocate a matrix for bootstrap sample estimates
    ncolM <- fit$TotalPara + 1
    ParaMatrix <- matrix(NA, nrow = Fnboots, ncol = ncolM)
    ParaMatrix <- as.data.frame(ParaMatrix)
    Realboot <- 0
    for (i in 1:Fnboots) {
      writeLines(paste0("Try ", i, " th sample now!"))
      bootsydata <- matrix(unlist(Data[1+(i-1)*3]), ncol = nycol)
      bootsydata <- as.data.frame(bootsydata)
      bootscdata <- matrix(unlist(Data[2+(i-1)*3]), ncol = nccol)
      bootscdata <- as.data.frame(bootscdata)
      bootsmdata <- matrix(unlist(Data[3*i]), ncol = 1)
      bootsmdata <- as.data.frame(bootsmdata)

      fit <- jmspline(ydata = bootsydata, cdata = bootscdata, mdata = bootsmdata,
                          sigmau_inv = sigmau_inv, tbtheta = tbtheta, tL = tL, tU = tU,
                          nbreak = nbreak, p01 = p01, p02 = p02, j_max = j_max,
                          k_max = k_max, quadpoint = quadpoint, maxiter = maxiter,
                          do.trace = FALSE, beta0init = tbeta0, beta1init = tbeta1,
                          sigmainit = tsigma2, thetainit = ttheta, sigmadinit = tsigmad,
                          gammainit = tgamma)
      
      if (fit$iter == maxiter) {
        ParaMatrix[i, ] <- NA
      } else {
        Realboot <- Realboot + 1
        writeLines(paste0(Realboot, " th sample's parameter estimates is collected!"))
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
        
      }
      
      if (Realboot == nboots) break
      
    }
    ParaMatrix <- ParaMatrix[complete.cases(ParaMatrix), ] 
    
  } else {
    ParaMatrixRaw <- parallel::mclapply(1:nboots, bootsfit, Data = Data, 
                                        nycol = nycol, nccol = nccol, 
                                        sigmau_inv = sigmau_inv, 
                                        tbtheta = tbtheta, tL = tL, tU = tU,
                                        nbreak = nbreak, p01 = p01, p02 = p02, 
                                        j_max = j_max,
                                        k_max = k_max, quadpoint = quadpoint, 
                                        maxiter = maxiter,
                                        beta0init = tbeta0, beta1init = tbeta1,
                                        sigmainit = tsigma2, thetainit = ttheta, 
                                        sigmadinit = tsigmad,
                                        gammainit = tgamma, mc.cores = ncores)
    
    ParaMatrix <- t(matrix(unlist(ParaMatrixRaw), nrow = (object$TotalPara + 1)))
    ParaMatrix <- ParaMatrix[complete.cases(ParaMatrix), ]
    
    FrowPara <- nrow(ParaMatrix)
    u <- FrowPara
    t <- nboots
    
    while (u < nboots && t < Fnboots) {
      nncores <- min((nboots - u), ncores)
      ParaMatrixRaw <- parallel::mclapply((t+1):(t + nboots - u), bootsfit, 
                                          Data = Data, 
                                          nycol = nycol, nccol = nccol, 
                                          sigmau_inv = sigmau_inv, 
                                          tbtheta = tbtheta, tL = tL, tU = tU,
                                          nbreak = nbreak, p01 = p01, p02 = p02, 
                                          j_max = j_max,
                                          k_max = k_max, quadpoint = quadpoint, 
                                          maxiter = maxiter,
                                          beta0init = tbeta0, beta1init = tbeta1,
                                          sigmainit = tsigma2, thetainit = ttheta, 
                                          sigmadinit = tsigmad,
                                          gammainit = tgamma, mc.cores = nncores)
      
      SubParaMatrix <- t(matrix(unlist(ParaMatrixRaw), nrow = (object$TotalPara + 1)))
      
      if (sum(is.na(SubParaMatrix)) == nrow(SubParaMatrix)*(object$TotalPara + 1)) {
        FrowPara <- 0
        u <- u + FrowPara
        t <- t + nboots - u
      } else if (nrow(SubParaMatrix) == 1) {
        FrowPara <- 1
        u <- u + FrowPara
        t <- t + nboots - u
        ParaMatrix <- rbind(ParaMatrix, SubParaMatrix)
      } else {
        SubParaMatrix <- SubParaMatrix[complete.cases(SubParaMatrix), ]
        FrowPara <- nrow(SubParaMatrix)
        u <- u + FrowPara
        t <- t + nboots - u
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
    
  }
  
  
  ##standard error:
  SEpara <- matrix(sapply(ParaMatrix[, c(1:(j_max*p_max), (j_max*p_max + 2):(j_max*p_max + j_max),
                                   (j_max*p_max + j_max*2 + q_b + k_max + 1):
                                     (j_max*p_max + j_max*2 + q_b + k_max + q_eta), 
                                   ncol(ParaMatrix))], sd, 2), nrow = 1)
  ##log standard error
  SElogpara <- vector()
  index <- c((j_max*p_max+j_max+1):(j_max*p_max+2*j_max), 
             (j_max*p_max+2*j_max+q_b+1):(j_max*p_max+2*j_max+q_b+k_max))
  count <- 1
  for (i in index) {
    SElogpara[count] <- sd(log(ParaMatrix[, i]))
    count <- count + 1
  }
  SElogpara <- matrix(SElogpara, nrow = 1)
  
  Table2 <- as.data.frame(matrix(0, nrow = (j_max*p_max + (j_max-1) + j_max + k_max + q_eta + 1), 
                                 ncol = 5))
  rownames(Table2) <- colnames(ParaMatrix)[c(1:(j_max*p_max), (j_max*p_max + 2):(j_max*p_max + j_max),
                                       (j_max*p_max + j_max + 1):(j_max*p_max + j_max*2),
                                       (j_max*p_max + j_max*2 + q_b + 1):
                                         (j_max*p_max + j_max*2 + q_b + k_max),
                                       (j_max*p_max + j_max*2 + q_b + k_max + 1):
                                         (j_max*p_max + j_max*2 + q_b + k_max + q_eta), 
                                       ncol(ParaMatrix))]
  pp <- 1
  for (t in 1:j_max) {
    for (u in 1:p_max) {
      Table2[pp, 1] = object$beta0_matrix[t, u]
      Table2[pp, 2] = SEpara[1, pp]
      pp = pp + 1
    }
  }
  
  for (t in 2:j_max) {
    Table2[pp, 1] = object$beta1_estimate[t]
    Table2[pp, 2] = SEpara[1, pp]
    pp = pp + 1
  }
  
  uu <- 1
  for (t in 1:j_max) {
    Table2[pp, 1] = object$sigma2_estimate[t]
    Table2[pp, 2] = SElogpara[1, uu]
    pp = pp + 1
    uu = uu + 1
  }
  
  for (t in 1:k_max) {
    Table2[pp, 1] = object$sigmad_estimate[t]
    Table2[pp, 2] = SElogpara[1, uu]
    pp = pp + 1
    uu = uu + 1
  }
  
  uu = pp - j_max - k_max
  for (t in 1:q_eta) {
    Table2[pp, 1] = object$eta_estimate[t]
    Table2[pp, 2] = SEpara[1, uu]
    pp = pp + 1
    uu = uu + 1
  }
  
  Table2[pp, 1] = object$gamma
  Table2[pp, 2] = SEpara[1, uu]
  
  Table2[, 3] <- Table2[, 1] - 1.96*Table2[, 2]
  Table2[, 4] <- Table2[, 1] + 1.96*Table2[, 2]
  
  index <- c((j_max*p_max + j_max):(j_max*p_max + j_max*2 + k_max - 1))
  for (i in index) {
    Table2[i, 3] <- exp(log(Table2[i, 1]) - 1.96*Table2[i, 2])
    Table2[i, 4] <- exp(log(Table2[i, 1]) + 1.96*Table2[i, 2])
    Table2[i, 5] <- NA
  }
  
  for (i in 1:nrow(Table2)) {
    if(!is.na(Table2[i, 5])) {
      zval = (Table2[i, 1]/Table2[i, 2])
      Table2[i, 5] = 2 * pnorm(-abs(zval))
    }
  }
  
  Table2 <- round(Table2, 3)
  colnames(Table2) <- c("coef", "SE", "95%Lower", "95%Upper", "p-value")
  
  if (print.bootspara == TRUE) {
    Result <- list(ParaMatrix, Table2)
    names(Result) <- c("ParaMatrix", "Table")
    return(Result)
  } else {
    return(Table2)
  }
  
}