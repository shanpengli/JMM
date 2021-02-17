##' @export
##'

bootsSE <- function(object, nboots = 100, print.para = FALSE, maxiter = 1000,
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
  
  Data = bootsdata_main(n, n1, tL, tU, q_eta, j_max, p01, p02, t_max, 
                            ydatanew, cdatanew, mdatanew, nboots)
  
  nycol <- dim(object$ydata)[2]
  nccol <- dim(object$cdata)[2]
  
  ## Allocate a matrix for bootstrap sample estimates
  ncolM <- j_max*p_max + 2*j_max + q_b + k_max + q_eta + k_max*q_b + 1
  ParaMatrix <- matrix(0, nrow = nboots, ncol = ncolM) 
  
  ## Initialize the parameter estimates
  tbeta0 <- as.data.frame(matrix(c(-0.01, -0.02), nrow = 2, ncol = 1))
  tbeta1 <- as.data.frame(c(1, 1.26))
  tsigma2 <- as.data.frame(c(3.2, 4.2))
  ttheta <- as.data.frame(c(15, -0.3, -7, 0.2, -4.7, 0.1, 1.4, -0.2, -0.2, 0))
  tsigmad <- as.data.frame(c(10.8, 6.4))
  tgamma <- 0.26
  
  
  
  
  if (ncores == 1)
  {
    for (i in 1:nboots) {
      writeLines(paste0("Running ", i, " th sample!"))
      bootsydata <- matrix(unlist(Data[1+(i-1)*3]), ncol = nycol)
      bootsydata <- as.data.frame(bootsydata)
      bootscdata <- matrix(unlist(Data[2+(i-1)*3]), ncol = nccol)
      bootscdata <- as.data.frame(bootscdata)
      bootsmdata <- matrix(unlist(Data[3*i]), ncol = 1)
      bootsmdata <- as.data.frame(bootsmdata)

      fit <- try(jmspline(ydata = bootsydata, cdata = bootscdata, mdata = bootsmdata,
                          sigmau_inv = sigmau_inv, tbtheta = tbtheta, tL = tL, tU = tU,
                          nbreak = nbreak, p01 = p01, p02 = p02, j_max = j_max,
                          k_max = k_max, quadpoint = quadpoint, maxiter = maxiter,
                          do.trace = FALSE, beta0init = tbeta0, beta1init = tbeta1,
                          sigmainit = tsigma2, thetainit = ttheta, sigmadinit = tsigmad,
                          gammainit = tgamma), silent = TRUE)
      
      if ('try-error' %in% class(fit)) {
        ParaMatrix[i, ] <- NA
      } else {
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
                pp = pp+1
              }
            }
            
            for (t in 1:j_max) {
              ParaMatrix[i, pp] = beta1[t]
              pp = pp+1
            }
            
            for (t in 1:j_max) {
              ParaMatrix[i, pp] = sigma2[t]
              pp = pp+1
            }
            
            for (t in 1:q_b) {
              ParaMatrix[i, pp] = theta[t]
              pp = pp+1
            }
            
            for (t in 1:k_max) {
              ParaMatrix[i, pp] = sigmad[t]
              pp = pp+1
            }
            
            for (t in 1:q_eta) {
              ParaMatrix[i, pp] = eta[t]
              pp = pp+1
            }
            
            for (t in 1:k_max) {
              for (u in 1:q_b) {
                ParaMatrix[i, pp] = btheta[u, t]
                pp = pp+1
              }
            }
            
            ParaMatrix[i, pp] = gamma
            
            if (print.para == TRUE) {
              print(ParaMatrix[i, ])
            }
        
      }
    }
    return(ParaMatrix)
  } else  {
    ParaMatrixRaw <- parallel::mclapply(1:nboots, bootsfit, mc.cores = ncores)
    
    return(ParaMatrixRaw)
  }
  
  
  
}