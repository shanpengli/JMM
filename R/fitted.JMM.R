##' @export
##'

fitted.JMM <- function(object) {
  if (!inherits(object, "JMM"))
    stop("Use only with 'JMM' objects.\n")
  
  ydata <- object$ydata
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
                             sigmau_invnew, thetanew, bthetanew, beta0new, beta1new)
  
  return(myresult)
  
  
}