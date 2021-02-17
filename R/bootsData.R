##' @export
##'

bootsData <- function(object, nboots = 100) {
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
  return(Data)
  
}