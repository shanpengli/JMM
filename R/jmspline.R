jmspline <- function(ydata, cdata, mdata, sigmau_inv, tbtheta, 
                     tL, tU, nbreak, p01, p02, j_max, k_max = 2,
                     point = 10, max.iter = 2500, do.trace = FALSE) {
  
  if (do.trace) {
    trace=1;
  } else {
    trace=0;
  }
  
  #Gaussian-Hermite quadrature nodes and weights
  #The dimension of xs/ws is half of the point value since they are symmetric
  
  if (point %% 2 == 1)
  {
    stop("Number of quadrature points can only be even!")
  }
  
  gq_vals <- statmod::gauss.quad(n = point, kind = "hermite")
  
  xs <- gq_vals$nodes[(point / 2 + 1) : point]
  
  ws <- gq_vals$weights[(point / 2 + 1) : point]
  
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
  ydatanew=tempdata(pattern = "", dataext = ".txt")
  writenh(ydata,ydatanew)
  mdatanew=tempdata(pattern = "", dataext = ".txt")
  writenh(mdata,mdatanew)
  cdatanew=tempdata(pattern = "", dataext = ".txt")
  writenh(cdata,cdatanew)
  sigmau_invnew=tempdata(pattern = "", dataext = ".txt")
  writenh(sigmau_inv,sigmau_invnew)
  tbthetanew=tempdata(pattern = "", dataext = ".txt")
  writenh(tbtheta,tbthetanew)
  
  myresult = jmspline_main(n, )
  
  
  
  
  
}