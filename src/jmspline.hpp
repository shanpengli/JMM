//
//  jmspline.hpp
//  JMM

#ifndef jmspline_hpp
#define jmspline_hpp

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <Rcpp.h>
using namespace Rcpp;

namespace jmsplinespace {
  
  int gsl_linalg_cholesky_decompn (gsl_matrix * A);
  
  int gsl_linalg_cholesky_solven (const gsl_matrix * LLT,
                                  const gsl_vector * b,
                                  gsl_vector * x);
  
  int inv_matrix(gsl_matrix *x_square);
  
  double HAZ(const gsl_matrix *H, const double t);
  
  double CH(const gsl_matrix *H, const double t);
  
  double MulVV(const gsl_vector *Z,const gsl_vector *beta);
  
  void MulV(const gsl_vector *Z,gsl_matrix *ZZ);
  
  void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta);
  
  void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB);
  
  void TransM(const gsl_matrix *A, gsl_matrix *B);
  
  double Abs(const double a, const double b);
  
  double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb);
  
  double DiffV(const gsl_vector *veca, const gsl_vector *vecb);
  
  double Min(const double t1, const double t2);
  
  void STAT(const gsl_matrix *store,int i,double *mean,double *sd);
  
  int GetN(double t);
  
  
  double SmoothHAZ(const gsl_matrix *H, double t);
  
  
  
  int EM(
      gsl_matrix *beta0,
      gsl_vector *beta1,
      gsl_vector *sigma,
      gsl_vector *theta,
      gsl_vector *sigmad,
      gsl_vector *eta,
      gsl_matrix *btheta,
      gsl_matrix *H01,
      double *gamma,
      const gsl_vector *pdim,
      const gsl_matrix *Bio,
      const gsl_matrix *C,
      const gsl_vector *M,
      gsl_bspline_workspace *bw,
      const gsl_matrix *sigmau_inv,
      const int quadpoint,
      const std::vector<double> xs,
      const std::vector<double> ws
  );
  
  
  
  int GetE(
      gsl_matrix *FUNA,
      gsl_matrix *FUNA2,
      gsl_matrix *FUNE,
      gsl_matrix *FUNAE,
      gsl_matrix *FUNA2E,
      gsl_matrix *FUNAEV,
      gsl_matrix *FUNA2EV,
      const gsl_matrix *beta0,
      const gsl_vector *beta1,
      const gsl_vector *sigma,
      const gsl_vector *theta,
      const gsl_vector *sigmad,
      const gsl_vector *eta,
      const gsl_matrix *btheta,
      const gsl_matrix *H01,
      const double gamma,
      const gsl_vector *pdim,
      const gsl_matrix *Bio,
      const gsl_matrix *C,
      const gsl_vector *M,
      gsl_bspline_workspace *bw,
      const gsl_matrix *sigmau_inv,
      const int quadpoint,
      const std::vector<double> xs,
      const std::vector<double> ws
  );
  
  
  double GetLogLike(
      const gsl_matrix *beta0,
      const gsl_vector *beta1,
      const gsl_vector *theta,
      const gsl_vector *sigmad,
      const gsl_vector *eta,
      const gsl_matrix *btheta,
      const gsl_matrix *H01,
      const gsl_vector *sigma,
      const double gamma,
      const gsl_vector *pdim,
      const gsl_matrix *Bio,
      const gsl_matrix *C,
      const gsl_vector *M,
      gsl_bspline_workspace *bw,
      const gsl_matrix *sigmau_inv,
      const int quadpoint,
      const std::vector<double> xs,
      const std::vector<double> ws
  );
  
  int Diffn(
      const gsl_matrix *prebeta0,
      const gsl_matrix *beta0,
      const gsl_vector *prebeta1,
      const gsl_vector *beta1,
      const gsl_vector *presigma,
      const gsl_vector *sigma,
      const gsl_vector *pretheta,
      const gsl_vector *theta,
      const gsl_vector *preeta,
      const gsl_vector *eta,
      const double pregamma,
      const double gamma
  );



  //declare before use it
  Rcpp::List jmspline_cmain(int n, int n_total, double tL, double tU,int p01, 
                            int p02, int q_b, int q_eta, int j_max, int t_max,
                            int nbreak, int k_max, int quadpoint, int maxiter, 
                            int trace, std::string ydatanew, std::string mdatanew,
                            std::string cdatanew, std::string sigmau_invnew,
                            std::string tbthetanew, 
                            std::vector<double> xs, 
                            std::vector<double> ws);

}


#endif /* jmspline_hpp */