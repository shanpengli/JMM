#ifndef Simdata_hpp
#define Simdata_hpp

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
#include <Rcpp.h>

using namespace Rcpp;

namespace Simdataspace {

  double MulVV(const gsl_vector *Z,const gsl_vector *beta);
  
  void MulV(const gsl_vector *Z,gsl_matrix *ZZ);
  
  void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta);
  
  void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB);
  
  
  double GETT(
      const gsl_rng *r,
      const int i,
      const double tl,
      const double tu,
      const int tnbreak,
      gsl_vector *B_spl,
      const gsl_vector *theta,
      const gsl_matrix *btheta,
      const gsl_vector *alphai,
      const double lambda0,
      const double gamma,
      const gsl_vector *eta,
      const gsl_matrix *C,
      const gsl_matrix *sigmau_inv,
      const int k_cubic);
  
  double Min(const double t1, const double t2);
  
  int GetN(double t);
  
  double GetMU(const gsl_vector *theta, 
               const gsl_matrix *btheta, 
               const gsl_vector *alphai, 
               gsl_vector *B_spl,
               const gsl_matrix *sigmau_inv);



  Rcpp::List Simdata_cmain(int n, int sim, int nbreak, double tL, double tU, int q_eta, 
                           int j_max, int p01, int p02, int t_max, int distr, double m_age,
                           double std_age, int k_max,
                           std::string tbthetanew, std::string sigmau_invnew,
                           std::string tbeta0new, std::string tbeta1new,
                           std::string tsigmanew, std::string tthetanew,
                           std::string tsigmadnew, std::string tetanew,
                           double tgamma, double lambda0);


}


#endif /* Simdata_hpp */