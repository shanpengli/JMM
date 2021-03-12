#ifndef getfitted_hpp
#define getfitted_hpp

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

namespace getfittedspace {

double MulVV(const gsl_vector *Z,const gsl_vector *beta);

void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta);

Rcpp::List getfitted_cmain(double tL, double tU, int nbreak, int k_max, int j_max, 
                           int p01, int p02,
                           std::string sigmau_invnew, 
                           std::string thetanew, 
                           std::string bthetanew,
                           std::string beta0new,
                           std::string beta1new,
                           std::string longvarnew);

}

#endif /* getfitted_hpp */
