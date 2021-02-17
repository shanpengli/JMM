#ifndef bootsdata_hpp
#define bootsdata_hpp

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

namespace bootsdataspace {

Rcpp::List bootsdata_cmain(int n, int n1, double tL, double tU, int q_eta, 
                           int j_max, int p01, int p02, int t_max,
                           std::string ydatanew, std::string cdatanew,
                           std::string mdatanew, int nboots);


}


#endif /* bootsdata_hpp */