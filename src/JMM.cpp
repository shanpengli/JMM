//
//  JMM.cpp
//  JMM

#include <stdio.h>
#include "jmspline.hpp"
#include "bootsdata.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List  jmspline_main(SEXP n, SEXP n_total, SEXP tL, SEXP tU, SEXP p01, 
                          SEXP p02, SEXP q_b, SEXP q_eta, SEXP j_max, 
                          SEXP t_max, SEXP nbreak, SEXP k_max, SEXP quadpoint, 
                          SEXP maxiter, SEXP trace, SEXP ydatanew, 
                          SEXP mdatanew, SEXP cdatanew, SEXP sigmau_invnew, 
                          SEXP tbthetanew, SEXP xs, SEXP ws, SEXP beta0initnew, 
                          SEXP beta1initnew, SEXP sigmainitnew, SEXP thetainitnew, 
                          SEXP sigmadinitnew, SEXP gammainit)
{
  Rcpp::List result;
  try {
    
    result=jmsplinespace::jmspline_cmain(as<int> (n), as<int> (n_total), 
                                         as<double> (tL), as<double> (tU), 
                                         as<int> (p01), as<int> (p02), 
                                         as<int>(q_b), as<int>(q_eta),
                                         as<int>(j_max), as<int>(t_max),
                                         as<int>(nbreak), as<int>(k_max), 
                                         as<int>(quadpoint), as<int>(maxiter),
                                         as<int> (trace), 
                                         as<std::string> (ydatanew),
                                         as<std::string> (mdatanew),
                                         as<std::string>(cdatanew),
                                         as<std::string>(sigmau_invnew),
                                         as<std::string>(tbthetanew),
                                         as<std::vector<double> >(xs),
                                         as<std::vector<double> >(ws),
                                         as<std::string> (beta0initnew),
                                         as<std::string> (beta1initnew),
                                         as<std::string>(sigmainitnew),
                                         as<std::string>(thetainitnew),
                                         as<std::string>(sigmadinitnew),
                                         as<double> (gammainit));
    if(Rf_isNull(result)){
      throw std::range_error("Possible files reading or format errors");
    }
    return result;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return R_NilValue;             // not reached
  
}


// [[Rcpp::export]]
Rcpp::List  bootsdata_main(SEXP n, SEXP n1, SEXP tL, SEXP tU, SEXP q_eta, SEXP j_max, 
                           SEXP p01, SEXP p02, SEXP t_max, 
                           SEXP ydatanew, SEXP cdatanew, SEXP mdatanew, SEXP nboots)
{
  Rcpp::List result;
  try {
    
    result=bootsdataspace::bootsdata_cmain(as<int> (n), as<int> (n1), 
                                         as<double> (tL), as<double> (tU),
                                         as<int>(q_eta), as<int>(j_max),
                                         as<int> (p01), as<int> (p02),
                                         as<int>(t_max),
                                         as<std::string> (ydatanew),
                                         as<std::string> (cdatanew),
                                         as<std::string>(mdatanew),
                                         as<int> (nboots));
    if(Rf_isNull(result)){
      throw std::range_error("Possible files reading or format errors");
    }
    return result;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return R_NilValue;             // not reached
  
}