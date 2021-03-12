//
//  JMM.cpp
//  JMM

#include <stdio.h>
#include "jmspline.hpp"
#include "bootsdata.hpp"
#include "Simdata.hpp"
#include "getfitted.hpp"
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
                          SEXP sigmadinitnew, SEXP gammainit, SEXP survvar, 
                          SEXP conversigmad)
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
                                         as<double> (gammainit),
                                         as<int> (survvar),
                                         as<int> (conversigmad));
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

// [[Rcpp::export]]
Rcpp::List  Simdata_main(SEXP n, SEXP sim, SEXP nbreak, SEXP tL, SEXP tU, SEXP q_eta, 
                         SEXP j_max, SEXP p01, SEXP p02, SEXP t_max, SEXP distr, SEXP m_age,
                         SEXP std_age, SEXP k_max,
                         SEXP tbthetanew, SEXP sigmau_invnew,
                         SEXP tbeta0new, SEXP tbeta1new,
                         SEXP tsigmanew, SEXP tthetanew,
                         SEXP tsigmadnew, SEXP tetanew,
                         SEXP tgamma, SEXP lambda0)
{
  Rcpp::List result;
  try {
    
    result=Simdataspace::Simdata_cmain(as<int> (n), as<int> (sim), 
                                         as<int> (nbreak),
                                         as<double> (tL), as<double> (tU), 
                                         as<int> (q_eta), as<int> (j_max), 
                                         as<int>(p01), as<int>(p02),
                                         as<int>(t_max), as<int>(distr),
                                         as<double>(m_age), as<double>(std_age), 
                                         as<int>(k_max),  
                                         as<std::string> (tbthetanew),
                                         as<std::string> (sigmau_invnew),
                                         as<std::string>(tbeta0new),
                                         as<std::string>(tbeta1new),
                                         as<std::string>(tsigmanew),
                                         as<std::string>(tthetanew),
                                         as<std::string>(tsigmadnew),
                                         as<std::string>(tetanew),
                                         as<double> (tgamma),
                                         as<double> (lambda0)
                                         );
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
Rcpp::List  getfitted_main(SEXP tL, SEXP tU, SEXP nbreak, SEXP k_max, SEXP j_max, SEXP p01, SEXP p02, 
                           SEXP sigmau_invnew, SEXP thetanew, SEXP bthetanew, SEXP beta0new, SEXP beta1new, 
                           SEXP longvarnew)
{
  Rcpp::List result;
  try {
    
    result=getfittedspace::getfitted_cmain(as<double> (tL), as<double> (tU), 
                                       as<int> (nbreak), as<int> (k_max), 
                                       as<int>(j_max), as<int>(p01), as<int>(p02),
                                       as<std::string> (sigmau_invnew),
                                       as<std::string> (thetanew),
                                       as<std::string> (bthetanew),
                                       as<std::string> (beta0new),
                                       as<std::string> (beta1new),
                                       as<std::string> (longvarnew)
    );
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