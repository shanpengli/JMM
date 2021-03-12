#include "getfitted.hpp"

namespace getfittedspace {

    double MulVV(const gsl_vector *Z,const gsl_vector *beta)
    {
      int p=Z->size;
      int i;
      double temp=0;
      
      for(i=0;i<p;i++)  temp+=gsl_vector_get(Z,i)*gsl_vector_get(beta,i);
      
      return (temp);
    }
    
    
    void MulM(const gsl_matrix *XX, const gsl_vector *X, gsl_vector *beta)
    {
      int p = XX->size1;
      int q = XX->size2;
      
      int i,j;
      double temp;
      
      for(i=0;i<p;i++)
      {
        temp=0;
        for(j=0;j<q;j++)  temp+=gsl_matrix_get(XX,i,j)*gsl_vector_get(X,j);
        gsl_vector_set(beta,i,temp);
      }
      
    }

    Rcpp::List getfitted_cmain(double tL, double tU, int nbreak, int k_max, int j_max,
                               int p01, int p02,
                               std::string sigmau_invnew, 
                               std::string thetanew, 
                               std::string bthetanew,
                               std::string beta0new,
                               std::string beta1new,
                               std::string longvarnew)
    {
      int k_cubic = 4;
      int q_b=k_cubic+nbreak-2;
      
      /* allocate space for estimated parameters */
      /*Define a vector of # of covariates for all biomarkers; 
       the current version can only allow two biomarkers*/
      gsl_vector *pmax = gsl_vector_calloc(2);
      gsl_vector_set(pmax, 0, p01);
      gsl_vector_set(pmax, 1, p02);
      int p_max = gsl_vector_max(pmax);
      gsl_matrix *beta0 = gsl_matrix_alloc(j_max,p_max);
      gsl_vector *beta1 = gsl_vector_alloc(j_max);
      
      /* read sigmau_inv matrix*/
      gsl_matrix *sigmau_inv=gsl_matrix_alloc(nbreak+k_cubic-2,nbreak+k_cubic-2);
      {
        FILE * f = fopen(sigmau_invnew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", sigmau_invnew.c_str());
          return R_NilValue;
        }
        
        
        int nrows=0;
        // Extract characters from file and store in character c
        for (char c = fgetc(f); c != EOF; c = fgetc(f))
          if (c == '\n')  nrows = nrows + 1;
          nrows=nrows+1;
          if (nbreak+k_cubic-2 == nrows)
          {   rewind(f);
            gsl_matrix_fscanf(f, sigmau_inv);
            fclose(f);
          }
          else
          {
            Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                    nbreak+k_cubic-2,
                    sigmau_invnew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
        
      /* read theta vector  */
      gsl_vector *theta = gsl_vector_alloc(q_b);
      {
        FILE * f = fopen(thetanew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", thetanew.c_str());
          return R_NilValue;
        }
        
        
        int nrows=0;
        // Extract characters from file and store in character c
        for (char c = fgetc(f); c != EOF; c = fgetc(f))
          if (c == '\n')  nrows = nrows + 1;
          nrows=nrows+1;
          if (q_b==nrows)
          {   rewind(f);
            gsl_vector_fscanf(f, theta);
            fclose(f);
          }
          else
          {
            Rprintf("Input subjects is %d, but the number of rows in %s is %d",q_b,
                    thetanew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      gsl_matrix *btheta = gsl_matrix_alloc(q_b,k_max);
      {
        FILE * f = fopen(bthetanew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", bthetanew.c_str());
          return R_NilValue;
        }
        
        
        int nrows=0;
        // Extract characters from file and store in character c
        for (char c = fgetc(f); c != EOF; c = fgetc(f))
          if (c == '\n')  nrows = nrows + 1;
          nrows=nrows+1;
          if (q_b == nrows)
          {   rewind(f);
            gsl_matrix_fscanf(f, btheta);
            fclose(f);
          }
          else
          {
            Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                    q_b,
                    bthetanew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      /* read beta0 matrix  */
      {
        FILE * f = fopen(beta0new.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", beta0new.c_str());
          return R_NilValue;
        }
        
        
        int nrows=0;
        // Extract characters from file and store in character c
        for (char c = fgetc(f); c != EOF; c = fgetc(f))
          if (c == '\n')  nrows = nrows + 1;
          nrows=nrows+1;
          if (j_max == nrows)
          {   rewind(f);
            gsl_matrix_fscanf(f, beta0);
            fclose(f);
          }
          else
          {
            Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                    j_max,
                    beta0new.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      /* read beta1 vector  */
      {
        FILE * f = fopen(beta1new.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", beta1new.c_str());
          return R_NilValue;
        }
        
        
        int nrows=0;
        // Extract characters from file and store in character c
        for (char c = fgetc(f); c != EOF; c = fgetc(f))
          if (c == '\n')  nrows = nrows + 1;
          nrows=nrows+1;
          if (j_max==nrows)
          {   rewind(f);
            gsl_vector_fscanf(f, beta1);
            fclose(f);
          }
          else
          {
            Rprintf("Input subjects is %d, but the number of rows in %s is %d",j_max,
                    beta1new.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      /* read longvar matrix*/
      gsl_matrix *longvarmatrix=gsl_matrix_alloc(j_max,p_max);
      {
        FILE * f = fopen(longvarnew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", longvarnew.c_str());
          return R_NilValue;
        }
        
        
        int nrows=0;
        // Extract characters from file and store in character c
        for (char c = fgetc(f); c != EOF; c = fgetc(f))
          if (c == '\n')  nrows = nrows + 1;
          nrows=nrows+1;
          if (j_max == nrows)
          {   rewind(f);
            gsl_matrix_fscanf(f, longvarmatrix);
            fclose(f);
          }
          else
          {
            Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                    j_max,
                    longvarnew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      double stepsize=0.01;
      int nstep=(tU-tL)/stepsize+1;
      double timex;
      
      gsl_bspline_workspace *bw=gsl_bspline_alloc(k_cubic,nbreak);
      gsl_bspline_knots_uniform(tL,tU,bw);
      
      
      gsl_vector *B_spl=gsl_vector_alloc(q_b);
      gsl_vector *helpb=gsl_vector_alloc(q_b);
      gsl_vector *htheta = gsl_vector_alloc(q_b);
      
      /** output the estimates  ***/
      NumericMatrix fittedvalue_matrix(nstep, 3 + k_max);
      
      int i,k,j;
      double sum=0;
      for(i=0;i<nstep;i++)
      {
        timex=(double)i*stepsize+tL;
        gsl_bspline_eval(timex,B_spl,bw);
        gsl_vector_memcpy(helpb,B_spl);
        MulM(sigmau_inv,helpb,B_spl);
        
        fittedvalue_matrix(i, 0) = timex;
        for (j=0;j<p01;j++) sum = sum + gsl_matrix_get(beta0, 0, j)*gsl_matrix_get(longvarmatrix, 0, j);
        fittedvalue_matrix(i, 1) = sum+MulVV(theta,B_spl)*gsl_vector_get(beta1, 0);
        sum = 0;
        for (j=0;j<p02;j++) sum = sum + gsl_matrix_get(beta0, 1, j)*gsl_matrix_get(longvarmatrix, 1, j);
        fittedvalue_matrix(i, 2) = sum+gsl_vector_get(beta1, 1)*MulVV(theta,B_spl);
        sum = 0;
        
        for(k=0;k<k_max;k++)
        {
          gsl_matrix_get_col(htheta,btheta,k);
          fittedvalue_matrix(i, 3+k) = MulVV(htheta,B_spl);
        }
        
      }
      
      Rcpp::List ret;
      ret["fittedvalue"] = fittedvalue_matrix;
      
      gsl_matrix_free(sigmau_inv);
      gsl_matrix_free(btheta);
      gsl_vector_free(theta);
      gsl_vector_free(htheta);
      gsl_vector_free(B_spl);
      gsl_vector_free(helpb);
      gsl_bspline_free(bw);
      gsl_matrix_free(beta0);
      gsl_vector_free(beta1);
      return ret;
      
      }


}