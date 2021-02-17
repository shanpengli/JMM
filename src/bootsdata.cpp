#include "bootsdata.hpp"

namespace bootsdataspace {

Rcpp::List bootsdata_cmain(int n, int n1, double tL, double tU, int q_eta, 
                           int j_max, int p01, int p02, int t_max,
                           std::string ydatanew, std::string cdatanew,
                          std::string mdatanew, int nboots)

{
  int pt=p01+p02;
  
  /* allocate space for data */
  gsl_vector *M= gsl_vector_alloc(n);                               /*** # obs per subject ***/
  gsl_matrix *C = gsl_matrix_alloc(n,q_eta+2);                      /*** data for event times ***/
  gsl_matrix *Bio = gsl_matrix_alloc(n1,j_max+pt+1);           /*** data for biomarkers ***/ 
  /* read matrices */
  /* read Y matrix  */
  {
    FILE * f = fopen(ydatanew.c_str(), "r");
    
    if (f == NULL)
    {
      Rprintf("File %s does not exist.\n", ydatanew.c_str());
      return R_NilValue;
    }
    
    
    int nrows=0;
    // Extract characters from file and store in character c
    for (char c = fgetc(f); c != EOF; c = fgetc(f))
      if (c == '\n')  nrows = nrows + 1;
      nrows=nrows+1;
      if (n1 == nrows)
      {   rewind(f);
        gsl_matrix_fscanf(f, Bio);
        fclose(f);
      }
      else
      {
        Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                n1,
                ydatanew.c_str(),nrows);
        fclose(f);
        return R_NilValue;
      }
  }
  
  /*read C matrix*/
  {
    FILE * f = fopen(cdatanew.c_str(), "r");
    
    if (f == NULL)
    {
      Rprintf("File %s does not exist.\n", cdatanew.c_str());
      return R_NilValue;
    }
    
    
    int nrows=0;
    // Extract characters from file and store in character c
    for (char c = fgetc(f); c != EOF; c = fgetc(f))
      if (c == '\n')  nrows = nrows + 1;
      nrows=nrows+1;
      if (n == nrows)
      {   rewind(f);
        gsl_matrix_fscanf(f, C);
        fclose(f);
      }
      else
      {
        Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                n,
                cdatanew.c_str(),nrows);
        fclose(f);
        return R_NilValue;
      }
  }
  
  /* read M vector  */
  {
    FILE * f = fopen(mdatanew.c_str(), "r");
    
    if (f == NULL)
    {
      Rprintf("File %s does not exist.\n", mdatanew.c_str());
      return R_NilValue;
    }
    
    
    int nrows=0;
    // Extract characters from file and store in character c
    for (char c = fgetc(f); c != EOF; c = fgetc(f))
      if (c == '\n')  nrows = nrows + 1;
      nrows=nrows+1;
      if (n==nrows)
      {   rewind(f);
        gsl_vector_fscanf(f, M);
        fclose(f);
      }
      else
      {
        Rprintf("Input subjects is %d, but the number of rows in %s is %d",n,
                mdatanew.c_str(),nrows);
        fclose(f);
        return R_NilValue;
      }
  }
  
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  int i,j,t,k,ni,u,v,a,loc,p,q, point;
  
  for(i=0;i<C->size1;i++)
  {
    if(gsl_matrix_get(C,i,0)>tU) 
    {
      gsl_matrix_set(C,i,0,tU);	
      gsl_matrix_set(C,i,1,0);
    }
  }
  
  double temp;
  
  int aa[1],bb[n];
  
  for(i=0;i<n;i++)  bb[i]=i;
  
  gsl_vector *BM= gsl_vector_alloc(n);                               
  gsl_matrix *BC = gsl_matrix_alloc(n,q_eta+2);                     
  gsl_matrix *HBio = gsl_matrix_alloc(t_max*n,j_max+pt+1);    
  
  int boots=0;
  
  Rcpp::List ret;
  
  char namey[30];
  char namec[30];
  char namem[30];
  
    do 
    {
    
    gsl_vector_set_zero(BM);
    gsl_matrix_set_zero(BC);
    gsl_matrix_set_zero(HBio);
    
    loc=0;
    for(i=0;i<n;i++)
    {
      gsl_ran_choose(r,aa,1,bb,n,sizeof(int));
      
      ni=(int)gsl_vector_get(M,aa[0]);
      gsl_vector_set(BM,i,(double)ni);
      for(j=0;j<C->size2;j++)  gsl_matrix_set(BC,i,j,gsl_matrix_get(C,aa[0],j));
      
      point=0;
      for(t=0;t<aa[0];t++)  point+=(int)gsl_vector_get(M,t);
      
      for(t=0;t<ni;t++)
      {
        for(j=0;j<Bio->size2;j++)  gsl_matrix_set(HBio,loc+t,j,gsl_matrix_get(Bio,point+t,j));
        
      }
      
      loc+=ni;
    }
    

    gsl_matrix *BBio=gsl_matrix_alloc(loc,Bio->size2);
    
    for(i=0;i<loc;i++)  
    {
      for(j=0;j<Bio->size2;j++)  gsl_matrix_set(BBio,i,j,gsl_matrix_get(HBio,i,j));
    }
    
    /** output the estimates  ***/
    NumericMatrix Bio_matrix(BBio->size1, BBio->size2);
    NumericMatrix C_matrix(BC->size1, BC->size2);
    NumericVector m_vec(BM->size);
    
    for (i=0;i<BBio->size1;i++) 
    {
      for (j=0;j<BBio->size2;j++) 
      {
        Bio_matrix(i, j) = gsl_matrix_get(BBio, i, j);
      }
    }
    
    for (i=0;i<BC->size1;i++) 
    {
      for (j=0;j<BC->size2;j++) 
      {
        C_matrix(i, j) = gsl_matrix_get(BC, i, j);
      }
    }
    for (i=0;i<BM->size;i++) m_vec(i) = gsl_vector_get(BM, i);
    
    sprintf(namey, "ydata_%d", boots);
    sprintf(namec, "cdata_%d", boots);
    sprintf(namem, "mdata_%d", boots);
    
    ret[namey] = Bio_matrix;
    ret[namec] = C_matrix;
    ret[namem] = m_vec;
    
    gsl_matrix_free(BBio);
    
    boots+=1;
    
  }while(boots<nboots);
  
  gsl_matrix_free(HBio);
  gsl_matrix_free(BC);
  gsl_vector_free(BM);
  
  gsl_matrix_free(Bio);
  gsl_matrix_free(C);
  gsl_vector_free(M);
  
  return ret;
  
  }

    
    



}



