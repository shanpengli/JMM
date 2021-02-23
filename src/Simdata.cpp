#include "Simdata.hpp"

namespace Simdataspace {

  double Min(const double t1, const double t2)
  {
    if(t1<t2) return t1;
    else return t2;
  }

  int GetN(double t)
  {
    
    return (int)(t/0.5+1);
    
  }

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
      const int k_cubic)
  {
    double stepsize=0.001,temp,temp1,temp2,temp3,bmt,cmf;
    int nstep=(int)(tu-tl)/stepsize;
    int k,j,p;
    
    gsl_bspline_workspace *bw=gsl_bspline_alloc(k_cubic,tnbreak);
    gsl_bspline_knots_uniform(tl,tu,bw);
    
    gsl_vector *helpb=gsl_vector_alloc(B_spl->size);
    
    cmf=gsl_ran_flat(r,0,1);
    
    temp=0;
    k=0;
    do
    {
      bmt=k*stepsize+tl;
      
      temp1=0;
      for(j=0;j<eta->size;j++)
      {
        /*if(gsl_matrix_get(C,i,2+j)<=bmt)*/
        temp1+=gsl_vector_get(eta,j)*gsl_matrix_get(C,i,2+j);
      }
      
      gsl_bspline_eval(bmt,B_spl,bw); 
      
      
      gsl_vector_memcpy(helpb,B_spl);
      MulM(sigmau_inv,helpb,B_spl);
      
      temp2=0;
      for(j=0;j<B_spl->size;j++)
        temp2+=gsl_vector_get(B_spl,j)*gsl_vector_get(theta,j);
      
      temp1+=temp2*gamma;
      
      temp2=0;
      for(p=0;p<alphai->size;p++)
      {
        temp3=0;
        for(j=0;j<B_spl->size;j++)
          temp3+=gsl_vector_get(B_spl,j)*gsl_matrix_get(btheta,j,p);
        temp3*=gsl_vector_get(alphai,p);
        temp2+=temp3;
      }
      
      temp1+=temp2*gamma;
      
      temp+=lambda0*exp(temp1)*stepsize;
      k+=1;
      
      if(bmt>=tu) return tu;
      
    }while((1-exp(0-temp))<cmf);
    
    gsl_vector_free(helpb);
    gsl_bspline_free(bw);
    
    return bmt;
    
  }

  double GetMU(const gsl_vector *theta, 
               const gsl_matrix *btheta, 
               const gsl_vector *alphai, 
               gsl_vector *B_spl,
               const gsl_matrix *sigmau_inv)
  {
    gsl_vector *helpv=gsl_vector_alloc(btheta->size1);
    gsl_vector *helpb=gsl_vector_alloc(B_spl->size);
    MulM(btheta,alphai,helpv);
    
    gsl_vector_memcpy(helpb,B_spl);
    MulM(sigmau_inv,helpb,B_spl);
    
    double temp=MulVV(B_spl,theta)+MulVV(B_spl,helpv);
    
    gsl_vector_free(helpv);
    gsl_vector_free(helpb);
    
    return temp;
    
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
  
  
  void MulV(const gsl_vector *Z,gsl_matrix *ZZ)
  {
    int p = Z->size;
    int i,j;
    
    for(i=0;i<p;i++)
    {
      for(j=0;j<p;j++) gsl_matrix_set(ZZ,i,j,gsl_vector_get(Z,i)*gsl_vector_get(Z,j));
    }
  }
  
  
  double MulVV(const gsl_vector *Z,const gsl_vector *beta)
  {
    int p=Z->size;
    int i;
    double temp=0;
    
    for(i=0;i<p;i++)  temp+=gsl_vector_get(Z,i)*gsl_vector_get(beta,i);
    
    return (temp);
  }
  
  
  void MulMM(const gsl_matrix *A,const gsl_matrix *B,gsl_matrix *AB)
  {
    int p=A->size1;
    int q=A->size2;
    int k=B->size2;
    
    int i,j,t;
    double temp;
    
    for(i=0;i<p;i++)  
    {
      for(j=0;j<k;j++)
      {
        temp=0;
        for(t=0;t<q;t++)  temp+=gsl_matrix_get(A,i,t)*gsl_matrix_get(B,t,j);
        gsl_matrix_set(AB,i,j,temp);
      }
    }
    
  }

  Rcpp::List Simdata_cmain(int n, int sim, int nbreak, double tL, double tU, int q_eta, 
                           int j_max, int p01, int p02, int t_max, int distr, double m_age,
                           double std_age, int k_max,
                           std::string tbthetanew, std::string sigmau_invnew,
                           std::string tbeta0new, std::string tbeta1new,
                           std::string tsigmanew, std::string tthetanew,
                           std::string tsigmadnew, std::string tetanew,
                           double tgamma, double lambda0
                           )
  {
    int k_cubic = 4;
    int n_total;
    int pt = p01+p02;
    int tq_b=k_cubic+nbreak-2;
    gsl_vector *pmax = gsl_vector_calloc(2);
    gsl_vector_set(pmax, 0, p01);
    gsl_vector_set(pmax, 1, p02);
    int p_max = gsl_vector_max(pmax);
    
    gsl_matrix *tsigmau_inv=gsl_matrix_alloc(nbreak+k_cubic-2,nbreak+k_cubic-2);
    
    /* read sigmau_inv matrix*/
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
          gsl_matrix_fscanf(f, tsigmau_inv);
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
    
    /* allocate space for data */
    gsl_vector *M= gsl_vector_alloc(n);                               /*** # obs per subject ***/
    gsl_matrix *C = gsl_matrix_alloc(n,q_eta+2);                /*** data for event times ***/
    gsl_matrix *Bio = gsl_matrix_alloc(n*t_max,j_max+pt+1);           /*** data for biomarkers ***/ 
    
    gsl_matrix_set_zero(Bio);
    gsl_matrix_set_zero(C);
    gsl_vector_set_zero(M);
    
    gsl_vector *pdim=gsl_vector_alloc(j_max);
    gsl_vector_set(pdim,0,p01);
    gsl_vector_set(pdim,1,p02);
    
    /* allocate space for true parameters */
    gsl_matrix *tbeta0 = gsl_matrix_alloc(j_max,p_max);
    gsl_vector *tbeta1 = gsl_vector_alloc(j_max),
               *tsigma = gsl_vector_alloc(j_max),
               *ttheta = gsl_vector_alloc(tq_b),
               *tsigmad = gsl_vector_alloc(k_max),
               *teta = gsl_vector_alloc(q_eta);
      
    gsl_matrix *tbtheta = gsl_matrix_alloc(tq_b,k_max);
    
    /* read tbtheta matrix*/
    {
      FILE * f = fopen(tbthetanew.c_str(), "r");
      
      if (f == NULL)
      {
        Rprintf("File %s does not exist.\n", tbthetanew.c_str());
        return R_NilValue;
      }
      
      
      int nrows=0;
      // Extract characters from file and store in character c
      for (char c = fgetc(f); c != EOF; c = fgetc(f))
        if (c == '\n')  nrows = nrows + 1;
        nrows=nrows+1;
        if (tq_b == nrows)
        {   rewind(f);
          gsl_matrix_fscanf(f, tbtheta);
          fclose(f);
        }
        else
        {
          Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                  tq_b,
                  tbthetanew.c_str(),nrows);
          fclose(f);
          return R_NilValue;
        }
    }
    
    /* read tbeta0 matrix  */
    {
      FILE * f = fopen(tbeta0new.c_str(), "r");
      
      if (f == NULL)
      {
        Rprintf("File %s does not exist.\n", tbeta0new.c_str());
        return R_NilValue;
      }
      
      
      int nrows=0;
      // Extract characters from file and store in character c
      for (char c = fgetc(f); c != EOF; c = fgetc(f))
        if (c == '\n')  nrows = nrows + 1;
        nrows=nrows+1;
        if (j_max == nrows)
        {   rewind(f);
          gsl_matrix_fscanf(f, tbeta0);
          fclose(f);
        }
        else
        {
          Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                  j_max,
                  tbeta0new.c_str(),nrows);
          fclose(f);
          return R_NilValue;
        }
    }
    
    /* read tbeta1 vector  */
    {
      FILE * f = fopen(tbeta1new.c_str(), "r");
      
      if (f == NULL)
      {
        Rprintf("File %s does not exist.\n", tbeta1new.c_str());
        return R_NilValue;
      }
      
      
      int nrows=0;
      // Extract characters from file and store in character c
      for (char c = fgetc(f); c != EOF; c = fgetc(f))
        if (c == '\n')  nrows = nrows + 1;
        nrows=nrows+1;
        if (j_max==nrows)
        {   rewind(f);
          gsl_vector_fscanf(f, tbeta1);
          fclose(f);
        }
        else
        {
          Rprintf("Input subjects is %d, but the number of rows in %s is %d",j_max,
                  tbeta1new.c_str(),nrows);
          fclose(f);
          return R_NilValue;
        }
    }
    
    /* read tsigma vector  */
    {
      FILE * f = fopen(tsigmanew.c_str(), "r");
      
      if (f == NULL)
      {
        Rprintf("File %s does not exist.\n", tsigmanew.c_str());
        return R_NilValue;
      }
      
      
      int nrows=0;
      // Extract characters from file and store in character c
      for (char c = fgetc(f); c != EOF; c = fgetc(f))
        if (c == '\n')  nrows = nrows + 1;
        nrows=nrows+1;
        if (j_max==nrows)
        {   rewind(f);
          gsl_vector_fscanf(f, tsigma);
          fclose(f);
        }
        else
        {
          Rprintf("Input subjects is %d, but the number of rows in %s is %d",j_max,
                  tsigmanew.c_str(),nrows);
          fclose(f);
          return R_NilValue;
        }
    }
    
    /* read theta vector  */
    {
      FILE * f = fopen(tthetanew.c_str(), "r");
      
      if (f == NULL)
      {
        Rprintf("File %s does not exist.\n", tthetanew.c_str());
        return R_NilValue;
      }
      
      
      int nrows=0;
      // Extract characters from file and store in character c
      for (char c = fgetc(f); c != EOF; c = fgetc(f))
        if (c == '\n')  nrows = nrows + 1;
        nrows=nrows+1;
        if (tq_b==nrows)
        {   rewind(f);
          gsl_vector_fscanf(f, ttheta);
          fclose(f);
        }
        else
        {
          Rprintf("Input subjects is %d, but the number of rows in %s is %d",tq_b,
                  tthetanew.c_str(),nrows);
          fclose(f);
          return R_NilValue;
        }
    }
    
    /* read tsigmad vector  */
    {
      FILE * f = fopen(tsigmadnew.c_str(), "r");
      
      if (f == NULL)
      {
        Rprintf("File %s does not exist.\n", tsigmadnew.c_str());
        return R_NilValue;
      }
      
      
      int nrows=0;
      // Extract characters from file and store in character c
      for (char c = fgetc(f); c != EOF; c = fgetc(f))
        if (c == '\n')  nrows = nrows + 1;
        nrows=nrows+1;
        if (k_max==nrows)
        {   rewind(f);
          gsl_vector_fscanf(f, tsigmad);
          fclose(f);
        }
        else
        {
          Rprintf("Input subjects is %d, but the number of rows in %s is %d",k_max,
                  tsigmadnew.c_str(),nrows);
          fclose(f);
          return R_NilValue;
        }
    }
    
    /* read teta vector  */
    {
      FILE * f = fopen(tetanew.c_str(), "r");
      
      if (f == NULL)
      {
        Rprintf("File %s does not exist.\n", tetanew.c_str());
        return R_NilValue;
      }
      
      
      int nrows=0;
      // Extract characters from file and store in character c
      for (char c = fgetc(f); c != EOF; c = fgetc(f))
        if (c == '\n')  nrows = nrows + 1;
        nrows=nrows+1;
        if (q_eta==nrows)
        {   rewind(f);
          gsl_vector_fscanf(f, teta);
          fclose(f);
        }
        else
        {
          Rprintf("Input subjects is %d, but the number of rows in %s is %d",q_eta,
                  tetanew.c_str(),nrows);
          fclose(f);
          return R_NilValue;
        }
    }
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    
    gsl_bspline_workspace *bw=gsl_bspline_alloc(k_cubic,nbreak);
    gsl_bspline_knots_uniform(tL,tU,bw);
    gsl_vector *B_spl=gsl_vector_alloc(nbreak+k_cubic-2);
    
    
    gsl_bspline_workspace *tbw=gsl_bspline_alloc(k_cubic,nbreak);
    gsl_bspline_knots_uniform(tL,tU,tbw);
    gsl_vector *tB_spl=gsl_vector_alloc(nbreak+k_cubic-2);
    
    int point,i,j,t,k,ni,iter,status,simn,loc,p,q;
    double crate,cratet,rate1,rate1t,temp,mu,/*dad,ar,age70,male,*/age,t1,censor;
    double max_censor=tU; /*probage=0.3,probmale=0.6,*/
    gsl_vector *alphai=gsl_vector_alloc(k_max);
    
    Rcpp::List ret;
    
    char namey[30];
    char namec[30];
    char namem[30];
    
    cratet=0;
    rate1t=0;
    
    for (simn=0;simn<sim;simn++) {
      gsl_vector_set_zero(M);
      gsl_matrix_set_zero(C); 
      gsl_matrix_set_zero(Bio);
      
      crate=0;
      rate1=0;
      point=0;
      
      for(i=0;i<n;i++)
      {
        
        if (distr == 1) {
          for(k=0;k<k_max;k++)    
          {
            gsl_vector_set(alphai,k,gsl_ran_gaussian(r,sqrt(gsl_vector_get(tsigmad,k))));
          }   
        } else {
          gsl_vector_set(alphai,0,gsl_ran_gamma(r,3,2)-6);
          gsl_vector_set(alphai,1,gsl_ran_gamma(r,2,2)-4);
        }

        age=gsl_ran_gaussian(r,std_age)+m_age;
        
        censor=gsl_ran_exponential(r,20);
        
        gsl_matrix_set(C,i,2,age);
        
        
        t1=GETT(r,i,tL,tU,nbreak,tB_spl,ttheta,tbtheta,alphai,lambda0,tgamma,teta,C,tsigmau_inv,k_cubic);
        
        
        if(t1<Min(censor,max_censor))
        {
          rate1+=1;
          gsl_matrix_set(C,i,0,t1);
          gsl_matrix_set(C,i,1,1);
          ni=GetN(t1);          
        }
        
        if(Min(censor,max_censor)<=t1)
        {
          crate+=1;
          gsl_matrix_set(C,i,0,Min(max_censor,censor));
          gsl_matrix_set(C,i,1,0);
          ni=GetN(Min(max_censor,censor));
        }  
        
        /*printf("ni=%d\n",ni);*/
        
        for(j=point;j<point+ni;j++)
        {
          
          gsl_matrix_set(Bio,j,2,age);
          /*
           gsl_matrix_set(Bio,j,3,male);
           */
          gsl_matrix_set(Bio,j,p01+3,age);
          /*
           gsl_matrix_set(Bio,j,p01+4,male);
           */             
          t1=(j-point)+gsl_ran_flat(r,-0.3,0.3);
          if(t1<tL) t1=tL;
          else if(t1>tU) t1=tU;
          gsl_bspline_eval(t1,tB_spl,tbw); 
          mu=GetMU(ttheta,tbtheta,alphai,tB_spl,tsigmau_inv);
          
          gsl_matrix_set(Bio,j,0,t1);
          
          temp=0;
          for(k=0;k<p01;k++)  temp+=gsl_matrix_get(Bio,j,k+2)*gsl_matrix_get(tbeta0,0,k);
          temp+=mu*gsl_vector_get(tbeta1,0)+gsl_ran_gaussian(r,sqrt(gsl_vector_get(tsigma,0)));
          gsl_matrix_set(Bio,j,1,temp);
          
          temp=0;
          for(k=0;k<p02;k++)  temp+=gsl_matrix_get(Bio,j,p01+1+k+2)*gsl_matrix_get(tbeta0,1,k);
          temp+=mu*gsl_vector_get(tbeta1,1)+gsl_ran_gaussian(r,sqrt(gsl_vector_get(tsigma,1)));
          gsl_matrix_set(Bio,j,p01+2,temp);
          
          
        }
        
        point+=ni;
        gsl_vector_set(M,i,(double)ni);
        
      }
      
      /** output the estimates  ***/
      NumericMatrix Bio_matrix(point, j_max+pt+1);
      NumericMatrix C_matrix(C->size1, C->size2);
      NumericVector m_vec(M->size);
      
      for (i=0;i<point;i++) 
      {
        for (j=0;j<(j_max+pt+1);j++) 
        {
          Bio_matrix(i, j) = gsl_matrix_get(Bio, i, j);
        }
      }
      
      for (i=0;i<C->size1;i++) 
      {
        for (j=0;j<C->size2;j++) 
        {
          C_matrix(i, j) = gsl_matrix_get(C, i, j);
        }
      }
      for (i=0;i<M->size;i++) m_vec(i) = gsl_vector_get(M, i);
      
      sprintf(namey, "ydata_%d", simn);
      sprintf(namec, "cdata_%d", simn);
      sprintf(namem, "mdata_%d", simn);
      
      ret[namey] = Bio_matrix;
      ret[namec] = C_matrix;
      ret[namem] = m_vec;
      
      cratet+=crate;
      rate1t+=rate1;
      
      }
    
    Rprintf("Average censoring rate=%f\n", cratet/((double)n*sim));
    Rprintf("Average event rate=%f\n", rate1t/((double)n*sim));
    
    gsl_matrix_free(Bio);
    gsl_matrix_free(C);
    gsl_vector_free(M);
    
    gsl_matrix_free(tbeta0);
    gsl_vector_free(tbeta1);
    gsl_vector_free(tsigma);
    gsl_vector_free(ttheta);
    gsl_vector_free(tsigmad);
    gsl_vector_free(teta);
    gsl_matrix_free(tbtheta);
    
    gsl_vector_free(B_spl);
    gsl_vector_free(pdim);
    gsl_vector_free(pmax);
    gsl_vector_free(alphai);
    gsl_bspline_free(bw);
    gsl_matrix_free(tsigmau_inv);
    
    gsl_vector_free(tB_spl); 
    gsl_bspline_free(tbw);
    
    
    return(ret);
    
  }


}