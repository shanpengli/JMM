//
//  jmspline.cpp
//  FastJM
//
//  Created by Shanpeng Li on 6/16/20.
//

#include "jmspline.hpp"

namespace jmsplinespace {

  int gsl_linalg_cholesky_decompn (gsl_matrix * A)
  {
    const size_t M = A->size1;
    const size_t N = A->size2;
    
    if (M != N)
    {
      return 100;
    }
    else
    {
      size_t i,j,k;
      int status = 0;
      
      /* Do the first 2 rows explicitly.  It is simple, and faster.  And
       * one can return if the matrix has only 1 or 2 rows.  
       */
      
      double A_00 = gsl_matrix_get (A, 0, 0);
      
      double L_00 = sqrt(A_00);
      
      if (A_00 <= 0)
      {
        return 100 ;
      }
      
      gsl_matrix_set (A, 0, 0, L_00);
      
      if (M > 1)
      {
        double A_10 = gsl_matrix_get (A, 1, 0);
        double A_11 = gsl_matrix_get (A, 1, 1);
        
        double L_10 = A_10 / L_00;
        double diag = A_11 - L_10 * L_10;
        double L_11 = sqrt(diag);
        
        if (diag <= 0)
        {
          return 100;
        }
        
        gsl_matrix_set (A, 1, 0, L_10);        
        gsl_matrix_set (A, 1, 1, L_11);
      }
      
      for (k = 2; k < M; k++)
      {
        double A_kk = gsl_matrix_get (A, k, k);
        
        for (i = 0; i < k; i++)
        {
          double sum = 0;
          
          double A_ki = gsl_matrix_get (A, k, i);
          double A_ii = gsl_matrix_get (A, i, i);
          
          gsl_vector_view ci = gsl_matrix_row (A, i);
          gsl_vector_view ck = gsl_matrix_row (A, k);
          
          if (i > 0) {
            gsl_vector_view di = gsl_vector_subvector(&ci.vector, 0, i);
            gsl_vector_view dk = gsl_vector_subvector(&ck.vector, 0, i);
            
            gsl_blas_ddot (&di.vector, &dk.vector, &sum);
          }
          
          A_ki = (A_ki - sum) / A_ii;
          gsl_matrix_set (A, k, i, A_ki);
        } 
        
        {
          gsl_vector_view ck = gsl_matrix_row (A, k);
          gsl_vector_view dk = gsl_vector_subvector (&ck.vector, 0, k);
          
          double sum = gsl_blas_dnrm2 (&dk.vector);
          double diag = A_kk - sum * sum;
          
          double L_kk = sqrt(diag);
          
          if (diag <= 0)
          {
            return 100;
          }
          
          gsl_matrix_set (A, k, k, L_kk);
        }
      }
      
      /* Now copy the transposed lower triangle to the upper triangle,
       * the diagonal is common.  
       */
      
      for (i = 1; i < M; i++)
      {
        for (j = 0; j < i; j++)
        {
          double A_ij = gsl_matrix_get (A, i, j);
          gsl_matrix_set (A, j, i, A_ij);
        }
      } 
      
      
    }
    
    return 0; 
  }
  
  
  int gsl_linalg_cholesky_solven (const gsl_matrix * LLT,
                                  const gsl_vector * b,
                                  gsl_vector * x)
  {
    if (LLT->size1 != LLT->size2)
    {
      return 100;
    }
    else if (LLT->size1 != b->size)
    {
      return 100;
    }
    else if (LLT->size2 != x->size)
    {
      return 100;
    }
    else
    {
      /* Copy x <- b */
      
      gsl_vector_memcpy (x, b);
      
      /* Solve for c using forward-substitution, L c = b */
      
      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);
      
      /* Perform back-substitution, U x = c */
      
      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);
      
      
      return 0;
    }
  }
  
  int gsl_linalg_cholesky_svx (const gsl_matrix * LLT,
                               gsl_vector * x)
  {
    if (LLT->size1 != LLT->size2)
    {
      return 100;
    }
    else if (LLT->size2 != x->size)
    {
      return 100;
    }
    else
    {
      /* Solve for c using forward-substitution, L c = b */
      
      gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, LLT, x);
      
      /* Perform back-substitution, U x = c */
      
      gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, LLT, x);
      
      return 0;
    }
  }

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
      const std::vector<double> ws,
      const int survvar
      )
    {
      int p_max=beta0->size2;
      int j_max=beta0->size1;
      int q_b=theta->size;
      int k_max=sigmad->size;
      int q_eta=eta->size;
      int a=H01->size2;
      int n=M->size;
      double nt;
      
      int i,j,t,k,u,m,p,q,point,loc;
      int status;
      double temp1,temp2,temp3,temp4,temp5,temp6;
      
      
      double pgamma=*gamma;
      
      gsl_matrix *FUNA=gsl_matrix_alloc(n,k_max);
      gsl_matrix *FUNA2=gsl_matrix_alloc(n,((k_max+1)*k_max)/2);
      gsl_matrix *FUNE=gsl_matrix_alloc(n,a);
      gsl_matrix *FUNAE=gsl_matrix_alloc(n,a);
      gsl_matrix *FUNA2E=gsl_matrix_alloc(n,a);
      gsl_matrix *FUNAEV=gsl_matrix_alloc(n,a*k_max);
      gsl_matrix *FUNA2EV=gsl_matrix_alloc(n,a*k_max);
      
      status=GetE(FUNA,FUNA2,FUNE,FUNAE,FUNA2E,FUNAEV,FUNA2EV,beta0,beta1,sigma,theta,sigmad,eta,btheta,H01,pgamma,pdim,Bio,C,M,bw,sigmau_inv, quadpoint, xs, ws);                       
      
      
      gsl_vector *B_spl=gsl_vector_alloc(q_b);
      gsl_vector *helpb=gsl_vector_alloc(q_b);
      
      nt=0;
      for(i=0;i<n;i++)
        nt+=gsl_vector_get(M,i);
      
      
      /* update sigmad */
      
      for(i=0;i<k_max;i++)
      {
        temp1=0;
        for(j=0;j<n;j++) temp1+=gsl_matrix_get(FUNA2,j,i);
        temp1=temp1/(double)n;
        gsl_vector_set(sigmad,i,temp1);
      }
      
      
      
      /* update sigma */
      
      gsl_vector_set_zero(sigma);
      point=0;
      for(i=0;i<n;i++)
      {
        u=(int)gsl_vector_get(M,i);
        for(t=0;t<u;t++) 
        {
          gsl_bspline_eval(gsl_matrix_get(Bio,point+t,0),B_spl,bw);
          gsl_vector_memcpy(helpb,B_spl);
          MulM(sigmau_inv,helpb,B_spl);
          
          for(j=0;j<j_max;j++)
          {
            loc=0;
            for(m=0;m<j;m++) loc+=(int)gsl_vector_get(pdim,m);
            
            temp1=gsl_matrix_get(Bio,point+t,loc+j+1);
            
            for(k=0;k<(int)gsl_vector_get(pdim,j);k++)
              temp1-=gsl_matrix_get(Bio,point+t,loc+j+2+k)*gsl_matrix_get(beta0,j,k);
            
            temp2=0;
            for(k=0;k<q_b;k++)
              temp2+=gsl_vector_get(B_spl,k)*gsl_vector_get(theta,k);
            
            temp1=temp1-temp2*gsl_vector_get(beta1,j);
            
            temp2=0;
            for(m=0;m<k_max;m++)
            {
              temp3=0;
              for(k=0;k<q_b;k++)
                temp3+=gsl_vector_get(B_spl,k)*gsl_matrix_get(btheta,k,m);
              temp3*=gsl_matrix_get(FUNA,i,m);
              temp2+=temp3;
            }
            temp2*=gsl_vector_get(beta1,j);
            
            temp3=temp1*temp1-2*temp1*temp2;
            
            temp1=0;
            for(m=0;m<k_max;m++)
            {
              temp2=0;
              for(k=0;k<q_b;k++)
                temp2+=gsl_vector_get(B_spl,k)*gsl_matrix_get(btheta,k,m);
              
              temp1+=gsl_pow_2(temp2*gsl_vector_get(beta1,j))*gsl_matrix_get(FUNA2,i,m);
              
              for(p=0;p<m;p++)
              {
                temp4=0;
                for(k=0;k<q_b;k++)
                  temp4+=gsl_vector_get(B_spl,k)*gsl_matrix_get(btheta,k,p);
                temp1+=2*gsl_pow_2(gsl_vector_get(beta1,j))*temp2*temp4*gsl_matrix_get(FUNA2,i,k_max+m*(m-1)/2+p);
              }
            }
            
            gsl_vector_set(sigma,j,gsl_vector_get(sigma,j)+temp3+temp1);
          }
        }
        
        point+=u;
      }
      
      gsl_vector_scale(sigma,1/nt);
      
      
      
      /* update beta0 */
      
      gsl_matrix_set_zero(beta0);
      
      for(j=0;j<j_max;j++)
      {
        p=(int)gsl_vector_get(pdim,j);
        gsl_vector *X=gsl_vector_alloc(p);
        gsl_matrix *XX=gsl_matrix_alloc(p,p);
        gsl_vector *SX=gsl_vector_alloc(p);
        gsl_matrix *SXX=gsl_matrix_alloc(p,p);
        
        gsl_vector_set_zero(SX);
        gsl_matrix_set_zero(SXX);
        
        loc=0;
        for(m=0;m<j;m++) loc+=(int)gsl_vector_get(pdim,m);
        
        point=0;
        for(i=0;i<n;i++)
        {
          u=(int)gsl_vector_get(M,i);
          for(t=0;t<u;t++) 
          {
            gsl_bspline_eval(gsl_matrix_get(Bio,point+t,0),B_spl,bw); 
            gsl_vector_memcpy(helpb,B_spl);
            MulM(sigmau_inv,helpb,B_spl);
            
            temp1=gsl_matrix_get(Bio,point+t,loc+j+1);
            temp2=0;
            for(k=0;k<q_b;k++)
              temp2+=gsl_vector_get(B_spl,k)*gsl_vector_get(theta,k);
            temp1-=temp2*gsl_vector_get(beta1,j);
            
            temp2=0;
            for(m=0;m<k_max;m++)
            {
              temp3=0;
              for(k=0;k<q_b;k++)
                temp3+=gsl_vector_get(B_spl,k)*gsl_matrix_get(btheta,k,m);
              
              temp2+=temp3*gsl_vector_get(beta1,j)*gsl_matrix_get(FUNA,i,m);
              
            }
            
            temp1=temp1-temp2;
            
            for(k=0;k<p;k++) 
              gsl_vector_set(X,k,gsl_matrix_get(Bio,point+t,loc+j+2+k));
            
            MulV(X,XX);
            gsl_vector_scale(X,temp1);
            
            gsl_vector_add(SX,X);
            gsl_matrix_add(SXX,XX);
          }
          point+=u;
        }
        
        inv_matrix(SXX);
        if(status==100) return status;
        
        MulM(SXX,SX,X);
        
        for(k=0;k<p;k++) gsl_matrix_set(beta0,j,k,gsl_vector_get(X,k));
        
        gsl_vector_free(X);
        gsl_vector_free(SX);
        gsl_matrix_free(XX);
        gsl_matrix_free(SXX);
        
      }
      
      
      
      /* update beta1 */
      
      for(j=0;j<j_max;j++)
      {
        if(j==0) gsl_vector_set(beta1,j,1);
        else
        {
          p=(int)gsl_vector_get(pdim,j);
          
          loc=0;
          for(m=0;m<j;m++) loc+=(int)gsl_vector_get(pdim,m);
          
          temp1=0;
          temp2=0;
          
          point=0;
          for(i=0;i<n;i++)
          {
            u=(int)gsl_vector_get(M,i);
            for(t=0;t<u;t++) 
            {
              gsl_bspline_eval(gsl_matrix_get(Bio,point+t,0),B_spl,bw); 
              gsl_vector_memcpy(helpb,B_spl);
              MulM(sigmau_inv,helpb,B_spl);
              
              temp3=gsl_matrix_get(Bio,point+t,loc+j+1);
              
              for(k=0;k<p;k++)
                temp3-=gsl_matrix_get(Bio,point+t,loc+j+2+k)*gsl_matrix_get(beta0,j,k);
              
              temp4=0;
              for(k=0;k<q_b;k++)
                temp4+=gsl_vector_get(B_spl,k)*gsl_vector_get(theta,k);
              
              temp6=0;
              for(m=0;m<k_max;m++)
              {
                temp5=0;
                for(k=0;k<q_b;k++)
                  temp5+=gsl_vector_get(B_spl,k)*gsl_matrix_get(btheta,k,m);
                
                temp6+=temp5*gsl_matrix_get(FUNA,i,m);
              }
              
              temp1+=temp3*(temp4+temp6);
              
              
              temp3=0;
              for(m=0;m<k_max;m++)
              {
                temp5=0;
                for(k=0;k<q_b;k++)
                  temp5+=gsl_vector_get(B_spl,k)*gsl_matrix_get(btheta,k,m);
                
                temp3+=gsl_pow_2(temp5)*gsl_matrix_get(FUNA2,i,m);
                
                for(k=0;k<m;k++)
                {
                  temp6=0;
                  for(q=0;q<q_b;q++)
                    temp6+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,k);
                  temp3+=2*temp5*temp6*gsl_matrix_get(FUNA2,i,k_max+m*(m-1)/2+k);
                }
              }
              
              
              temp5=0;
              for(m=0;m<k_max;m++)
              {
                temp6=0;
                for(k=0;k<q_b;k++)
                  temp6+=gsl_vector_get(B_spl,k)*gsl_matrix_get(btheta,k,m);
                
                temp5+=temp6*gsl_matrix_get(FUNA,i,m);
              }
              
              
              temp2+=temp4*temp4+2*temp4*temp5+temp3;
            }
            
            point+=u;
          }
          
          gsl_vector_set(beta1,j,temp1/temp2);
        }
      }
      
      
      
      /* update H01 */
      
      
      for(p=0;p<a;p++)
      {
        gsl_bspline_eval(gsl_matrix_get(H01,0,p),B_spl,bw); 
        gsl_vector_memcpy(helpb,B_spl);
        MulM(sigmau_inv,helpb,B_spl);
        
        temp1=0;
        
        for(j=0;j<n;j++)
        {
          if(gsl_matrix_get(C,j,0)>=gsl_matrix_get(H01,0,p))
          {   
            temp2=0;
            for(u=0;u<q_eta;u++)  
            {
              temp2+=gsl_vector_get(eta,u)*gsl_matrix_get(C,j,2+u);
            }
            
            temp3=0;
            for(u=0;u<q_b;u++)  
              temp3+=gsl_vector_get(B_spl,u)*gsl_vector_get(theta,u);
            
            temp2+=temp3*(*gamma);              
            
            temp1+=exp(temp2)*gsl_matrix_get(FUNE,j,p);                
          }
        }
        
        gsl_matrix_set(H01,2,p,gsl_matrix_get(H01,1,p)/temp1);
        
      }
      
      
      
      
      /* update eta */
      if (survvar == 1) {
        gsl_vector *W=gsl_vector_alloc(q_eta);
        gsl_matrix *WW=gsl_matrix_alloc(q_eta,q_eta);
        gsl_vector *SW=gsl_vector_alloc(q_eta);
        gsl_matrix *SWW=gsl_matrix_alloc(q_eta,q_eta);
        
        gsl_vector_set_zero(SW);
        gsl_matrix_set_zero(SWW);
        
        
        for(j=0;j<n;j++)
        {
          for(u=0;u<q_eta;u++)    gsl_vector_set(W,u,gsl_matrix_get(C,j,2+u));
          
          if(gsl_matrix_get(C,j,1)==1)    gsl_vector_add(SW,W);
          
          for(p=0;p<a;p++)  
          {
            if(gsl_matrix_get(C,j,0)>=gsl_matrix_get(H01,0,p))
            {
              for(u=0;u<q_eta;u++)    gsl_vector_set(W,u,gsl_matrix_get(C,j,2+u));
              MulV(W,WW);
              
              gsl_bspline_eval(gsl_matrix_get(H01,0,p),B_spl,bw); 
              gsl_vector_memcpy(helpb,B_spl);
              MulM(sigmau_inv,helpb,B_spl);
              
              temp1=MulVV(eta,W);
              temp2=0;
              for(u=0;u<q_b;u++)    temp2+=gsl_vector_get(B_spl,u)*gsl_vector_get(theta,u);
              
              temp1+=temp2*(*gamma);
              
              gsl_vector_scale(W,gsl_matrix_get(H01,2,p)*exp(temp1)*gsl_matrix_get(FUNE,j,p));
              gsl_matrix_scale(WW,gsl_matrix_get(H01,2,p)*exp(temp1)*gsl_matrix_get(FUNE,j,p));
              
              gsl_vector_sub(SW,W);
              gsl_matrix_add(SWW,WW);
            }
          }
        }
        
        
        status=inv_matrix(SWW);
        if(status==100) return status;
        
        MulM(SWW,SW,W);
        
        for(p=0;p<q_eta;p++)   gsl_vector_set(eta,p,gsl_vector_get(eta,p)+gsl_vector_get(W,p));
        
        gsl_vector_free(W);
        gsl_matrix_free(WW);
        gsl_vector_free(SW);
        gsl_matrix_free(SWW);
        }
      
      /* update gamma  */
      
      temp1=0;
      temp2=0;
      
      for(j=0;j<n;j++)
      {
        if(gsl_matrix_get(C,j,1)==1)
        {
          gsl_bspline_eval(gsl_matrix_get(C,j,0),B_spl,bw);
          gsl_vector_memcpy(helpb,B_spl);
          MulM(sigmau_inv,helpb,B_spl);
          
          temp3=0;
          for(i=0;i<q_b;i++)    temp3+=gsl_vector_get(B_spl,i)*gsl_vector_get(theta,i);
          
          for(k=0;k<k_max;k++)
          {
            temp5=0;
            for(u=0;u<q_b;u++)    temp5+=gsl_vector_get(B_spl,u)*gsl_matrix_get(btheta,u,k);
            
            temp5*=gsl_matrix_get(FUNA,j,k);
            temp3+=temp5;
          }
          
          temp1+=temp3;
        }
        
        
        for(p=0;p<a;p++)  
        {
          if(gsl_matrix_get(C,j,0)>=gsl_matrix_get(H01,0,p))
          {
            gsl_bspline_eval(gsl_matrix_get(H01,0,p),B_spl,bw); 
            gsl_vector_memcpy(helpb,B_spl);
            MulM(sigmau_inv,helpb,B_spl);
            
            temp3=0;
            for(u=0;u<q_eta;u++)    temp3+=gsl_vector_get(eta,u)*gsl_matrix_get(C,j,2+u);
            
            temp4=0;
            for(u=0;u<q_b;u++)    temp4+=gsl_vector_get(B_spl,u)*gsl_vector_get(theta,u);
            
            temp3+=temp4*(*gamma);
            
            temp1-=gsl_matrix_get(H01,2,p)*exp(temp3)*(temp4*gsl_matrix_get(FUNE,j,p)+gsl_matrix_get(FUNAE,j,p));
            temp2+=gsl_matrix_get(H01,2,p)*exp(temp3)*(temp4*temp4*gsl_matrix_get(FUNE,j,p)
                                                         +2*temp4*gsl_matrix_get(FUNAE,j,p)+gsl_matrix_get(FUNA2E,j,p));
            
          }
        }
      }
      
      (*gamma)=(*gamma)+temp1/temp2;
      
      
      /* update theta */
      
      
      gsl_vector *Z=gsl_vector_alloc(q_b);
      gsl_matrix *ZZ=gsl_matrix_alloc(q_b,q_b);
      gsl_vector *SZ=gsl_vector_alloc(q_b);
      gsl_matrix *SZZ=gsl_matrix_alloc(q_b,q_b);
      
      gsl_vector_set_zero(SZ);
      gsl_matrix_set_zero(SZZ);
      
      point=0;
      for(i=0;i<n;i++)
      {
        u=(int)gsl_vector_get(M,i);
        for(t=0;t<u;t++) 
        {
          gsl_bspline_eval(gsl_matrix_get(Bio,point+t,0),B_spl,bw); 
          gsl_vector_memcpy(helpb,B_spl);
          MulM(sigmau_inv,helpb,B_spl);
          
          for(j=0;j<j_max;j++)
          {
            loc=0;
            for(m=0;m<j;m++) loc+=(int)gsl_vector_get(pdim,m);
            
            temp1=gsl_matrix_get(Bio,point+t,loc+j+1);
            
            for(k=0;k<(int)gsl_vector_get(pdim,j);k++)
              temp1-=gsl_matrix_get(Bio,point+t,loc+j+2+k)*gsl_matrix_get(beta0,j,k);
            
            temp2=0;
            for(k=0;k<q_b;k++)
              temp2+=gsl_vector_get(B_spl,k)*gsl_vector_get(theta,k);
            
            temp1=temp1-temp2*gsl_vector_get(beta1,j);
            
            temp2=0;
            for(m=0;m<k_max;m++)
            {
              temp3=0;
              for(k=0;k<q_b;k++)    temp3+=gsl_vector_get(B_spl,k)*gsl_matrix_get(btheta,k,m);
              temp3*=gsl_matrix_get(FUNA,i,m);
              temp2+=temp3;
            }
            temp2*=gsl_vector_get(beta1,j);
            
            temp1=(temp1-temp2)/gsl_vector_get(sigma,j);
            
            gsl_vector_memcpy(Z,B_spl);
            MulV(Z,ZZ);
            
            gsl_vector_scale(Z,temp1*gsl_vector_get(beta1,j));
            gsl_vector_add(SZ,Z);
            
            gsl_matrix_scale(ZZ,gsl_vector_get(beta1,j)*gsl_vector_get(beta1,j)/gsl_vector_get(sigma,j));
            gsl_matrix_add(SZZ,ZZ);
          }
        }
        
        point+=u;
        
        
        if(gsl_matrix_get(C,i,1)==1)
        {
          gsl_bspline_eval(gsl_matrix_get(C,i,0),B_spl,bw); 
          gsl_vector_memcpy(helpb,B_spl);
          MulM(sigmau_inv,helpb,B_spl);
          
          gsl_vector_scale(B_spl,*gamma);
          gsl_vector_add(SZ,B_spl);
        }
        
        for(p=0;p<a;p++)  
        {
          if(gsl_matrix_get(C,i,0)>=gsl_matrix_get(H01,0,p))
          {
            gsl_bspline_eval(gsl_matrix_get(H01,0,p),B_spl,bw);
            gsl_vector_memcpy(helpb,B_spl);
            MulM(sigmau_inv,helpb,B_spl);
            
            gsl_vector_memcpy(Z,B_spl);
            MulV(Z,ZZ);
            
            temp1=0;
            for(u=0;u<q_eta;u++)    temp1+=gsl_vector_get(eta,u)*gsl_matrix_get(C,i,2+u);
            
            temp2=0;
            for(u=0;u<q_b;u++)    temp2+=gsl_vector_get(B_spl,u)*gsl_vector_get(theta,u);
            
            temp1+=temp2*(*gamma);
            
            gsl_vector_scale(Z,gsl_matrix_get(H01,2,p)*exp(temp1)*gsl_matrix_get(FUNE,i,p)*(*gamma));
            gsl_matrix_scale(ZZ,gsl_matrix_get(H01,2,p)*exp(temp1)*gsl_matrix_get(FUNE,i,p)*(*gamma)*(*gamma));
            
            gsl_vector_sub(SZ,Z);
            gsl_matrix_add(SZZ,ZZ);
          }
        }
      }
      
      status=inv_matrix(SZZ);
      if(status==100) return status;
      
      MulM(SZZ,SZ,Z);
      for(p=0;p<q_b;p++)   gsl_vector_set(theta,p,gsl_vector_get(theta,p)+gsl_vector_get(Z,p));
      
      
      
      
      /* update btheta */
      
      gsl_matrix *prebtheta=gsl_matrix_alloc(btheta->size1,btheta->size2);
      
      do
      {
        gsl_matrix_memcpy(prebtheta,btheta);
        
        for(k=0;k<k_max;k++)
        {
          gsl_vector_set_zero(SZ);
          gsl_matrix_set_zero(SZZ);
          
          point=0;
          for(i=0;i<n;i++)
          {
            u=(int)gsl_vector_get(M,i);
            for(t=0;t<u;t++) 
            {
              for(j=0;j<j_max;j++)
              {
                gsl_bspline_eval(gsl_matrix_get(Bio,point+t,0),B_spl,bw); 
                gsl_vector_memcpy(helpb,B_spl);
                MulM(sigmau_inv,helpb,B_spl);
                
                gsl_vector_memcpy(Z,B_spl);
                MulV(Z,ZZ);
                
                loc=0;
                for(m=0;m<j;m++) loc+=(int)gsl_vector_get(pdim,m);
                
                temp1=gsl_matrix_get(Bio,point+t,loc+j+1);
                
                for(q=0;q<(int)gsl_vector_get(pdim,j);q++)
                  temp1-=gsl_matrix_get(Bio,point+t,loc+j+2+q)*gsl_matrix_get(beta0,j,q);
                
                temp2=0;
                for(q=0;q<q_b;q++)
                  temp2+=gsl_vector_get(B_spl,q)*gsl_vector_get(theta,q);
                
                temp1-=temp2*gsl_vector_get(beta1,j);
                temp1=temp1*gsl_vector_get(beta1,j)*gsl_matrix_get(FUNA,i,k)/gsl_vector_get(sigma,j);
                
                temp2=0;
                for(m=0;m<k_max;m++)
                {
                  temp3=0;
                  for(q=0;q<q_b;q++)
                    temp3+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,m);
                  
                  if(m==k)
                    temp3*=gsl_matrix_get(FUNA2,i,k);
                  if(m>k)
                    temp3*=gsl_matrix_get(FUNA2,i,k_max+m*(m-1)/2+k);
                  if(m<k)
                    temp3*=gsl_matrix_get(FUNA2,i,k_max+k*(k-1)/2+m);
                  
                  temp2+=temp3;
                }
                
                temp2*=gsl_pow_2(gsl_vector_get(beta1,j))/gsl_vector_get(sigma,j);
                
                gsl_vector_scale(Z,temp1-temp2);
                gsl_vector_add(SZ,Z);
                
                gsl_matrix_scale(ZZ,gsl_matrix_get(FUNA2,i,k)*gsl_pow_2(gsl_vector_get(beta1,j))/gsl_vector_get(sigma,j));
                gsl_matrix_add(SZZ,ZZ);
                
              }
            }
            
            point+=u;
            
            
            if(gsl_matrix_get(C,i,1)==1)
            {
              gsl_bspline_eval(gsl_matrix_get(C,i,0),B_spl,bw); 
              gsl_vector_memcpy(helpb,B_spl);
              MulM(sigmau_inv,helpb,B_spl);
              
              gsl_vector_scale(B_spl,(*gamma)*gsl_matrix_get(FUNA,i,k));
              gsl_vector_add(SZ,B_spl);
            }
            
            for(p=0;p<a;p++)  
            {
              if(gsl_matrix_get(C,i,0)>=gsl_matrix_get(H01,0,p))
              {
                gsl_bspline_eval(gsl_matrix_get(H01,0,p),B_spl,bw); 
                gsl_vector_memcpy(helpb,B_spl);
                MulM(sigmau_inv,helpb,B_spl);
                
                gsl_vector_memcpy(Z,B_spl);
                MulV(Z,ZZ);
                
                temp1=0;
                for(u=0;u<q_eta;u++)
                {
                  temp1+=gsl_vector_get(eta,u)*gsl_matrix_get(C,i,2+u);
                }
                
                temp2=0;
                for(u=0;u<q_b;u++)
                  temp2+=gsl_vector_get(B_spl,u)*gsl_vector_get(theta,u);
                
                temp1+=temp2*(*gamma);
                
                gsl_vector_scale(Z,gsl_matrix_get(H01,2,p)*exp(temp1)*gsl_matrix_get(FUNAEV,i,p*k_max+k)*(*gamma));
                gsl_matrix_scale(ZZ,gsl_matrix_get(H01,2,p)*exp(temp1)*gsl_matrix_get(FUNA2EV,i,p*k_max+k)*(*gamma)*(*gamma));
                
                gsl_vector_sub(SZ,Z);
                gsl_matrix_add(SZZ,ZZ);
              }
            }
          }
          
          
          status=inv_matrix(SZZ);
          if(status==100) return status;
          
          MulM(SZZ,SZ,Z);
          
          for(p=0;p<q_b;p++)   gsl_matrix_set(btheta,p,k,gsl_matrix_get(btheta,p,k)+gsl_vector_get(Z,p));
        }
        
        
      }while(DiffM(prebtheta,btheta)==1);
      
      
      
      gsl_vector_free(Z);
      gsl_matrix_free(ZZ);
      gsl_vector_free(SZ);
      gsl_matrix_free(SZZ);
      gsl_matrix_free(prebtheta);
      
      
      for(k=0;k<k_max;k++)
      {
        temp1=0;
        for(q=0;q<q_b;q++)
          temp1+=gsl_pow_2(gsl_matrix_get(btheta,q,k));
        
        gsl_vector_set(sigmad,k,gsl_vector_get(sigmad,k)*temp1);
        
        for(q=0;q<q_b;q++)
          gsl_matrix_set(btheta,q,k,gsl_matrix_get(btheta,q,k)/sqrt(temp1));
      }
      
      
      
      
      gsl_matrix_free(FUNA);
      gsl_matrix_free(FUNA2); 
      gsl_matrix_free(FUNE); 
      gsl_matrix_free(FUNAE);
      gsl_matrix_free(FUNA2E);
      gsl_matrix_free(FUNAEV);
      gsl_matrix_free(FUNA2EV);
      
      gsl_vector_free(B_spl);
      gsl_vector_free(helpb);
      
      
      return 1;
      
    }





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
    )
    {
      int p_max=beta0->size2;
      int j_max=beta0->size1;
      int q_b=theta->size;
      int k_max=sigmad->size;
      int q_eta=eta->size;
      int a=H01->size2;
      int n=M->size;
      
      int i,j,t,k,u,m,p,q,point,loc,status;
      double temp1,temp2,temp3,temp4,temp5,temp6,temp7;
      
      
      gsl_matrix_set_zero(FUNA);
      gsl_matrix_set_zero(FUNA2);
      gsl_matrix_set_zero(FUNE);
      gsl_matrix_set_zero(FUNAE);
      gsl_matrix_set_zero(FUNA2E);
      gsl_matrix_set_zero(FUNAEV);
      gsl_matrix_set_zero(FUNA2EV);
      
      int da1,da2,da3,da4;
      
      
      gsl_vector *xi = gsl_vector_alloc(quadpoint);
      gsl_vector *wi = gsl_vector_alloc(quadpoint);
      gsl_vector *ti=gsl_vector_alloc(k_max);
      gsl_vector *B_spl=gsl_vector_alloc(q_b);
      gsl_vector *helpb=gsl_vector_alloc(q_b);
      
      
      for(i=0;i<quadpoint/2;i++)   gsl_vector_set(xi,i,xs[i]);
      for(i=0;i<quadpoint/2;i++)   gsl_vector_set(wi,i,ws[i]);
      
      for(i=0;i<quadpoint/2;i++)   gsl_vector_set(xi,quadpoint-i-1, 0-gsl_vector_get(xi,i));
      for(i=0;i<quadpoint/2;i++)   gsl_vector_set(wi,quadpoint-i-1, gsl_vector_get(wi,i));
      
      
      
      point=0;    
      for(i=0;i<n;i++)
      {
        temp1=0;
        
        u=(int)gsl_vector_get(M,i);
        
        for(da1=0;da1<quadpoint;da1++)
        {
          for(da2=0;da2<quadpoint;da2++)
          {
            /*
             for(da3=0;da3<quadpoint;da3++)
             {
             */
            gsl_vector_set(ti,0,gsl_vector_get(xi,da1)*sqrt(2*gsl_vector_get(sigmad,0)));
            gsl_vector_set(ti,1,gsl_vector_get(xi,da2)*sqrt(2*gsl_vector_get(sigmad,1)));
            /*
             gsl_vector_set(ti,2,gsl_vector_get(xi,da3)*sqrt(2*gsl_vector_get(sigmad,2)));
             */
            temp2=exp(10);
            
            for(t=0;t<u;t++) 
            {
              gsl_bspline_eval(gsl_matrix_get(Bio,point+t,0),B_spl,bw); 
              gsl_vector_memcpy(helpb,B_spl);
              MulM(sigmau_inv,helpb,B_spl);
              
              for(j=0;j<j_max;j++)
              {
                loc=0;
                for(m=0;m<j;m++) loc+=(int)gsl_vector_get(pdim,m);
                
                temp3=gsl_matrix_get(Bio,point+t,loc+j+1);
                
                for(q=0;q<(int)gsl_vector_get(pdim,j);q++)
                  temp3-=gsl_matrix_get(Bio,point+t,loc+j+2+q)*gsl_matrix_get(beta0,j,q);
                
                temp4=0;
                for(q=0;q<q_b;q++)
                  temp4+=gsl_vector_get(B_spl,q)*gsl_vector_get(theta,q);
                
                temp3-=temp4*gsl_vector_get(beta1,j);
                
                temp4=0;
                for(m=0;m<k_max;m++)
                {
                  temp5=0;
                  for(q=0;q<q_b;q++)
                    temp5+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,m);
                  temp5*=gsl_vector_get(ti,m);
                  
                  temp4+=temp5;
                }
                
                temp4*=gsl_vector_get(beta1,j);
                temp3-=temp4;
                
                temp2*=exp(-1/(2*gsl_vector_get(sigma,j))*gsl_pow_2(temp3))/sqrt(2*M_PI*gsl_vector_get(sigma,j));
              }
            }
            
            if(gsl_matrix_get(C,i,1)==1)  
            {
              gsl_bspline_eval(gsl_matrix_get(C,i,0),B_spl,bw);
              gsl_vector_memcpy(helpb,B_spl);
              MulM(sigmau_inv,helpb,B_spl);
              
              temp3=0;
              for(q=0;q<q_eta;q++)
              {
                temp3+=gsl_vector_get(eta,q)*gsl_matrix_get(C,i,2+q);
              }
              
              temp4=0;
              for(q=0;q<q_b;q++)
                temp4+=gsl_vector_get(B_spl,q)*gsl_vector_get(theta,q);
              temp4*=gamma;
              
              temp5=0;
              for(m=0;m<k_max;m++)
              {
                temp6=0;
                for(q=0;q<q_b;q++)
                  temp6+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,m);
                temp6*=gsl_vector_get(ti,m);
                
                temp5+=temp6;
              }
              
              temp5*=gamma;
              
              temp2*=HAZ(H01,gsl_matrix_get(C,i,0))*exp(temp3+temp4+temp5);
            }
            
            
            temp7=0;
            for(p=0;p<a;p++)  
            {
              if(gsl_matrix_get(C,i,0)>=gsl_matrix_get(H01,0,p))
              {
                gsl_bspline_eval(gsl_matrix_get(H01,0,p),B_spl,bw); 
                gsl_vector_memcpy(helpb,B_spl);
                MulM(sigmau_inv,helpb,B_spl);
                
                temp3=0;
                for(q=0;q<q_eta;q++)
                {
                  temp3+=gsl_vector_get(eta,q)*gsl_matrix_get(C,i,2+q);
                }
                
                temp4=0;
                for(q=0;q<q_b;q++)
                  temp4+=gsl_vector_get(B_spl,q)*gsl_vector_get(theta,q);
                temp4*=gamma;
                
                temp5=0;
                for(m=0;m<k_max;m++)
                {
                  temp6=0;
                  for(q=0;q<q_b;q++)
                    temp6+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,m);
                  temp6*=gsl_vector_get(ti,m);
                  
                  temp5+=temp6;
                }
                temp5*=gamma;
                
                temp7+=gsl_matrix_get(H01,2,p)*exp(temp3+temp4+temp5);
              }
            }
            
            temp2*=exp(0-temp7);
            /*
             temp2*=gsl_vector_get(wi,da1);
             */
            temp2*=gsl_vector_get(wi,da1)*gsl_vector_get(wi,da2);
            /*
             temp2*=gsl_vector_get(wi,da1)*gsl_vector_get(wi,da2)*gsl_vector_get(wi,da3);
             */
            
            temp1+=temp2;
            
            for(m=0;m<k_max;m++)
            {
              gsl_matrix_set(FUNA,i,m,gsl_matrix_get(FUNA,i,m)+temp2*gsl_vector_get(ti,m));
              
              for(k=0;k<=m;k++)
              {
                if(k==m) gsl_matrix_set(FUNA2,i,k,gsl_matrix_get(FUNA2,i,k)+temp2*gsl_pow_2(gsl_vector_get(ti,k)));
                else
                {
                  gsl_matrix_set(FUNA2,i,k_max+m*(m-1)/2+k,gsl_matrix_get(FUNA2,i,k_max+m*(m-1)/2+k)+temp2*
                    gsl_vector_get(ti,k)*gsl_vector_get(ti,m));
                }
              }
            }
            
            for(p=0;p<a;p++)
            {
              gsl_bspline_eval(gsl_matrix_get(H01,0,p),B_spl,bw); 
              gsl_vector_memcpy(helpb,B_spl);
              MulM(sigmau_inv,helpb,B_spl);
              
              temp3=0;
              for(m=0;m<k_max;m++)
              {
                temp4=0;
                for(q=0;q<q_b;q++)
                  temp4+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,m);
                temp3+=temp4*gsl_vector_get(ti,m);
              }
              
              gsl_matrix_set(FUNE,i,p,gsl_matrix_get(FUNE,i,p)+temp2*exp(temp3*gamma));
              gsl_matrix_set(FUNAE,i,p,gsl_matrix_get(FUNAE,i,p)+temp2*exp(temp3*gamma)*temp3);
              gsl_matrix_set(FUNA2E,i,p,gsl_matrix_get(FUNA2E,i,p)+temp2*exp(temp3*gamma)*temp3*temp3);
              
              for(m=0;m<k_max;m++)
              {
                gsl_matrix_set(FUNAEV,i,p*k_max+m,gsl_matrix_get(FUNAEV,i,p*k_max+m)+temp2*exp(temp3*gamma)*
                  gsl_vector_get(ti,m));
                gsl_matrix_set(FUNA2EV,i,p*k_max+m,gsl_matrix_get(FUNA2EV,i,p*k_max+m)+temp2*exp(temp3*gamma)*
                  gsl_pow_2(gsl_vector_get(ti,m)));
              }
            }
            
          }
        }
        /*
             }
         
         */
        
        if(temp1==0) return 100;
        
        
        for(k=0;k<FUNA->size2;k++)
          gsl_matrix_set(FUNA,i,k,gsl_matrix_get(FUNA,i,k)/temp1);
        
        for(k=0;k<FUNA2->size2;k++)
          gsl_matrix_set(FUNA2,i,k,gsl_matrix_get(FUNA2,i,k)/temp1);
        
        for(k=0;k<FUNE->size2;k++)
        {
          gsl_matrix_set(FUNE,i,k,gsl_matrix_get(FUNE,i,k)/temp1);
          gsl_matrix_set(FUNAE,i,k,gsl_matrix_get(FUNAE,i,k)/temp1);
          gsl_matrix_set(FUNA2E,i,k,gsl_matrix_get(FUNA2E,i,k)/temp1);
        }
        
        for(k=0;k<FUNAEV->size2;k++)
        {
          gsl_matrix_set(FUNAEV,i,k,gsl_matrix_get(FUNAEV,i,k)/temp1);
          gsl_matrix_set(FUNA2EV,i,k,gsl_matrix_get(FUNA2EV,i,k)/temp1);
        }
        
        point+=u;
      }
      
      
      
      
      
      gsl_vector_free(xi);
      gsl_vector_free(wi);
      gsl_vector_free(ti);
      gsl_vector_free(B_spl);
      gsl_vector_free(helpb);
      
      
      
      return 0;
    }


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
        )
    {
      int p_max=beta0->size2;
      int j_max=beta0->size1;
      int q_b=theta->size;
      int k_max=sigmad->size;
      int q_eta=eta->size;
      int a=H01->size2;
      int n=M->size;
      
      int i,j,t,k,u,m,p,q,point,loc,status;
      double temp1,temp2,temp3,temp4,temp5,temp6,temp7;
      
      
      int da1,da2,da3,da4;
      
      
      gsl_vector *xi = gsl_vector_alloc(quadpoint);
      gsl_vector *wi = gsl_vector_alloc(quadpoint);
      gsl_vector *ti=gsl_vector_alloc(k_max);
      gsl_vector *B_spl=gsl_vector_alloc(q_b);
      gsl_vector *helpb=gsl_vector_alloc(q_b);
      
      
      for(i=0;i<quadpoint/2;i++)   gsl_vector_set(xi,i,xs[i]);
      for(i=0;i<quadpoint/2;i++)   gsl_vector_set(wi,i,ws[i]);
      
      for(i=0;i<quadpoint/2;i++)   gsl_vector_set(xi,quadpoint-i-1, 0-gsl_vector_get(xi,i));
      for(i=0;i<quadpoint/2;i++)   gsl_vector_set(wi,quadpoint-i-1, gsl_vector_get(wi,i));
      
      double loglike=0;
      
      point=0;    
      for(i=0;i<n;i++)
      {
        temp1=0;
        
        u=(int)gsl_vector_get(M,i);
        
        for(da1=0;da1<quadpoint;da1++)
        {
          for(da2=0;da2<quadpoint;da2++)
          {
            /*
             for(da3=0;da3<quadpoint;da3++)
             {
             */
            
            gsl_vector_set(ti,0,gsl_vector_get(xi,da1)*sqrt(2*gsl_vector_get(sigmad,0)));
            gsl_vector_set(ti,1,gsl_vector_get(xi,da2)*sqrt(2*gsl_vector_get(sigmad,1)));
            /*
             gsl_vector_set(ti,2,gsl_vector_get(xi,da3)*sqrt(2*gsl_vector_get(sigmad,2)));
             */
            temp2=1/sqrt(M_PI);
            
            for(t=0;t<u;t++) 
            {
              gsl_bspline_eval(gsl_matrix_get(Bio,point+t,0),B_spl,bw); 
              gsl_vector_memcpy(helpb,B_spl);
              MulM(sigmau_inv,helpb,B_spl);
              
              for(j=0;j<j_max;j++)
              {
                loc=0;
                for(m=0;m<j;m++) loc+=(int)gsl_vector_get(pdim,m);
                
                temp3=gsl_matrix_get(Bio,point+t,loc+j+1);
                
                for(q=0;q<(int)gsl_vector_get(pdim,j);q++)
                  temp3-=gsl_matrix_get(Bio,point+t,loc+j+2+q)*gsl_matrix_get(beta0,j,q);
                
                temp4=0;
                for(q=0;q<q_b;q++)
                  temp4+=gsl_vector_get(B_spl,q)*gsl_vector_get(theta,q);
                
                temp3-=temp4*gsl_vector_get(beta1,j);
                
                temp4=0;
                for(m=0;m<k_max;m++)
                {
                  temp5=0;
                  for(q=0;q<q_b;q++)
                    temp5+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,m);
                  temp5*=gsl_vector_get(ti,m);
                  
                  temp4+=temp5;
                }
                
                temp4*=gsl_vector_get(beta1,j);
                temp3-=temp4;
                
                temp2*=exp(-1/(2*gsl_vector_get(sigma,j))*gsl_pow_2(temp3))
                  /sqrt(2*M_PI*gsl_vector_get(sigma,j));
              }
            }
            
            
            
            if(gsl_matrix_get(C,i,1)==1)  
            {
              gsl_bspline_eval(gsl_matrix_get(C,i,0),B_spl,bw);
              gsl_vector_memcpy(helpb,B_spl);
              MulM(sigmau_inv,helpb,B_spl);
              
              temp3=0;
              for(q=0;q<q_eta;q++)
              {
                temp3+=gsl_vector_get(eta,q)*gsl_matrix_get(C,i,2+q);
                
              }
              
              temp4=0;
              for(q=0;q<q_b;q++)
                temp4+=gsl_vector_get(B_spl,q)*gsl_vector_get(theta,q);
              temp4*=gamma;
              
              temp5=0;
              for(m=0;m<k_max;m++)
              {
                temp6=0;
                for(q=0;q<q_b;q++)
                  temp6+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,m);
                temp6*=gsl_vector_get(ti,m);
                
                temp5+=temp6;
              }
              
              temp5*=gamma;
              
              temp2*=HAZ(H01,gsl_matrix_get(C,i,0))*exp(temp3+temp4+temp5);
              
            }
            
            
            
            temp7=0;
            for(p=0;p<a;p++)  
            {
              if(gsl_matrix_get(C,i,0)>=gsl_matrix_get(H01,0,p))
              {
                gsl_bspline_eval(gsl_matrix_get(H01,0,p),B_spl,bw); 
                gsl_vector_memcpy(helpb,B_spl);
                MulM(sigmau_inv,helpb,B_spl);
                
                temp3=0;
                for(q=0;q<q_eta;q++)
                {
                  temp3+=gsl_vector_get(eta,q)*gsl_matrix_get(C,i,2+q);
                }
                
                temp4=0;
                for(q=0;q<q_b;q++)
                  temp4+=gsl_vector_get(B_spl,q)*gsl_vector_get(theta,q);
                temp4*=gamma;
                
                temp5=0;
                for(m=0;m<k_max;m++)
                {
                  temp6=0;
                  for(q=0;q<q_b;q++)
                    temp6+=gsl_vector_get(B_spl,q)*gsl_matrix_get(btheta,q,m);
                  temp6*=gsl_vector_get(ti,m);
                  
                  temp5+=temp6;
                }
                temp5*=gamma;
                
                temp7+=gsl_matrix_get(H01,2,p)*exp(temp3+temp4+temp5);
              }
            }
            
            temp2*=exp(0-temp7);
            /*
             temp2*=gsl_vector_get(wi,da1);
             */
            temp2*=gsl_vector_get(wi,da1)*gsl_vector_get(wi,da2);
            /*
             temp2*=gsl_vector_get(wi,da1)*gsl_vector_get(wi,da2)*gsl_vector_get(wi,da3);
             
             */
            temp1+=temp2;
            
            
          }
        }
        /*
             }
         */
        
        loglike+=log(temp1);
        
        point+=u;
      }
      
      
      
      
      
      gsl_vector_free(xi);
      gsl_vector_free(wi);
      gsl_vector_free(ti);
      gsl_vector_free(B_spl);
      gsl_vector_free(helpb);
      
      
      return loglike;
    }


    double SmoothHAZ(const gsl_matrix *H, double t)
    {
      int a=H->size2;
      int i;
      double temp1=0,temp2=0;
      double bw=1;
      
      for(i=0;i<a;i++)
      {
        temp1+=gsl_matrix_get(H,2,i)*gsl_ran_gaussian_pdf((t-gsl_matrix_get(H,0,i))/bw,1);
        temp2+=gsl_ran_gaussian_pdf((t-gsl_matrix_get(H,0,i))/bw,1);
      }
      
      return (temp1/temp2);
    }


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
        const double gamma,
        const gsl_vector *presigmad,
        const gsl_vector *sigmad,
        const int conversigmad
    )
    {    
      
      double epsilon=0.0001;  
      
      if (conversigmad == 0) {
        
        if(DiffM(prebeta0,beta0)>epsilon || DiffV(prebeta1,beta1)>epsilon || DiffV(presigma,sigma)>epsilon ||
           DiffV(pretheta,theta)>epsilon || DiffV(preeta,eta)>epsilon || Abs(pregamma,gamma)>epsilon)
          
          return 1;
        
        else return 0;
      } else{
        
        if(DiffM(prebeta0,beta0)>epsilon || DiffV(prebeta1,beta1)>epsilon || DiffV(presigma,sigma)>epsilon ||
           DiffV(pretheta,theta)>epsilon || DiffV(preeta,eta)>epsilon || Abs(pregamma,gamma)>epsilon 
             || DiffV(presigmad,sigmad)>epsilon)
          
          return 1;
        
        else return 0;
        
        }
      
      
    }
    
    
    
    
    
    int inv_matrix(gsl_matrix *x_square)
    {
      int i,j;
      int k = x_square->size1;
      
      int status;
      
      gsl_vector *temp_vector=gsl_vector_alloc(k),
        *solution=gsl_vector_alloc(k);
      gsl_matrix *out = gsl_matrix_alloc(k,k);
      
      for(i=0;i<k;i++)
      {
        for(j=0;j<k;j++) gsl_matrix_set(out,i,j,gsl_matrix_get(x_square,i,j));
      }
      
      status=gsl_linalg_cholesky_decompn(out);
      if(status==100) 
      {
        gsl_vector_free(temp_vector);
        gsl_vector_free(solution);
        gsl_matrix_free(out);
        
        return status;
      }
      
      for (i = 0; i < k; i++)
      {
        gsl_vector_set_all(temp_vector,0);
        gsl_vector_set(temp_vector,i,1);
        
        status=gsl_linalg_cholesky_solven(out, temp_vector, solution);
        
        if(status==100) 
        {
          gsl_vector_free(temp_vector);
          gsl_vector_free(solution);
          gsl_matrix_free(out);
          
          return status;
        }
        
        gsl_matrix_set_col(x_square,i,solution);
      }
      
      gsl_vector_free(temp_vector);
      gsl_vector_free(solution);
      gsl_matrix_free(out);
      
      return 0;
      
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
    
    
    
    
    
    void TransM(const gsl_matrix *A, gsl_matrix *B)
    {
      int rowa = A->size1;
      int cola = A->size2;
      
      int i, j;
      
      for(i=0;i<rowa;i++)
      {
        for(j=0;j<cola;j++)
        {
          gsl_matrix_set(B,j,i,gsl_matrix_get(A,i,j));
        }
      }
      
    }   
    
    
    
    double CH(const gsl_matrix *H, double t)
    {
      int a=H->size2;
      int i;
      
      double ch;
      
      if(t<gsl_matrix_get(H,0,0)) ch=0;
      else
      {
        ch=0;
        i=0;
        do
        {
          ch+=gsl_matrix_get(H,2,i);
          i+=1;
        }while(i<=a-1 && t>=gsl_matrix_get(H,0,i));
      }
      
      return (ch);
    }
    
    
    double HAZ(const gsl_matrix *H, double t)
    {
      int a=H->size2;
      int i;
      double temp=0;
      
      for(i=0;i<a;i++)
      {
        if(t==gsl_matrix_get(H,0,i)) temp=gsl_matrix_get(H,2,i);
      }
      
      return (temp);
    }       
    
    
    
    
    double Min(const double t1, const double t2)
    {
      if(t1<t2) return t1;
      else return t2;
    }
    
    
    
    
    double Abs(const double a, const double b)
    {
      
      if (a>=b) return a-b;
      else return b-a;
    }
    
    
    double DiffV(const gsl_vector *veca, const gsl_vector *vecb)
    {
      int k=veca->size;
      int i;
      double diff=0;
      
      for(i=0;i<k;i++)
      {
        if(Abs(gsl_vector_get(veca,i),gsl_vector_get(vecb,i))>diff) diff=Abs(gsl_vector_get(veca,i),gsl_vector_get(vecb,i));
      }
      
      return (diff);   
    }
    
    
    double DiffM(const gsl_matrix *matrixa, const gsl_matrix *matrixb)
    {
      int nrow=matrixa->size1, ncol=matrixa->size2;
      int i, j;
      double diff=0;
      
      double epsilon=0.0001;  
      
      for(i=0;i<nrow;i++)
      {
        for(j=0;j<ncol;j++)
        {
          if(Abs(gsl_matrix_get(matrixa,i,j),gsl_matrix_get(matrixb,i,j))>diff) 
            diff=Abs(gsl_matrix_get(matrixa,i,j),gsl_matrix_get(matrixb,i,j));
        }
      }
      
      if(diff>epsilon) return 1;
      else return 0;
      
    }
    
    
    
    void STAT(const gsl_matrix *store,int i,double *mean,double *sd)
    {
      int n=store->size1;
      int j;
      
      *mean=0, *sd=0;
      
      for(j=0;j<n;j++)   *mean+=gsl_matrix_get(store,j,i);
      *mean=*mean/(double)n;
      
      for(j=0;j<n;j++)   *sd+=gsl_pow_2(gsl_matrix_get(store,j,i)-*mean);
      *sd = sqrt(*sd/(double)(n-1));
      
    }
    
    
    
    
    int GetN(double t)
    {
      
      return (int)(t+1);
      
    }
    
    
    Rcpp::List jmspline_cmain(int n, int n_total, double tL, double tU,int p01, 
                              int p02, int q_b, int q_eta, int j_max, int t_max,
                              int nbreak, int k_max, int quadpoint, int maxiter, 
                              int trace, std::string ydatanew, std::string mdatanew,
                              std::string cdatanew, std::string sigmau_invnew,
                              std::string tbthetanew, 
                              std::vector<double> xs, 
                              std::vector<double> ws,
                              std::string beta0initnew, std::string beta1initnew,
                              std::string sigmainitnew, std::string thetainitnew,
                              std::string sigmadinitnew,double gamma,
                              int survvar, int conversigmad)
    {
      
      /* allocate space for data */
      int k_cubic = 4;
      int pt = p01+p02;
      /*Define a vector of # of covariates for all biomarkers; 
       the current version can only allow two biomarkers*/
      gsl_vector *pmax = gsl_vector_calloc(2);
      gsl_vector_set(pmax, 0, p01);
      gsl_vector_set(pmax, 1, p02);
      int p_max = gsl_vector_max(pmax);
      
      gsl_matrix *sigmau_inv=gsl_matrix_alloc(nbreak+k_cubic-2,nbreak+k_cubic-2);
      gsl_vector *M= gsl_vector_alloc(n);                               /*# obs per subject */
      gsl_matrix *C = gsl_matrix_alloc(n,q_eta+2);                      /*data for event times */
      gsl_matrix *Bio = gsl_matrix_alloc(n_total,j_max+pt+1);           /*data for biomarkers */ 
      gsl_matrix *btheta = gsl_matrix_alloc(q_b,k_max);
      
      /* allocate space for estimated parameters */    
      gsl_matrix *beta0 = gsl_matrix_alloc(j_max,p_max);
      gsl_vector *beta1 = gsl_vector_alloc(j_max),
        *sigma = gsl_vector_alloc(j_max),
        *theta = gsl_vector_alloc(q_b),
        *sigmad = gsl_vector_alloc(k_max),
        *eta = gsl_vector_alloc(q_eta);
      
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
          if (n_total == nrows)
          {   rewind(f);
            gsl_matrix_fscanf(f, Bio);
            fclose(f);
          }
          else
          {
            Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                    n_total,
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
          if (q_b == nrows)
          {   rewind(f);
            gsl_matrix_fscanf(f, btheta);
            fclose(f);
          }
          else
          {
            Rprintf("Input oberservations is %d, but the number of rows in %s is %d",
                    q_b,
                    tbthetanew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      /* read beta0 matrix  */
      {
        FILE * f = fopen(beta0initnew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", beta0initnew.c_str());
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
                    beta0initnew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      /* read beta1 vector  */
      {
        FILE * f = fopen(beta1initnew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", beta1initnew.c_str());
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
                    beta1initnew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      /* read sigma vector  */
      {
        FILE * f = fopen(sigmainitnew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", sigmainitnew.c_str());
          return R_NilValue;
        }
        
        
        int nrows=0;
        // Extract characters from file and store in character c
        for (char c = fgetc(f); c != EOF; c = fgetc(f))
          if (c == '\n')  nrows = nrows + 1;
          nrows=nrows+1;
          if (j_max==nrows)
          {   rewind(f);
            gsl_vector_fscanf(f, sigma);
            fclose(f);
          }
          else
          {
            Rprintf("Input subjects is %d, but the number of rows in %s is %d",j_max,
                    sigmainitnew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      /* read theta vector  */
      {
        FILE * f = fopen(thetainitnew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", thetainitnew.c_str());
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
                    thetainitnew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      /* read sigmad vector  */
      {
        FILE * f = fopen(sigmadinitnew.c_str(), "r");
        
        if (f == NULL)
        {
          Rprintf("File %s does not exist.\n", sigmadinitnew.c_str());
          return R_NilValue;
        }
        
        
        int nrows=0;
        // Extract characters from file and store in character c
        for (char c = fgetc(f); c != EOF; c = fgetc(f))
          if (c == '\n')  nrows = nrows + 1;
          nrows=nrows+1;
          if (k_max==nrows)
          {   rewind(f);
            gsl_vector_fscanf(f, sigmad);
            fclose(f);
          }
          else
          {
            Rprintf("Input subjects is %d, but the number of rows in %s is %d",k_max,
                    sigmadinitnew.c_str(),nrows);
            fclose(f);
            return R_NilValue;
          }
      }
      
      gsl_vector *pdim=gsl_vector_alloc(j_max);
      gsl_vector_set(pdim,0,p01);
      gsl_vector_set(pdim,1,p02);
      
      int loca=p01+p02;
        
        /* allocate space for pre parameters */            
        gsl_matrix *prebeta0 = gsl_matrix_alloc(j_max,p_max);
        gsl_vector *prebeta1 = gsl_vector_alloc(j_max),
                   *presigma = gsl_vector_alloc(j_max),
                   *pretheta = gsl_vector_alloc(q_b),
                   *presigmad = gsl_vector_alloc(k_max),
                   *preeta = gsl_vector_alloc(q_eta);
          
        gsl_matrix *prebtheta = gsl_matrix_alloc(q_b,k_max);
          
        double pregamma;   
        
        // /* assign the start values of parameters */
        // 
        // gsl_matrix_set(beta0,0,0,0);
        // gsl_matrix_set(beta0,1,0,0);
        // 
        // gsl_vector_set(beta1,0,1);
        // gsl_vector_set(beta1,1,1.28);
        // 
        // 
        // gsl_vector_set(sigma,0,3.1);
        // gsl_vector_set(sigma,1,3.8);
        // 
        // 
        // gsl_vector_set(theta,0,15);
        // gsl_vector_set(theta,1,-0.1);
        // gsl_vector_set(theta,2,-6);
        // gsl_vector_set(theta,3,-0.1);
        // gsl_vector_set(theta,4,1.1);
        // gsl_vector_set(theta,5,0.1);
        // gsl_vector_set(theta,6,-0.24);
        // gsl_vector_set(theta,7,-0.4);
        // gsl_vector_set(theta,8,-1.6);
        // gsl_vector_set(theta,9,-0.1);
        // 
        // gsl_vector_set(sigmad,0,8);
        // gsl_vector_set(sigmad,1,5);
        // 
        // gsl_vector_set_zero(eta);
        // 
        // gamma=0.3;
        gsl_vector_set_zero(eta);
        
        // const gsl_rng_type * T;
        // gsl_rng * r;
        // gsl_rng_env_setup();
        // T = gsl_rng_default;
        // r = gsl_rng_alloc (T);
        
        
        gsl_bspline_workspace *bw=gsl_bspline_alloc(k_cubic,nbreak);
        gsl_bspline_knots_uniform(tL,tU,bw);
        gsl_vector *B_spl=gsl_vector_alloc(nbreak+k_cubic-2);
        
        
        gsl_matrix *V=gsl_matrix_alloc(k_max,k_max);
        gsl_vector *S=gsl_vector_alloc(q_b);
        gsl_vector *work=gsl_vector_alloc(q_b);
        gsl_matrix *bthetat=gsl_matrix_alloc(k_max,q_b);
        gsl_matrix *hbtheta=gsl_matrix_alloc(q_b,k_max);
        gsl_matrix *reta=gsl_matrix_alloc(q_b,q_b);
        gsl_matrix *veta=gsl_matrix_alloc(q_b,q_b);
        
        int point,i,j,t,k,ni,iter,status,loc,p,q;
        
        for(i=0;i<C->size1;i++)
        {
          if(gsl_matrix_get(C,i,0)>tU) 
          {
            gsl_matrix_set(C,i,0,tU);	
            gsl_matrix_set(C,i,1,0);
          }
          
        }
        
        gsl_matrix * FH01 = gsl_matrix_alloc(2,n);         /* baseline hazard function for event risk */
        
        gsl_matrix_set_zero(FH01);
        
        
        /* find # of events */
        int u, v, a=0;
        
        for(j=0;j<n;j++)
        {
          u=0;
          if(gsl_matrix_get(C,j,1)==1) 
          {
            if(a>=1)
            {
              for(t=0;t<a;t++)
              {
                if(gsl_matrix_get(FH01,0,t)==gsl_matrix_get(C,j,0))
                {
                  gsl_matrix_set(FH01,1,t,gsl_matrix_get(FH01,1,t)+1);
                  u=1;
                }
              }
              if(u==0)
              {
                t=-1;
                
                do
                {
                  t=t+1;
                }while(t<=a-1 && gsl_matrix_get(FH01,0,t) < gsl_matrix_get(C,j,0));
                
                if(t==a)
                {
                  gsl_matrix_set(FH01,0,t,gsl_matrix_get(C,j,0));
                  gsl_matrix_set(FH01,1,t,1);
                }
                
                if(t<a)
                {
                  for(v=a-1;v>=t;v--)
                  {
                    gsl_matrix_set(FH01,0,v+1,gsl_matrix_get(FH01,0,v));
                    gsl_matrix_set(FH01,1,v+1,gsl_matrix_get(FH01,1,v));
                  }
                  gsl_matrix_set(FH01,0,t,gsl_matrix_get(C,j,0));
                  gsl_matrix_set(FH01,1,t,1);
                }
                
                a=a+1;
              }
            }
            
            if(a==0)
            {
              gsl_matrix_set(FH01,0,0,gsl_matrix_get(C,j,0));
              gsl_matrix_set(FH01,1,0,1);
              a=a+1;
            }
          }
        }
        
        
        if(a==0) 
        {
          Rprintf("No failure time information for risk 1; Program exits\n");
        }
        
        
        gsl_matrix * H01 = gsl_matrix_alloc(3,a);                 /* baseline hazard function */
        
        for(i=0;i<3;i++)
        {
          if(i<=1)
          {
            for(j=0;j<a;j++)    gsl_matrix_set(H01,i,j, gsl_matrix_get(FH01,i,j));
          }
          if(i==2)   
          {
            for(j=0;j<a;j++)    gsl_matrix_set(H01,i,j,0.0001);
          }
        }
        
        
        gsl_matrix * preH01 = gsl_matrix_alloc(3,a);               
        
        iter=0;
        do
        {
          
          iter=iter+1;
          
          /* store previous parameter estimates */
          
          gsl_matrix_memcpy(prebeta0,beta0);
          gsl_vector_memcpy(prebeta1,beta1);
          gsl_vector_memcpy(presigma,sigma);
          gsl_vector_memcpy(pretheta,theta);
          gsl_vector_memcpy(presigmad,sigmad);
          gsl_vector_memcpy(preeta,eta);
          gsl_matrix_memcpy(prebtheta,btheta);
          gsl_matrix_memcpy(preH01,H01);
          
          pregamma=gamma;
          
          
          /* run EM to get updated parameter estimates */
          
          status=EM(beta0,beta1,sigma,theta,sigmad,eta,btheta,H01,&gamma,pdim,
                    Bio,C,M,bw,sigmau_inv, quadpoint, xs, ws, survvar);
          
          
          if (trace) {
            Rprintf("status=%d  iter=%d\n",status,iter);
            
            Rprintf("beta0 = \n");
            for (j=0;j<j_max;j++)
            {
              for(i=0;i<p_max;i++)  Rprintf("%f     ", gsl_matrix_get(beta0,j,i));
              Rprintf("\n");                     
            }
            
            
            Rprintf("beta1 = \n");
            for (i=0;i<beta1->size;i++)
            {
              Rprintf("%f     ", gsl_vector_get(beta1,i)); 
              
            }
            Rprintf("\n");
            
            
            Rprintf("sigma = \n");
            for (i=0;i<sigma->size;i++)
            {
              Rprintf("%f     ", gsl_vector_get(sigma,i)); 
              
            }
            Rprintf("\n");
            
            
            Rprintf("theta = \n");
            for (i=0;i<theta->size;i++)
            {
              Rprintf("%f     ", gsl_vector_get(theta,i)); 
              
            }
            Rprintf("\n");
            
            Rprintf("sigmad = \n");
            for (i=0;i<sigmad->size;i++)
            {
              Rprintf("%f     ", gsl_vector_get(sigmad,i)); 
              
            }
            Rprintf("\n");
            
            Rprintf("eta = \n");
            for (i=0;i<eta->size;i++)
            {
              Rprintf("%f     ", gsl_vector_get(eta,i)); 
              
            }
            Rprintf("\n");
            
            Rprintf("btheta = \n");
            for(k=0;k<q_b;k++)
            {
              for(i=0;i<k_max;i++)   Rprintf("%f   ",gsl_matrix_get(btheta,k,i));
              Rprintf("\n");
            }
            
            
            Rprintf("gamma = %f\n",gamma);
            }
          
        }while(Diffn(prebeta0,beta0,prebeta1,beta1,presigma,sigma,pretheta,
                     theta,preeta,eta,pregamma,gamma,presigmad,sigmad,conversigmad)==1 && status != 100 && 
                     iter<maxiter);
        
        if(status==100) 
        {
          printf("program stops because of error\n");
        }
        
        
        if(iter==maxiter) 
        {
          printf("program stops because of nonconvergence\n");
        }
        
        /** output the estimates  ***/
        NumericMatrix Bio_matrix(Bio->size1, Bio->size2);
        NumericMatrix C_matrix(C->size1, C->size2);
        NumericVector m_vec(M->size);
        NumericMatrix beta0_matrix(j_max, p_max);
        NumericVector beta1_estimate(beta1->size);
        NumericVector sigma2_estimate(sigma->size);
        NumericVector theta_estimate(theta->size);
        NumericVector sigmad_estimate(sigmad->size);
        NumericVector eta_estimate(eta->size);
        NumericMatrix btheta_matrix(q_b, k_max);
        NumericMatrix BaselineHazard(H01->size1, H01->size2);
        double loglike = 0.0;
        
        gsl_matrix *FUNA=gsl_matrix_alloc(M->size,k_max);
        NumericMatrix FUNA_matrix(FUNA->size1, FUNA->size2);
        
        for (i=0;i<Bio->size1;i++) 
        {
          for (j=0;j<Bio->size2;j++) 
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
        
        if(status==1 && iter<maxiter)
        {
          
          gsl_matrix_set_zero(V);
          for(k=0;k<k_max;k++) gsl_matrix_set(V,k,k,gsl_vector_get(sigmad,k));
          
          TransM(btheta,bthetat);
          MulMM(btheta,V,hbtheta);
          MulMM(hbtheta,bthetat,reta);
          
          gsl_linalg_SV_decomp(reta,veta,S,work);
          
          for(k=0;k<k_max;k++)
          {
            gsl_vector_set(sigmad,k,gsl_vector_get(S,k));
            
            for(i=0;i<q_b;i++)
              gsl_matrix_set(btheta,i,k,gsl_matrix_get(reta,i,k));
          }
          
          for (i=0;i<j_max;i++) 
          {
            for (j=0;j<p_max;j++) 
            {
              beta0_matrix(i, j) = gsl_matrix_get(beta0, i, j);
            }
          }
          
          for (i=0;i<beta1->size;i++) beta1_estimate(i) = gsl_vector_get(beta1, i);
          for (i=0;i<sigma->size;i++) sigma2_estimate(i) = gsl_vector_get(sigma, i);
          for (i=0;i<theta->size;i++) theta_estimate(i) = gsl_vector_get(theta, i);
          for (i=0;i<sigmad->size;i++) sigmad_estimate(i) = gsl_vector_get(sigmad, i);
          for (i=0;i<eta->size;i++) eta_estimate(i) = gsl_vector_get(eta, i);
          
          for (i=0;i<q_b;i++) 
          {
            for (j=0;j<k_max;j++) 
            {
              btheta_matrix(i, j) = gsl_matrix_get(btheta, i, j);
            }
          }
          
          loglike = GetLogLike(beta0,beta1,theta,sigmad,eta,btheta,H01,sigma,gamma,
                               pdim,Bio,C,M,bw,sigmau_inv, quadpoint, xs, ws);
          
          for (i=0;i<H01->size1;i++) 
          {
            for (j=0;j<H01->size2;j++) 
            {
              BaselineHazard(i, j) = gsl_matrix_get(H01, i, j);
            }
          }
          
          
          gsl_matrix *FUNA2=gsl_matrix_alloc(M->size,((k_max+1)*k_max)/2);
          gsl_matrix *FUNE=gsl_matrix_alloc(M->size,H01->size2);
          gsl_matrix *FUNAE=gsl_matrix_alloc(M->size,H01->size2);
          gsl_matrix *FUNA2E=gsl_matrix_alloc(M->size,H01->size2);
          gsl_matrix *FUNAEV=gsl_matrix_alloc(M->size,H01->size2*k_max);
          gsl_matrix *FUNA2EV=gsl_matrix_alloc(M->size,H01->size2*k_max);
          
          
          GetE(FUNA,FUNA2,FUNE,FUNAE,FUNA2E,FUNAEV,FUNA2EV,beta0,beta1,sigma,
               theta,sigmad,eta,btheta,H01,gamma,pdim,Bio,C,M,bw,sigmau_inv,quadpoint, xs, ws);  
          
          
          
          for (i=0;i<FUNA->size1;i++) 
          {
            for (j=0;j<FUNA->size2;j++) 
            {
              FUNA_matrix(i, j) = gsl_matrix_get(FUNA, i, j);
            }
          }
          
          
          gsl_matrix_free(FUNA2); 
          gsl_matrix_free(FUNE); 
          gsl_matrix_free(FUNAE);
          gsl_matrix_free(FUNA2E);
          gsl_matrix_free(FUNAEV);
          gsl_matrix_free(FUNA2EV);
          
          
        }
        
        int TotalPara = 0;
        TotalPara = (beta0->size1)*(beta0->size2) + beta1->size + sigma->size +
          theta->size + sigmad->size + eta->size + (btheta->size1)*(btheta->size2);
        
        double AIC = 0.0;
        
        AIC = 2*(TotalPara - loglike);
        
        Rcpp::List ret;
        ret["ydata"] = Bio_matrix;
        ret["cdata"] = C_matrix;
        ret["N1"] = m_vec;
        ret["beta0_matrix"] = beta0_matrix;
        ret["beta1_estimate"] = beta1_estimate;
        ret["sigma2_estimate"] = sigma2_estimate;
        ret["theta_estimate"] = theta_estimate;
        ret["sigmad_estimate"] = sigmad_estimate;
        ret["eta_estimate"] = eta_estimate;
        ret["btheta_matrix"] = btheta_matrix;
        ret["gamma"] = gamma;
        ret["BaselineHazard"] = BaselineHazard;
        ret["NumEventTime"] = a;
        ret["FUNA_matrix"] = FUNA_matrix;
        ret["iter"] = iter;
        ret["loglike"] = loglike;
        ret["AIC"] = AIC;
        ret["TotalPara"] = TotalPara;
        ret["SurvivalVar"] = survvar;
        
        
        gsl_matrix_free(FUNA);
        
        gsl_matrix_free(FH01);
        gsl_matrix_free(H01);
        gsl_matrix_free(preH01);
        
        gsl_matrix_free(Bio);
        gsl_matrix_free(C);
        gsl_vector_free(M);
        
        gsl_matrix_free(beta0);
        gsl_vector_free(beta1);
        gsl_vector_free(sigma);
        gsl_vector_free(theta);
        gsl_vector_free(sigmad);
        gsl_vector_free(eta);
        gsl_matrix_free(btheta);
        
        gsl_matrix_free(prebeta0);
        gsl_vector_free(prebeta1);
        gsl_vector_free(presigma);
        gsl_vector_free(pretheta);
        gsl_vector_free(presigmad);
        gsl_vector_free(preeta);
        gsl_matrix_free(prebtheta);
        
        gsl_vector_free(B_spl);
        gsl_vector_free(pdim);  
        gsl_bspline_free(bw);
        gsl_matrix_free(sigmau_inv);
        
        gsl_matrix_free(V);
        gsl_vector_free(S);
        gsl_vector_free(work);
        gsl_matrix_free(bthetat);
        gsl_matrix_free(hbtheta);
        gsl_matrix_free(reta);
        gsl_matrix_free(veta);

        return ret;
  }

}