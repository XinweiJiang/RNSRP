#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mex.h>
#include <math.h>
#include "matrix.h"

#define delta 1e-12

/*
 Euclidean Projection onto l1 Ball (eplb)
 
        min  1/2 ||x- y||_2^2
        s.t. ||x||_1 <= z
 
which is converted to the following zero finding problem
 
        f(lambda)= sum( max( |y|-lambda,0) )-z=0
 
 Usage:
 [x, lambda, iter_step]=eplb(y, n, z, lambda0);
 
 */

void eplb(double * x, double *root, int * steps, double * v,
int n, double z, double lambda0)
{
    
    int i, flag=0;
    int rho_1, rho_2, rho, rho_T, rho_S;
    int L_i, R_i, V_i;
    double lambda_1, lambda_2, lambda_T, lambda_S, lambda;
    double s_1, s_2, s, s_T, s_S, v_max, v_V_i;
    double f_lambda_1, f_lambda_2, f_lambda, f_lambda_T, f_lambda_S;
    int *V, *L, *R, *dp;
    int iter_step=0;
    
    V=(int *)malloc(sizeof(int)*n);
    L=(int *)malloc(sizeof(int)*n);
    R=(int *)malloc(sizeof(int)*n);
    
    if( (V==NULL) || (L==NULL) || (R==NULL) ){
        printf("\n Memory allocation failure!");
        exit (-1);
    }
    
    /* find the maximal absolute value in v
     * and copy the (absolute) values from v to x
     */
    s_1=x[0]=v_max=abs(v[0]);
    for (i=1;i<n; i++){
        if (v[i]>0)
            x[i]=v[i];
        else
            x[i]=-v[i];
        /* x[i]= fabs(v[i]); */
        
        if (x[i]> v_max)
            v_max=x[i];
                
        s_1+=x[i];
        /* s_1 = ||v||_1=||x||_1*/
    }
    
    /* If ||v||_1 <= z, then v is the solution  */
    if (s_1 <= z){
        flag=1;
        lambda=0;
    }
    
    if (!flag){
    /* The interval where the root lies  */
        lambda_1=v_max-z; lambda_2=v_max;
        
    /*-------------------------------------------------------------------
                               Initialization
     *-------------------------------------------------------------------
     */
        if ( (lambda0<lambda_2) && (lambda0> lambda_1) )
            lambda=lambda0;
        else
            lambda=lambda_1;
        
        R_i=L_i=0; s=0;
        for(i=0;i<n;i++){
            if (x[i] > lambda){
                R[ R_i ]= i; R_i++;
                s+=x[i];
            }
            else{
                L[ L_i ]=i; L_i++;
            }
        }
        rho=R_i;
        f_lambda=s-rho* lambda -z;
        
        /* lambda is the solution  */
        if ( fabs(f_lambda) <= delta )
            flag=1;
        
        if (f_lambda >0){
            dp=V;  V=R;  R=dp; V_i=R_i;
            
         /* update lambda_1 with lambda */
            lambda_1=lambda;
            rho_1=rho; s_1=s; f_lambda_1=f_lambda;
            
        /* compute the needed values for lambda_2 */
            rho_2=0; s_2=0; f_lambda_2=-z;
        }
        else{
            dp=V; V=L; L=dp; V_i=L_i;
            
        /* update lambda_2 with lambda */
            lambda_2=lambda;
            rho_2=rho; s_2=s; f_lambda_2=f_lambda;
            
        /* update the values for lambda_1 
         * s_1 has been computed before 
         */
            lambda_1=-delta; rho_1=n;
            f_lambda_1=s_1-rho_1* lambda_1 -z;
        }
     /*-------------------------------------------------------------------
                          End of initialization
      *--------------------------------------------------------------------
      */       
        
    }/* end of if(!flag) */
    
    while (!flag){
        iter_step++;
        
        /* compute lambda_T  */
        lambda_T=lambda_1 + f_lambda_1 /rho_1;
        if(rho_2 !=0){
            if (lambda_2 + f_lambda_2 /rho_2 >	lambda_T)
                lambda_T=lambda_2 + f_lambda_2 /rho_2;
        }
        
        /* compute lambda_S */
        lambda_S=lambda_2 - f_lambda_2 *(lambda_2-lambda_1)/(f_lambda_2-f_lambda_1);
        
        if (fabs(lambda_T-lambda_S) <= delta){
            lambda=lambda_T;
            break;
        }
        
        /* set lambda as the middle point of lambda_T and lambda_S */
        lambda=(lambda_T+lambda_S)/2;
        
        s_T=s_S=s=0;
        rho_T=rho_S=rho=0;
        R_i=L_i=0;
        
        for (i=0;i<V_i;i++){
            v_V_i=x[V[i]];
            
            if (v_V_i>lambda){
                if (v_V_i > lambda_S){
                    s_S+=v_V_i;
                    rho_S++;
                }
                else{
                    s+=v_V_i; R[R_i]=V[i]; R_i++;
                }
            }
            else{
                if (v_V_i > lambda_T){
                    s_T+=v_V_i; L[L_i]=V[i]; L_i++;
                }
            }
        }
        
        rho=R_i; rho_T=L_i;
        s_S+=s_2; rho_S+=rho_2;
        s+=s_S; rho+=rho_S;
        s_T+=s; rho_T+=rho;
        f_lambda_S=s_S-rho_S*lambda_S-z;
        f_lambda=s-rho*lambda-z;
        f_lambda_T=s_T-rho_T*lambda_T-z;
        
        if ( fabs(f_lambda)< delta ){
            flag=1;
            break;
        }
        if ( fabs(f_lambda_S)< delta ){
            lambda=lambda_S; flag=1;
            break;
        }
        if ( fabs(f_lambda_T)< delta ){
            lambda=lambda_T; flag=1;
            break;
        }
        
        /*
        printf("\n\n f_lambda_1=%5.6f, f_lambda_2=%5.6f, f_lambda=%5.6f",f_lambda_1,f_lambda_2, f_lambda);
        printf("\n lambda_1=%5.6f, lambda_2=%5.6f, lambda=%5.6f",lambda_1, lambda_2, lambda);
        printf("\n rho_1=%d, rho_2=%d, rho=%d ",rho_1, rho_2, rho);
         */
        
        if (f_lambda <0){
            lambda_2=lambda;  s_2=s;  rho_2=rho;
            f_lambda_2=f_lambda;            
            
            lambda_1=lambda_T; s_1=s_T; rho_1=rho_T;
            f_lambda_1=f_lambda_T;
            
            dp=V; V=L; L=dp; V_i=L_i;
        }
        else{
            lambda_1=lambda;  s_1=s; rho_1=rho;
            f_lambda_1=f_lambda;
            
            lambda_2=lambda_S; s_2=s_S; rho_2=rho_S;
            f_lambda_2=f_lambda_S;
            
            dp=V; V=R; R=dp; V_i=R_i;
        }
        
        if (V_i==0){
            //printf("\n rho=%d, rho_1=%d, rho_2=%d",rho, rho_1, rho_2);
            lambda=(s - z)/ rho;
            break;
        }
    }/* end of while */
    
    free(V); free(R);free(L);
    
    for(i=0;i<n;i++){        
        if (x[i] > lambda)
            x[i]-=lambda;
        else
            x[i]=0;
        
        if (v[i] < 0)
            x[i]=-x[i];
    }
    *root=lambda;
    *steps=iter_step;
}


void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* v=            mxGetPr(prhs[0]);
    int     n=            mxGetScalar(prhs[1]);
    double  z=            mxGetScalar(prhs[2]);
    double  lambda0=      mxGetScalar(prhs[3]);
    
    double *x, *lambda;
    int *iter_step;
    /* set up output arguments */
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[2] = mxCreateNumericMatrix(1,1, mxINT32_CLASS, 0);
    
    x=mxGetPr(plhs[0]);
    lambda=mxGetPr(plhs[1]);
    iter_step=mxGetPr(plhs[2]);
    
    eplb(x, lambda, iter_step, v, n, z, lambda0);
}

