#include <float.h> //DBL_EPSILON
#include <R_ext/Lapack.h>
#include <R.h>
#include <Rinternals.h>
#include "matrix.h"
//#include <R_ext/Applic.h>
//#include <R_ext/BLAS.h>
//#include <R_ext/RS.h> //R definitions for 'extending' R, registering functions,...

#define PRINTF Rprintf
#define MAX(A,B)    ((A) > (B) ? (A) : (B))
#define MIN(A,B)    ((A) < (B) ? (A) : (B))



//    sresid.var=array(NA, dim=c(nvar, nvar, ndead, ndead))
//    for (u in 1:ndead) for (m in u:ndead)   sresid.var[,,u,m] <- (u==m)*imat[,,u] - imat[,,u] %*% v.inv %*% imat[,,m]  
//    for (u in 1:(ndead-1)) for (m in (u+1):ndead)   sresid.var[,,m,u] <- t(sresid.var[,,u,m])
//    test.var=matrix(0,nvar,nvar)
//    for (u in 1:ndead) for (m in 1:ndead)   test.var=test.var+xx[u]*xx[m]*sresid.var[,,u,m]
            

SEXP compute_var(SEXP _nvar, SEXP _xx, SEXP _imat, SEXP _vinv){
    
    int p=asInteger(_nvar);
    int n=length(_xx);     
    double *xx=REAL(_xx), *imat=REAL(_imat), *vinv=REAL(_vinv);
    
    SEXP _ans=PROTECT(allocMatrix(REALSXP, p, p)); 
    double *test_var=REAL(_ans);

    //PRINTF("%f %f \n",dil_r,sigma_sq);
    //PRINTF("%i %i \n",n,p);
    
    double * tmp1 = (double *) malloc(p * p * sizeof(double));
    double * tmp   = (double *) malloc(p * p * sizeof(double));
    int i,j; // index 1:p
    for (i=0; i<p; i++) for (j=0; j<p; j++) test_var[i+j*p]=0;     

    double w; // weight
    int u,m; // index 1:n
    for (u=0;u<n;u++) for (m=u;m<n;m++) {
        w=xx[u]*xx[m];
        // imat[,,u] %*% v.inv %*% imat[,,m]  is stored in tmp
        matprod(imat+u*p*p, &p, &p,  vinv, &p, &p, tmp1);
        matprod(tmp1, &p, &p,  imat+m*p*p, &p, &p, tmp);
        // no need to create sresid.var, directly use sresid.var[,,u,m] to compute test_var
        for (i=0; i<p; i++) for (j=0; j<p; j++) {
            test_var[i+j*p] += w*(-tmp[i+j*p]);
            if (u==m) {
               test_var[i+j*p] += w*imat[i+j*p+u*p*p];
            } else {
               test_var[i+j*p] += w*(-tmp[j+i*p]);
            }
        }
    }
    
            
    free(tmp);
    free(tmp1);
    UNPROTECT(1);
    return _ans;
}

