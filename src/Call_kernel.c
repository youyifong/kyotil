// Author: Krisztian Sebestyen ksebestyen@gmail.com

#include "kernels.h"
#include "matrix.h"
#include <Rinternals.h>
#include <strings.h>
#include <float.h>



#define	ABS(A)  	((A) > 0 ? (A) : (-(A)))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))
#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define EQ(x,y,eps)  ((ABS((x)-(y))) <= (eps) ? 1 : 0) 
#define CX(i,j,ld) (((j)*(ld))+(i)) //column-major index
#define RX(i,j,ld) (((i)*(ld))+(j)) //row-major index
#define I(x,y) ((x) == (y))
#define IBS2(x,y) (0.5*(2.0 - (ABS((x)-(y)))))
#define IBS(x,y,order) (((order) - (ABS((x)-(y)))) / (order))
// IBS2(x,y)   = (2.0 - |x - y|)/2.0, x,y in {0,1,2} for a measure
// K_IBS2(i,j) = 0.5*sum(w[] * (2 - |x[i,]-y[j,]|)) / sum(w)
// Hamming similarity measure
// I(x,y) = IBS1(x,y) = sum(w * (1 - |x-y|)) / sum(w) for x,y in {0,1}

// (euclidean distance)^2
// dist <- matrix(rowSums((X1[rep(1:nrow(X1),nrow(X2)),,drop=F] - X2[rep(1:nrow(X2),each=nrow(X1)),,drop=F])^2),nrow = nrow(X))
SEXP Call_edist2(SEXP _x1, SEXP _x2,SEXP _dist)
{
     int nr1=nrows(_x1), nc1=ncols(_x1), nr2=nrows(_x2), nc2=ncols(_x2);
     double *x1=REAL(_x1), *x2=REAL(_x2), *dist=REAL(_dist);

	double s,ss;
	int i1,i2,j;
	int nc = MIN(nc1,nc2);	
	memset(dist,0,(size_t)(nr1 * nr2 * sizeof(double)));
	for(i2 = 0;i2 < nr2;i2++){
		for(i1 = 0;i1 < nr1;i1++){
			ss = 0.0;
			for(j = 0;j < nc;j++){
				s = x1[CX(i1,j,nr1)] - x2[CX(i2,j,nr2)];
				s *=s;
				ss += s;
			}
			dist[CX(i1,i2,nr1)] = ss;		
		}
	}	
	
	return R_NilValue;
}

SEXP Call_hammingSim_kernel(SEXP _x,SEXP _y,SEXP _weights,SEXP _K)
{
	int nrx=nrows(_x), ncx=ncols(_x), nry = nrows(_y), ncy = nrows(_y);
	double *x=REAL(_x), *y=REAL(_y), *K=REAL(_K);
	double* weights = isReal(_weights) ? REAL(_weights) : NULL;

	int i,j,k;
	int nc = MIN(ncx,ncy);
	double sum_w;
	if(!weights){
		sum_w = (double)nc;
		for(i = 0;i < nrx;i++){
			for(j = 0;j < nry;j++){
				int s = 0;
				for(k = 0;k < nc;k++) s += I(x[CX(i,k,nrx)],y[CX(j,k,nry)]);
				K[CX(i,j,nrx)] = (double)s / sum_w;
			}
		}
	}else{
		sum_w = 0.0;
		for(k = 0;k < nc;k++) sum_w += weights[k];
		for(i = 0;i < nrx;i++){
			for(j = 0;j < nry;j++){
				double s = 0.0;
				for(k = 0;k < nc;k++) s += weights[k] * (double)(I( x[CX(i,k,nrx)],y[CX(j,k,nry)]));
				K[CX(i,j,nrx)] = s / sum_w;
			}
		}	
	}	
	
	return R_NilValue;
}


SEXP Call_ibs2_kernel(SEXP _x,SEXP _y,SEXP _weights,SEXP _K){
	int nrx=nrows(_x), ncx=ncols(_x), nry = nrows(_y), ncy = nrows(_y);
	double *x=REAL(_x), *y=REAL(_y), *K=REAL(_K);
	double* weights = isReal(_weights) ? REAL(_weights) : NULL;

	int i,j,k;
	int nc = MIN(ncx,ncy);
	double s;
	double denom = (double)nc;		
	if(!weights){
		for(i = 0;i < nrx;i++){
			for(j = 0;j < nry;j++){
				s = 0.0;
				for(k = 0;k < nc;k++) s += IBS(x[CX(i,k,nrx)],y[CX(j,k,nry)],2);				
				K[CX(i,j,nrx)] = s;
			}
		}
	}else{
		double sum_w = 0.0;for(k = 0;k < nc;k++)sum_w += weights[k];
		denom = sum_w;
		for(i = 0;i < nrx;i++){
			for(j = 0;j < nry;j++){
				s = 0.0;				
				for(k = 0;k < nc;k++) s += weights[k] * IBS(x[CX(i,k,nrx)],y[CX(j,k,nry)],2);
				K[CX(i,j,nrx)] = s;
			}
		}	
	}
	for(i = 0;i < nrx;i++){
		for(j = 0;j < nry;j++){	
			K[CX(i,j,nrx)] /= denom;
		}
	}
	
	return R_NilValue;
}


SEXP Call_ibs2X_kernel(SEXP _x,SEXP _weights,SEXP _K)
{
	int nr=nrows(_x), nc=ncols(_x);
	double *x=REAL(_x), *K=REAL(_K);
	double* weights = isReal(_weights) ? REAL(_weights) : NULL;
	 
	int i,j,k;
	double s;
	double denom = (double)nc;
	if(!weights){
		for(i = 0;i < nr;i++){
			for(j = i;j < nr;j++){		
				s = 0.0;
				for(k = 0;k < nc;k++) s += IBS(x[CX(i,k,nr)],x[CX(j,k,nr)],2.0);
				K[CX(i,j,nr)] = s;
			}
		}
	}else{
		double sum_w = 0.0;for(k = 0;k < nc;k++)sum_w += weights[k];
		denom = sum_w;	
		for(i = 0;i < nr;i++){
			for(j = i;j < nr;j++){		
				s = 0.0;
				for(k = 0;k < nc;k++) s += weights[k] * IBS(x[CX(i,k,nr)],x[CX(j,k,nr)],2.0);
				K[CX(i,j,nr)] = s;
			}
		}	
	}
	for(i = 0;i < nr;i++){
		for(j = i;j < nr;j++){	
			K[CX(i,j,nr)] /= denom;
		}
	}
	for(i = 0;i < nr - 1;i++)
		for(j = i+1;j < nr;j++)
			K[CX(j,i,nr)] = K[CX(i,j,nr)];
			
	return R_NilValue;
}


SEXP Call_getKernel(SEXP _x1, SEXP _x2, SEXP _kernel, SEXP _para,SEXP _K) { 
    
	// Rprintf("x1(%p) x2(%p) kernel(%p) para(%p) K(%p)\n",_x1,_x2,_kernel,_para,_K);
    int nr1=nrows(_x1), nc1=ncols(_x1), nr2=nrows(_x2), nc2=ncols(_x2);
    double *x1=REAL(_x1), *x2=REAL(_x2), *K=REAL(_K);
    
	
	KERNEL_TYPE kernel = (KERNEL_TYPE)(*INTEGER(_kernel));
	if(kernel == LINEAR){
		tcrossprod(x1,&nr1,&nc1,x2,&nr2,&nc2,K);	
	}else if(kernel == EUCLIDEAN){
		Call_edist2(_x1,_x2,_K);
	}else if(kernel == POLYNOMIAL){
		double para = *REAL(_para); // do not allow 'NULL' parameter
		tcrossprod(x1,&nr1,&nc1,x2,&nr2,&nc2,K);	
		for(int i = 0;i < (nr1) * (nr2);i++) K[i] = pow(K[i] + 1.0,para);
	}else if(kernel == RBF){
		double para = *REAL(_para); // do not allow 'NULL' parameter
		Call_edist2(_x1,_x2,_K);
		for(int i = 0;i < nr1*nr2;i++) 
			K[i] = exp(-para * K[i]);
		for(int i = 0;i < (nr1) * (nr2);i++)
			if(EQ(K[i],0.0,DBL_EPSILON))
				K[i] = 0.0;			
	}else if(kernel == IBS){
		Call_ibs2_kernel(_x1,_x2,_para,_K);
	}else if(kernel == HAMMING){
		Call_hammingSim_kernel(_x1,_x2,_para,_K);
	}
    return R_NilValue;
}
