// Author Krisztian Sebestyen ksebestyen@gmail.com
#include "kernels.h"
#include "matrix.h"
#include <strings.h>
#include <math.h>
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


// IBS2 kernel
// input matrix x in {0,1,2}
// x_{nr,nc} in column-major order, type double in {0,1,2}
// weights NULL or length nc
// K_{nr,nr} accessed in column-major order since only
// upper-triangular region is computed and copied to lower-triangular region
	 
void ibs2X_kernel(
int* _nrx, int* _ncx,double* x,
double* weights,double* K)
{
	int nr = *_nrx;
	int nc = *_ncx;

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
			
}


// x,y input matrices of type double containing values in {0,1,2}
// x_{nrx,ncx} in column-major order
// y_{nry,ncy} in column-major order
// weights NULL or length MIN(ncx,ncy)
// K_{nrx,nry} accessed in column-major 
// if 'x' and 'y' have differing column dimensions, the first 1..min(ncol(x),ncol(y)) are used


void ibs2_kernel(
int* _nrx, int* _ncx,double* x,
int* _nry, int* _ncy,double* y,
double* weights,        
double* K){
	int nrx = *_nrx;
	int nry = *_nry;
	int ncx = *_ncx;
	int ncy = *_ncy;

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
}






// Hamming *similarity* measure for binary data {0,1} (as opposed to hamming dissimilarity) 

void hammingSim_kernel(int* _nrx, int* _ncx,double* x,int* _nry,int* _ncy,double* y,double* weights,double* K)
{
	int nrx = *_nrx;
	int nry = *_nry;
	int ncx = *_ncx;
	int ncy = *_ncy;

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
	
}


// (euclidean distance)^2
// dist <- matrix(rowSums((X1[rep(1:nrow(X1),nrow(X2)),,drop=F] - X2[rep(1:nrow(X2),each=nrow(X1)),,drop=F])^2),nrow = nrow(X))
void edist2(int* _nr1,int* _nc1,double* x1,int* _nr2,int* _nc2,double* x2,double* dist)
{
	int nr1 = *_nr1;
	int nr2 = *_nr2;
	int nc1 = *_nc1;
	int nc2 = *_nc2;	

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
}




void getKernel(int* nr1,int* nc1,double* x1,int* nr2,int* nc2,double* x2,int* _kernel,double* para,double* K){


	KERNEL_TYPE kernel = (KERNEL_TYPE)*_kernel;
	if(kernel == LINEAR){
		tcrossprod(x1,nr1,nc1,x2,nr2,nc2,K);	
		return;
	}else if(kernel == EUCLIDEAN){
		edist2(nr1,nc1,x1,nr2,nc2,x2,K);
		return;
	}else if(kernel == POLYNOMIAL){
		tcrossprod(x1,nr1,nc1,x2,nr2,nc2,K);	
		for(int i = 0;i < (*nr1) * (*nr2);i++) 
			K[i] = pow(K[i] + 1.0,para[0]);
		return;
	}else if(kernel == RBF){
		edist2(nr1,nc1,x1,nr2,nc2,x2,K);
		for(int i = 0;i < (*nr1) * (*nr2);i++) 
			K[i] = exp(-para[0] * K[i]);
		for(int i = 0;i < (*nr1) * (*nr2);i++)
			if(EQ(K[i],0.0,DBL_EPSILON))
				K[i] = 0.0;			
		return;		
	}else if(kernel == IBS){
		ibs2_kernel(nr1,nc1,x1,nr2,nc2,x2,para,K);	
		return;
	}else if(kernel == HAMMING){
		hammingSim_kernel(nr1,nc1,x1,nr2,nc2,x2,para,K);	
		return;
	}
}


// getK is a more compact interface called from C code
// all kernel functions should be based on KERNEL_ARG and KERNEL_PARAM so that
// kernels with multi-dimensional parameters can be conveniently added without
// having to change all other C-code which use on getK or getKernel
void getK(KERNEL_ARG* data,KERNEL_PARAM* param,double* K){
	getKernel(&(data->nr1),&(data->nc1),data->x1,&(data->nr2),&(data->nc2),data->x2,
	(int*)&(param->type),&(param->param),K);
	return;
}

// for(i2 in 1:nrow(X2))
	// for(i1 in 1:nrow(X))
		// cat(i1," ",i2,"\n")

//rowSums((X[rep(1:nrow(X),nrow(X)),,drop=F] - X[rep(1:nrow(X),each=nrow(X)),,drop=F])^2)
// getK=function(X,kernel,para=NULL,X2=NULL){
    // kernel=substr(kernel,1,1)
    // if (kernel=="r" | kernel=="e") {
        // if (!is.null(X2)) {
            // aux = X[rep(1:nrow(X),nrow(X2)),,drop=F] - X2[rep(1:nrow(X2),each=nrow(X)),,drop=F]
            // dist.mat = matrix(rowSums(aux^2), nrow=nrow(X))
// #            aux=X2[rep(1:nrow(X2),nrow(X)),] - X[rep(1:nrow(X),each=nrow(X2)),]
// #            dist.mat = matrix(rowSums(aux^2), nrow=nrow(X2))
        // } else {
            // dist.mat = as.matrix(dist(X))^2
        // }
    // }
    
    // if (is.null(X2)) X2=X
    // switch(kernel, 
        // p=(tcrossprod(X,X2)+1)^para, # polynomial
        // r=exp(-para*dist.mat), # rbf
        // e=dist.mat, # Euclidean edist
        // l=tcrossprod(X,X2), # linear
        // i = ibs(X,X2), # IBS 
		// stop(kernel %+% " kernel not supported")
    // )
// }


