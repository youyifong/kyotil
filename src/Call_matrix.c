// Author: Krisztian Sebestyen ksebestyen@gmail.com

#include "matrix.h"
#include <Rinternals.h>
#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))


// R's .Call interface
SEXP Call_upper_trap(SEXP _x,SEXP _diag, SEXP _k){
	int nr = nrows(_x);
	int nc = ncols(_x);
	// no need to protect since 'rrbind()' is native C 
    SEXP _ans = PROTECT(allocMatrix(REALSXP, nr, nc));
	memset(REAL(_ans),0,nr * nc * sizeof(double));
	upper_trap(nrows(_x) , ncols(_x),REAL(_x) , isReal(_diag) ? REAL(_diag) : NULL , *INTEGER(_k) , REAL(_ans));
	UNPROTECT(1);
	return _ans;
}
SEXP Call_lower_trap(SEXP _x,SEXP _diag, SEXP _k){
	int nr = nrows(_x);
	int nc = ncols(_x);
    SEXP _ans=PROTECT(allocMatrix(REALSXP, nr, nc));
	memset(REAL(_ans),0,nr * nc * sizeof(double));
	lower_trap(nrows(_x) , ncols(_x),REAL(_x) , isReal(_diag) ? REAL(_diag) : NULL , *INTEGER(_k) , REAL(_ans));
	UNPROTECT(1);
	return _ans;
}

SEXP Call_rrbind(SEXP _x,SEXP _times, SEXP _each,SEXP _vec_each){
	int nr = nrows(_x);
	int nc = ncols(_x);
	int nr_each = 0;
	int nry = nr;
	int times = isInteger(_times) ? *INTEGER(_times) : 0;
	int each = isInteger(_each) ? *INTEGER(_each) : 0;
	int* vec_each = (_vec_each == R_NilValue) ? NULL : INTEGER(_vec_each);
	
	if((nr == 0) || (nc == 0)) return R_NilValue;
	
	times = MAX(times,0);
	each = MAX(each,0);
	if((times == 0) && (each == 0) && (!vec_each)){
		Rprintf("rrbind: Error, both 'times' and 'each' are 0\n");
		return R_NilValue;
	}
	
// imitate r's 'rep' for a matrix	
// exactly one of times or each or vec_each can be executed in rrbind/rcbind C code	
// (otherwise an extra buffer would be needed because C code uses memcpy(), and
// memcpy() requires disjoint buffers which will fail for rrbind in column-major order when 'each' is done first
// )
	if((times > 0) && (each <= 1)) each = 0;
	if((each > 0) && (times <= 1)) times = 0;
	
	if(vec_each){
		each = 0;
		nr_each = 0;
		for(int i = 0;i < nr;i++)nr_each += vec_each[i];
	}else if(each > 0) nr_each = each * nr; 
	
	if(nr_each > 1) nry = nr_each;
	if(times > 1) nry *= times;

	// Rprintf("nrx(%i) ncx(%i) nr_each(%i) nry(%i) each(%i) times(%i)\n",nr,nc,nr_each,nry,each,times);
	
	SEXP _ans=PROTECT(allocMatrix(REALSXP, nry, nc));
	memset(REAL(_ans),0,nry * nc * sizeof(double));
	
	if((times == 0) || ((each == 0) && !vec_each)){
		rrbind(REAL(_x),nr,nc,times,each,vec_each,REAL(_ans));
		UNPROTECT(1);
		return _ans;
	}

	// do 'each' first, no need to protect since only native C is executed
    SEXP _ans_each=PROTECT(allocMatrix(REALSXP, nr_each, nc));
	memset(REAL(_ans_each),0,nr_each * nc * sizeof(double));
	
	// execute 'each' replication 
	rrbind(REAL(_x),nr,nc,0,each,vec_each,REAL(_ans_each));
	// execute 'times' replication 
	rrbind(REAL(_ans_each),nr_each,nc,times,0,NULL,REAL(_ans));
	
	UNPROTECT(2);
	return _ans;	
}
SEXP Call_rcbind(SEXP _x,SEXP _times, SEXP _each,SEXP _vec_each){
	int nr = nrows(_x);
	int nc = ncols(_x);
	int nc_each = 0;
	int ncy = nc;	
	int times = isInteger(_times) ? *INTEGER(_times) : 0;
	int each = isInteger(_each) ? *INTEGER(_each) : 0;
	int* vec_each = (_vec_each == R_NilValue) ? NULL : INTEGER(_vec_each);
  	double* x = REAL(_x);

	if((nr == 0) || (nc == 0)) return R_NilValue;
	
	
	times = MAX(times,0);
	each = MAX(each,0);
	if((times == 0) && (each == 0) && (!vec_each)){
		Rprintf("rcbind: Error, both 'times' and 'each' are 0\n");
		return R_NilValue;
	}
	
// imitate r's 'rep' for a matrix	
// exactly one of times or each or vec_each can be executed in rrbind/rcbind C code	
// (otherwise an extra buffer would be needed because C code uses memcpy(), and
// memcpy() requires disjoint buffers which will fail for rcbind in column-major order when 'each' is done first
// )
	if((times > 0) && (each <= 1)) each = 0;
	if((each > 0) && (times <= 1)) times = 0;
	
	if(vec_each){
		each = 0;
		nc_each = 0;
		for(int i = 0;i < nc;i++)nc_each += vec_each[i];
	}else if(each > 0) nc_each = each * nc; 
	 
	if(nc_each > 1) ncy = nc_each;
	if(times > 1) ncy *= times;
	
	SEXP _ans=PROTECT(allocMatrix(REALSXP, nr, ncy));
	double* ans = REAL(_ans);
	memset(ans,0,nr * ncy * sizeof(double));
	
	if((times == 0) || ((each == 0) && !vec_each)){
		rcbind(x,nr,nc,times,each,vec_each,ans);
		UNPROTECT(1);		
		return _ans;
	}

    SEXP _ans_each=PROTECT(allocMatrix(REALSXP, nr, nc_each));
	double* ans_each = REAL(_ans_each);
	memset(ans_each,0,nr * nc_each * sizeof(double));
	
	// execute 'each' replication 
	rcbind(x,nr,nc,0,each,vec_each,ans_each);
	// execute 'times' replication 
	rcbind(ans_each,nr,nc_each,times,0,NULL,ans);
		
	UNPROTECT(2);		
	return _ans;	
}


SEXP Call_dxd(SEXP _d1, SEXP _x, SEXP _d2){
	int n=length(_d1);
	SEXP _ans=PROTECT(allocMatrix(REALSXP, n, n));
	double *d1=REAL(_d1), *d2=REAL(_d2), *x=REAL(_x), *ans=REAL(_ans);

	int i,j;
	for(j = 0;j < n;j++){
		for(i = 0;i < n;i++){
			*ans = d1[i] * (*x) * d2[j];
			ans++; x++;
		}
	}
	
	UNPROTECT(1);
    return _ans;
}


SEXP xdx(SEXP _x, SEXP _d){
	int n = length(_d);
	int p = ncols(_x); // _x is n by p in R
	SEXP _ans=PROTECT(allocMatrix(REALSXP, p, p));
	double *d=REAL(_d), *x=REAL(_x), *ans=REAL(_ans);

	int i,j,k;
	for(j = 0;j < p;j++)
		for(i = 0;i < p;i++)
			ans[j+i*p]=0;
			
	for(k=0; k<n; k++)
    	for(j = 0;j < p;j++){
    		for(i = 0;i < p;i++){
                // row k col j of x in R is row j col k of x in C
    			ans[j+i*p] += d[k] * x[j*n+k] * x[i*n+k];
    		}
    	}
	
	UNPROTECT(1);
    return _ans;
}

