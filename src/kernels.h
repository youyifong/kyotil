#ifndef KERNEL_h
	#define KERNEL_h
	typedef enum {LINEAR = 0,EUCLIDEAN = 1,POLYNOMIAL = 2,RBF = 3,IBS = 4,HAMMING=5}KERNEL_TYPE;
	typedef struct kernel_param_struct{
		KERNEL_TYPE type;
		double param;
		int n;          // length of 'params', 'ncx' for IBS kernel, 
		double* params; // length 'ncx' weights for IBS kernel
	}KERNEL_PARAM;
	typedef struct kernel_arg_struct{
		int nr1;
		int nc1;
		double* x1;
		int nr2;
		int nc2;
		double* x2;
	}KERNEL_ARG;
#endif	

void getK(KERNEL_ARG*,KERNEL_PARAM*,double* K);
void getKernel(int* _nr1,int* nc1,double* x1,int* _nr2,int* nc2,double* x2,int* kernel,double* para,double* K);
void edist2(int* _nr1,int* nc1,double* x1,int* _nr2,int* nc2,double* x2,double* dist);
void ibs2_kernel(int* _nrx, int* _ncx,double* x,int* _nry, int* _ncy,double* y,double* weights,double* K);
void ibs2X_kernel(int* _nrx, int* _ncx,double* x,double* weights,double* K);
void hammingSim_kernel(int* _nrx, int* _ncx,double* x,int* _nry, int* _ncy,double* y,double* weights,double* K);

