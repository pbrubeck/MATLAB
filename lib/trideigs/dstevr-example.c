#include <stdlib.h>
#include <stdio.h>

/* DSTEVR prototype */
extern void dstevr(char* JOBZ, char *RANGE, int* N, double* D, double* E, double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z, int* LDZ, int* ISUPPZ, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);
/* Auxiliary routines prototypes */
extern void printMatrix( char* desc, int m, int n, double* a, int lda );

/* Main program */
int main() {
	/* Locals */
	int N, i;
	printf("Enter matrix size N=");
	scanf("%d",&N);
	double *D=(double*)malloc(N*sizeof(double));
	double *E=(double*)malloc(N*sizeof(double));
	for(i=0; i<N; i++){
		int k=i+1;
		D[i]=0;
		E[i]=k/sqrt(4*k*k-1);
	}
	
	double VL=0.0, VU=0.0, ABSTOL=0.0;
	int IL=0, IU=0, M=N, LDZ=N, LWORK=20*N, LIWORK=10*N, INFO;
	
	int *ISUPPZ=(int*)malloc(2*M*sizeof(int));
	int *IWORK=(int*)malloc(LIWORK*sizeof(int));
	double *W=(double*)malloc(N*sizeof(double));
	double *Z=(double*)malloc(LDZ*M*sizeof(double));
	double *WORK=(double*)malloc(LWORK*sizeof(double));
	
	dstevr("V","A",&N,D,E,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z,&LDZ,ISUPPZ,WORK,&LWORK,IWORK,&LIWORK,&INFO);
	
	/* Check for convergence */
	if(INFO>0){
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
	/* Print eigenvalues and eigenvectors*/
	printMatrix("Eigenvalues", 1, M, W, 1);
	printMatrix("Eigenvectors (stored columnwise)", LDZ, M, Z, LDZ);
	
	/* Free workspace */
	free((void*)WORK);
	free((void*)IWORK);
	exit(0);
} /* End of DSTEVR Example */

/* Auxiliary routine: printing a matrix */
void printMatrix(char* desc, int m, int n, double* a, int lda){
	int i, j;
	printf("\n %s\n", desc);
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			printf(" % .4f", a[i+j*lda]);
		}
		printf("\n");
	}
}
