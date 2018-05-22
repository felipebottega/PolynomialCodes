#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "lapacke.h"

/* Global variables */
#define NUMBER_OF_TESTS 100
#define n 1500
	
/* Auxiliary routines prototypes  */
void print_vector( char* desc, lapack_int k, lapack_complex_double* sol );
double normal_distribution( void );
	
/* Main program */
int main(){
	srand(time(NULL));
	
        /* Local variables */
	lapack_int i, j, t, lda, ldvl, ldvr, lwork, ilo, ihi, info;
	char balanc, jobvl, jobvr, sense;
	double re, im, abnrm;
	balanc = 'B';
	jobvl = 'N';
	jobvr = 'N';
	sense = 'N';

	/* Polynomial degree verification */
	if( n < 1 ){ printf("Error: n must be bigger than 0.\n"); exit( 1 ); };
	lda = n;
	ldvl = n;
	ldvr = n;

	 /* Arrays */
	lapack_complex_double pol[n+1], a[n*n], w[n], vl[n], vr[n], work[1], sol[n*NUMBER_OF_TESTS];
	double scale[n], rconde[n], rcondv[n], rwork[2*n-1];
	
	for(t = 1; t <= NUMBER_OF_TESTS; t++){
		/* Creation of random polynomial and its companion matrix */
		for(i = 0; i < n; i++) { 
			re = normal_distribution();
			im = normal_distribution();
                	pol[i] = re + im*I;
		}	
		pol[n] = 1;
		for(i = 0; i < n; i++) { for(j = 0; j < n; j++) a[i*n+j] = 0; };
		for(i = 1; i < n; i++) a[i*n+i-1] = 1;
		for(i = 0; i < n; i++) a[(i+1)*n-1] = -pol[i]; 
		
		/* Computation of the companion matrix's eigenvalues */
		info = LAPACKE_zgeevx(LAPACK_COL_MAJOR, balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv);

		/* Vector with all computed zeros */
		if( info != 0 ) t--; 
		else{ for(i = 0; i < n; i++) sol[(t-1)*n+i] = w[i]; }
	
	}

	/* Print the computed zeros */	
	
        exit( 0 );
} 

/* Auxiliary routine: printing a vector */
void print_vector( char* desc, lapack_int k, lapack_complex_double* sol ) {
	lapack_int i;
	printf( "\n %s\n", desc );
	        for(i = 0; i < k; i++){
		    if( cimag(sol[i]) >= 0 ) printf( " %6.4lf+%6.4lfi\n", creal(sol[i]), cimag(sol[i]) );
			else printf( " %6.4lf%6.4lfi\n", creal(sol[i]), cimag(sol[i]) );
	        }
	        printf( "\n" );
}

/* Generates a pseudo-random number with distribution N(0,1) */
double normal_distribution( void ){
	lapack_int i, s;
	double mean, sum = 0;
		
	/* Sum of 1000 values with uniform distribution in [-1,1] */
	for(i = 1; i <= 1000; i++){
		s = rand();
		sum = sum + 2.0*((double)s/(double)(RAND_MAX+1.0))-1.0;
	}
	/* Computes the arithmetic mean of these 1000 random values. This mean is a random value with distribution N(0,1). */
	mean = sum/1000;

	return mean;
}
