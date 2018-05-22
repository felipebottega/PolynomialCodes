#include <stdlib.h>
#include <stdio.h>
#include "lapacke.h"
	
/* Auxiliary routine prototype */
void print_pol( char* desc, lapack_int n, lapack_complex_double* pol);
void print_matrix( char* desc, lapack_int n, lapack_complex_double* a );
void print_vector(char* desc, lapack_int n, lapack_complex_double* w);
	
/* Main program */
int main() {
	printf("This program calculates the zeros of a complex monic polynomial p=z^n + a_(n-1)*z^(n-1) + ... + a_(1)*z + a_0 using the QR iteration.\n");
        /* Variáveis */
	lapack_int i, j, n, lda, ldvl, ldvr, lwork, ilo, ihi, info;
	char balanc, jobvl, jobvr, sense;
	double re, im, abnrm;
	balanc = 'B';
	jobvl = 'N';
	jobvr = 'N';
	sense = 'N';

	/* Definition of the degree of the polynomial */
	printf("Enter the value of the degree n of the polynomial.\n");
	printf("n = ");
	scanf("%d",&n);
	if( n < 1 ){ printf("Error: n must be greater than 0.\n"); exit( 1 ); };
	lda = n;
	ldvl = n;
	ldvr = n;

	 /* Arrays */
	lapack_complex_double pol[n+1], a[n*n], w[n], vl[n], vr[n], work[1];
	double scale[n], rconde[n], rcondv[n], rwork[2*n-1];

	/* Definition of the polynomial and the companion matrix */
	printf("\nEnter the values of the coefficients of p. You must enter two real numbers separated by a space. The first number is the real part of the coefficient and the second is the imaginary part.\n");
	for(i = 0; i < n; i++) { 
		printf("a_(%d) = ", i);
		scanf("%lf %lf",&re, &im);
                pol[i] = re + im*I;
	}	
	pol[n] = 1;
	for(i = 0; i < n; i++) { for(j = 0; j < n; j++) a[i*n+j] = 0; };
	for(i = 1; i < n; i++) a[i*n+i-1] = 1;
	for(i = 0; i < n; i++) a[(i+1)*n-1] = -pol[i]; 

	/* Print polynomial */
	print_pol("p = ", n, pol);
       
	/* Print companion matrix */
 	print_matrix( "A = ", n, a );

	/* Computation of the eigenvalues of the companion matrix */
	info = LAPACKE_zgeevx(LAPACK_COL_MAJOR, balanc, jobvl, jobvr, sense, n, a, lda, w, vl, ldvl, vr, ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv);

	/* Print solution */
	if( info != 0 ){ printf("The algorithm failed to compute the zeros.\n"); exit( 1 ); }
	else print_vector("Zeros of p:", n, w);

        exit( 0 );
} 

/* Auxiliary routine to print the polynomial */
void print_pol( char* desc, lapack_int n, lapack_complex_double* pol) {
	lapack_int i;
	printf("p = ");
	for(i = n; i >= 1; i--) printf( "(%6.4lf + %6.4lfi)z^%d + ", creal(pol[i]), cimag(pol[i]), i );
	printf( "%6.4lf + %6.4lfi\n", creal(pol[0]), cimag(pol[0]) );
}
	
/* Auxiliary routine to print the matrix */
void print_matrix( char* desc, lapack_int n, lapack_complex_double* a ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
			 if( i*n+j == (i+1)*n-1 && cimag(a[i*n+j] ) >= 0) printf( "     %6.4lf + %6.4lfi", creal(a[i*n+j]), cimag(a[i*n+j]) );
			 else if( i*n+j == (i+1)*n-1 && cimag(a[i*n+j] ) < 0 ) printf( "     %6.4lf %6.4lfi", creal(a[i*n+j]), cimag(a[i*n+j]) );
			 else	printf( " %6.0lf", creal(a[i*n+j]) );
		}
                printf( "\n" );
        }
}

/* Auxiliary routine to print the vector */
	void print_vector( char* desc, lapack_int n, lapack_complex_double* w ) {
	        lapack_int i;
	        printf( "\n %s\n", desc );
	        for(i = 0; i < n; i++) {
			if( cimag(w[i]) >= 0 ) printf( " %6.4lf + %6.4lfi\n", creal(w[i]), cimag(w[i]) );
			else printf( " %6.4lf %6.4lfi\n", creal(w[i]), cimag(w[i]) );
		}
	        printf( "\n" );
	}
