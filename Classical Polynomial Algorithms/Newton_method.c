#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<complex.h>
#include<math.h>

#define B 0.177124344467704 //B = (3-sqrt(7))/2
#define C 0.219223593595584 //C = (5-sqrt(17))/4

int d;
double complex z, Dp;

double Gamma(double repol[], double impol[]);

int main(){
	int i, j, cont, b, t=-1;
	double g, max, re, im, radius, theta, R=0, e_out=0.0000000001;
	double complex D, p;
	srand(time(NULL));

	printf("\nThis program calculates the zeros of a complex monic polynomial\np(z) = a_0 + a_1*z + ... + a_(d-1)*z^(d-1) + z^d.\n\nEnter the value of the degree of the polynomial.\nd = ");
	scanf("%d",&d);
	double complex pol[d+1], sol[d];
	double gam[d], repol[d+1], impol[d+1];
	for(i=0; i<d; i++)	{
		sol[i] = 0;
		gam[i] = 0;
	}
	printf("Enter the values of the coefficients below. You must enter two real numbers separated by a space. The first number is the real part of the coefficient and the second is the imaginary part.\n");
	for(i=0; i<d; i++){ 
		printf("a_(%d) = ", i);
		scanf("%lf %lf",&re, &im);
                pol[i] = re + im*I;
		repol[i] = re;
		impol[i] = im;
	}
	pol[d] = 1;
	repol[d] = 1;
	impol[d] = 0;
	printf("\nYour polynomial is\n%lf + %lfi + ", creal(pol[0]), cimag(pol[0]));
	for(i=1; i<d; i++){ 
		printf("(%lf + %lfi)z^%d + ", creal(pol[i]), cimag(pol[i]), i);
	}
        printf("z^%d.\n\n",d); //Up to this point, we have the description of the program and the creation of the polynomial by the user.

	for(i=0; i<d; i++){
		R = R + cabs(pol[i]);
	}
	if(R < 1){
		R = 1;
	}
	printf("All the zeros of this polynomial are in the complex disc centered on the origin with radius %lf.\n", R); //The program determines the region where the zeros will be searched, since all zeros are contained in the disk centered on the origin with radius R = max{1, |a_0| + |a_1| + ... + |a_(d-1)|}. 
		
	printf("Computed zeros:\n\n");	

	for(i=0; i<1000000*R; i++){
		re = R*(float)rand()/RAND_MAX + rand()/RAND_MAX;
       		im = R*(float)rand()/RAND_MAX + rand()/RAND_MAX;
       		z = re + im*I; //z is a guess to start Newton's iteration. If z is an approximate zero, it will converge to the desired precision in a few iterations. The idea is to make a lot of random guesses and apply the iterations to each of them.
		
		b = 1;
		for(j=0; j<=t; j++){
			if(cabs(z-sol[t]) <= B/gam[t]){
				b = 0;
				break;
			}
		} //Checks if z is an approximate zero of some already computed zero. In the positive case, we do not do Newton's iteration and generate another point.

		if(b){
			cont = 0;
			Dp = 0;
			p = pol[0];
			for(j=1; j<=d; j++){
				Dp = Dp + j*pol[j]*cpow(z,j-1);
				p = p + pol[j]*cpow(z,j);
			}
			while(cabs(p) > e_out && cont < 10){
				if(Dp != 0){
					z = z - cpow(Dp,-1)*p;
					if(cabs(z)>R){
						break;
					}
				}
				else{
					break;
				}
				Dp = 0;
				p = pol[0];
				for(j=1; j<=d; j++){
					Dp = Dp + j*pol[j]*cpow(z,j-1);
					p = p + pol[j]*cpow(z,j);
				}
				cont++;			
			} //Newton's iteration starting from z. At most 100 iterations are done (there is nothing so special in this value).
			
			if(cabs(p) <= e_out){	
				 //Computes gamma (p, z), where z is the candidate to the computed solution by Newton iteration.
				g = Gamma(repol, impol);    		
				for(j=0; j<=t; j++){	
					if(gam[j] > g ){
						max = gam[j];
					}
					else{
						max = g;
					}
					//It verifies if the computed solution is very close to another one already computed. If the difference is less than this bound, then it is the same solution.
					if(cabs(z-sol[j]) < C/max){	
						b = 0;
						break;
					}
				}
				if(b){
					t++;
					sol[t] = z;								
					gam[t] = g; 
					if(cimag(z) < 0){
						printf("%.16lf %.16lfi, Gamma = %.16lf\n", creal(z), cimag(z), gam[t]);
					}
					else{
						printf("%.16lf +%.16lfi, Gamma = %.16lf\n", creal(z), cimag(z), gam[t]);
					}
				}
			} //Checks if z is a candidate for zero. If so, we test if z is too close to any zero already computed. If this happens, we consider it to be the same zero and discard it. Otherwise, we include z in the solution vector.
		}

		if(t == d-1){
			break;
		} //Checks whether the solution vector already contains all zeros.
	}

	printf("\nTotal of computed zeros: %d\nTotal of tested points: %d\n\n",t+1, i);
		
	return 0;
}

double Gamma(double repol[], double impol[]){	
	int i, j, k;
	double aux, re, im, g = 0;
	double complex D, DDp[d+1];
	
	for(i=0; i<=d-1; i++){	
		re = repol[i+1];
		im = impol[i+1];
		DDp[i] = (i+1)*(re+im*I);
	}
	DDp[d] = 0; //DDp[i] is the vector of the coordinates of the first derivative of p.
		
	for(j=2; j<=d; j++){
		D = 0;
		for(i=0; i<=d-j; i++){
			DDp[i] = (i+1)*DDp[i+1];
			D = D + DDp[i]*cpow(z,i);
		} //For each j, D is the jth derivative of p in z.
				
		k=1;
		for(i=2; i<=j; i++){
			k=k*i;
		}
		aux = pow( cabs(cpow(Dp,-1)*D/k), 1.0/(j-1.0) ); //Candidate to gamma(p,z).
		if(g < aux){
			g = aux;
		}
	}
	
	return g;
}
