#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

// Mitchell Miller
// PHYS 440 Lab 1
// Exc. 2.1 1/19/12
// Hybrid Method

double func(double x){

	// Specify desired function
	double funcValue;
	funcValue = (1./128.)*(6435.*x*x*x*x*x*x*x*x-12012.*x*x*x*x*x*x+6930.*x*x*x*x-1260.*x*x+35.);
	return(funcValue);

}

double deriv(double x){

	// Provide the derivative of your function
	double derValue;
	derValue = (1./128.)*(6435.*8.*x*x*x*x*x*x*x-12012.*6.*x*x*x*x*x+6930.*4.*x*x*x-1260.*2.*x);
	return(derValue);

}

double hybrid(double left,double right){

	// Initialize variables
	double fLeft = func(left);
	double fRight = func(right);
	double best, fBest, dBest;
	double funcX, derX;
	double delta;
	int counter;
	double error = 1;
	double tol = 0.00000005;

	if(fabs(fLeft <= fabs(fRight))){

		best = left;
		fBest = fLeft;

	}
	else{

		best = right;
		fBest = fRight;

	}

	dBest = deriv(best);

	while(error > tol){

		counter++;
		if((dBest*(best-left)-fBest) * (dBest*(best-right)-fBest) <= 0){

			delta = -1.*(fBest) / (dBest);
			best = best + delta;

		}
		else{

			delta = (right - left) / 2.0;
			best = (left + right) / 2.0;

		}

		error = fabs(delta / best);

		if(error <= tol){

			// printf("The function zero = %f\n", best);

		}
		else{

			fBest = func(best);
			dBest = deriv(best);

			if(fLeft*fBest <= 0){

				right = best;
				fRight = fBest;

			}
			else{

				left = best;
				fLeft = fBest;

			}

		}

	}

	return(best);

}

int main(){

	double rootOne, rootTwo, rootThree, rootFour;
	rootOne = hybrid(0.1,0.3);
	printf("1st Root = %f\n", rootOne);
	rootTwo = hybrid(0.4,0.6);
	printf("2nd Root = %f\n", rootTwo);
	rootThree = hybrid(0.7,0.9);
	printf("3rd Root = %f\n", rootThree);
	rootFour = hybrid(0.95,1.1);
	printf("4th Root = %f\n", rootFour);

}
