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
	funcValue = sin(x) - x/2.;
	return(funcValue);

}

double deriv(double dLeft, double dRight){

	// Provide the approximate derivative of your function
	double derValue;
	derValue = (func(dRight)-func(dLeft))/(dRight-dLeft);
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

	dBest = deriv(left,right);

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
			dBest = deriv(left,right);

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

	double rootOne;
	rootOne = hybrid(pi / 2., pi);
	printf("1st Root = %f\n", rootOne);

}
