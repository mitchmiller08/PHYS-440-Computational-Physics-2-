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
	double c = sqrt(0.2641);
	double v = 10.;
	double a = 10.;
	double b = 1.;
	funcValue = sqrt(x-v)*sin(a*c*sqrt(x))*cos(b*c*sqrt(x-v))+sqrt(x)*cos(a*c*sqrt(x))*sin(b*c*sqrt(x-v));
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

	double rootOne, rootTwo;
	rootOne = hybrid(11.5,12.5);
	printf("1st Root = %f\n", rootOne);
	rootTwo = hybrid(15.5,16.5);
	printf("2nd Root = %f\n", rootTwo);

}
