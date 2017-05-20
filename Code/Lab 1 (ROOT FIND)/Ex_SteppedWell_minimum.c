#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

// Mitchell Miller
// PHYS 440 Lab 1
// Exc. 2.1 1/19/12
// Hybrid Method
// Calculate smallest width of 'a' to have a bound state

double func(double x, double a){

	// Specify desired function
	double funcValue;
	double c = sqrt(0.2641);
	double v = 10.;
	funcValue = (v-2.0*x)*sin(a*c*sqrt(x))+2.0*x*sqrt(v/x-1.0)*cos(a*c*x);
	return(funcValue);

}

double deriv(double dLeft, double dRight, double w){

	// Provide the approximate derivative of your function
	double derValue;
	derValue = (func(dRight,w)-func(dLeft,w))/(dRight-dLeft);
	return(derValue);

}

int main(){

	double width = 1.365;
	int index=0;

	while(index < 10){	

		// Initialize variables
		double left = 5.;
		double right = 10.;
		double best=0., fBest=0., dBest=0.;
		double funcX=0., derX=0.;
		double delta=0.;
		double error = 1;
		double tol = 0.00000005;
		double fLeft = func(left, width);
		double fRight = func(right, width);

		if(fabs(fLeft <= fabs(fRight))){

			best = left;
			fBest = fLeft;

		}
		else{

			best = right;
			fBest = fRight;

		}

		dBest = deriv(left,right,width);

		while(error > tol){

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

				fBest = func(best,width);
				dBest = deriv(left,right,width);

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

		printf("Root = %f\n", best);
		printf("Well Size = %f\n", width);
		width = width - 0.001;
		index++;

	}

}
