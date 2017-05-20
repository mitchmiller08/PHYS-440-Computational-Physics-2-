#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

// Mitchell Miller
// PHYS 440 Lab 1
// Exc. 2.1 1/18/12
// Biseciton Method

double function(double x){

	// Specify desired function
	return(cos(x)-x);

}

int main() {

	// Initialize variables
	double left = 0;
	double right = pi/2;
	double funcLeft, funcRight;
	double mid, funcMid;
	double error;
	double tol = 0.00000005;
	int counter = 0;
	double product;

	funcRight = function(right);
	funcLeft = function(left);
	error = 10*tol;

	while(error > tol){

		mid = (right + left) / 2;
		funcMid = function(mid);
		product = funcLeft * funcMid;
		counter++;
		
		if(product <= 0.0){	// Restricts zero to left half

			right = mid;
			funcRight = funcMid;

		}
		else{			// Restricts zero to right half

			left = mid;
			funcLeft = funcMid;

		}

		error = fabs((right - left)/mid);

		// Print zero and number of loops
		printf("The function zero = %f\n", mid);
		printf("Number of loops = %i\n",counter);
		printf("Error = %e\n", error);

	}

}
