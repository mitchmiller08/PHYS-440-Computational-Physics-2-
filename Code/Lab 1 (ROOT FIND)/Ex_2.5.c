#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

// Mitchell Miller
// PHYS 440 Lab 1
// Exc. 2.1 1/19/12
// Newton - Raphson Method

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

int main() {

	// Initialize variables
	double xValue = 0.167;
	double funcX, derX;
	double delta;
	double error = 1;
	double tol = 0.00000005;
	int counter = 0;

	while(error > tol){

		funcX = func(xValue);
		derX = deriv(xValue);
		delta = (-1) * funcX / derX;
		xValue = xValue + delta;
		error = fabs(delta / xValue);
		counter++;
		
		printf("The function zero = %f\n", xValue);
		printf("Number of loops = %i\n",counter);
		printf("Error = %e\n", error);

	}

}
