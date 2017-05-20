#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

// Mitchell Miller
// PHYS 440 Lab 1
// Exc. 2.1 1/18/12
// Newton - Raphson Method

double func(double x){

	// Specify desired function
	return(cos(x) - x);

}

double deriv(double x){

	// Provide the derivative of your function
	return(-sin(x) - 1.0);

}

int main() {

	// Initialize variables
	double xValue = 0.8;
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
