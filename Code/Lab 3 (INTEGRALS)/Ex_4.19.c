#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 3
// Exc. 4.19 2/28/12
// Gaussian Quadrature Method

double func(double x, double y){

	double funcValue;
	funcValue = pow(e,-x*y);
	return(funcValue);

}

double mapX(double x){

	double mapValue;
	double a = 0;
	double b = 2;
	mapValue = (b-a)/2.*x + (a+b)/2.;
	return(mapValue);

}

double mapY(double y){

	double mapValue;
	double a = 0;
	double b = 1;
	mapValue = (b-a)/2.*y + (a+b)/2.;
	return(mapValue);

}


double FofY(double y){

	double integral=0;
	int counter=0;
	int order=5;
	double x[5] = {0.9061798459386640,0.5384693101056831,0.0	       ,-0.5384693101056831,-0.9061798459386640};
	double w[5] = {0.2369268850561891,0.4786286704993665,0.5688888888888889, 0.4786286704993665, 0.2369268850561891};

	for(counter=0;counter<order;counter++){

		integral = integral + func(mapX(x[counter]),mapY(y))*w[counter];

	}

	integral = 0.5*integral;
	return(integral);

}

double quad5(){

	double integral=0;
	int counter=0;
	int order=5;
	double x[5] = {0.9061798459386640,0.5384693101056831,0.0	       ,-0.5384693101056831,-0.9061798459386640};
	double w[5] = {0.2369268850561891,0.4786286704993665,0.5688888888888889, 0.4786286704993665, 0.2369268850561891};

	for(counter=0;counter<order;counter++){

		integral = integral + FofY(x[counter])*w[counter];

	}

	return(integral);

}


int main(){

	double integral2,integral3,integral4,integral5;

	integral5 = quad5();

	printf("N = 5\t%.15lf\n",integral5);

}
