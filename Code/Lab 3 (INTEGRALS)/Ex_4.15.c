#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

// Mitchell Miller
// PHYS 440 Lab 3
// Exc. 4.15 2/28/12
// Gaussian Quadrature Method

double func(double x){

	double funcValue;
	funcValue = x*x*x*x*x*x*x;
	return(funcValue);

}

double map(double x){

	double mapValue;
	double a = 0;
	double b = 1;
	mapValue = (b-a)/2.*x + (a+b)/2.;
	return(mapValue);

}

double quad2(){

	double integral=0;
	int counter=0;
	int order=2;
	double x[2] = {0.5773502691896258,-0.5773502691896258};
	double w[2] = {1.0		 ,1.0		     };

	for(counter=0;counter<order;counter++){

		integral = integral + func(map(x[counter]))*w[counter];

	}
	
	integral = 0.5*integral;
	return(integral);

}

double quad3(){

	double integral=0;
	int counter=0;
	int order=3;
	double x[3] = {0.7745966692414834,0.0		    ,-0.7745966692414834};
	double w[3] = {0.5555555555555556,0.8888888888888889, 0.5555555555555556};

	for(counter=0;counter<order;counter++){

		integral = integral + func(map(x[counter]))*w[counter];

	}

	integral = 0.5*integral;
	return(integral);

}

double quad4(){

	double integral=0;
	int counter=0;
	int order=4;
	double x[4] = {0.8611363115940526,0.3399810435848563,-0.3399810435848563,-0.8611363115940526};
	double w[4] = {0.3478548451374539,0.6521451548625461, 0.6521451548625461, 0.3478548451374539};

	for(counter=0;counter<order;counter++){

		integral = integral + func(map(x[counter]))*w[counter];

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

		integral = integral + func(map(x[counter]))*w[counter];

	}

	integral = 0.5*integral;
	return(integral);

}

int main(){

	double integral2,integral3,integral4,integral5;

	integral2 = quad2();
	integral3 = quad3();
	integral4 = quad4();
	integral5 = quad5();

	printf("N = 2\t%.15lf\n",integral2);
	printf("N = 3\t%.15lf\n",integral3);
	printf("N = 4\t%.15lf\n",integral4);
	printf("N = 5\t%.15lf\n",integral5);

}
