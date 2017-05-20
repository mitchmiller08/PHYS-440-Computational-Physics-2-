#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 3
// Fig 4.12 2/28/12
// Monte Carlo Method

double func(double x){

	double funcValue;
	funcValue = pow(e,x);
	return(funcValue);

}

double monteCarlo(){

	double sum=0.;
	int counter1,counter2;
	double N;
	double monte,x,error;

	for(counter1=1;counter1<401;counter1++){

		for(counter2=1;counter2<11;counter2++){

			x = drand48();
			sum = sum + func(x);

		}

		N = counter1 * 10.;

		monte = sum/N;

		error = fabs(monte - (e-1.))/(e-1.);

		//printf("N = %lf\t%lf\t%lf\n",N,monte,error);

	}

	return(monte);

}	

int main(){

	FILE *outfile;
	outfile = fopen("Fig_4.12.400.dat","w");
	double integral;
	int counter;
	
	for(counter=0;counter<10000;counter++){
	
		integral = monteCarlo();
		fprintf(outfile,"%lf\n",integral);
		
	}

}	
