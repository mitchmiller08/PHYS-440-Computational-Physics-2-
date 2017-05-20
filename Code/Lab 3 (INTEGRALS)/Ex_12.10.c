#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 3
// Ex 12.10 2/28/12
// Monte Carlo Method

double sumRand(){

	double sum=0.;
	int counter;

	for(counter=0;counter<12;counter++){

		sum = sum + drand48();

	}

	return(sum);

}

int main(){

	int counter;
	double sum;

	FILE *outfile;
	outfile = fopen("Ex_12.10.dat","w");

	for(counter=0;counter<10000;counter++){

		sum = sumRand();
		fprintf(outfile,"%lf\n",sum);

	}

}
