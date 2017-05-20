#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define e 2.71828182845

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.13 2/17/12

double function(double x){
	//Specify your desired function
	return(x*pow(e,x));
	//return(cos(x));
}

double forwardDiff(double xValue, double stepSize){
	double differentiatedValue;

	differentiatedValue = (function(xValue+stepSize) - function(xValue)) / stepSize;
	return(differentiatedValue);

}

double centralDiff(double xValue, double stepSize){
	double differentiatedValue;

	differentiatedValue = (function(xValue+stepSize)-function(xValue-stepSize)) / (2.*stepSize);
	return(differentiatedValue);

}

double extrapolatedDiff(double xValue, double stepSize){
	double differentiatedValue;

	differentiatedValue = (function(xValue) - function(xValue-stepSize)) / stepSize;
	return(differentiatedValue);

}

double forwardSecDiff(double xValue, double stepSize){

	double differentiatedValue;

	differentiatedValue = ((function(xValue) - 2.0*function(xValue+stepSize) + function(xValue+2.0*stepSize)) / (stepSize * stepSize));

	return(differentiatedValue);

}

double centralSecDiff(double xValue, double stepSize){

	double differentiatedValue;

	differentiatedValue = ((function(xValue+stepSize) - 2.0*function(xValue) + function(xValue-stepSize)) / (stepSize * stepSize));

	return(differentiatedValue);

}

int main(){
	//Initialize variables
	double stepSize=0.5;
	int counter;
	double forSolution,cenSolution,extSolution,for2Solution,cen2Solution;
	double forError,cenError,extError,for2Error,cen2Error;
	double firstDeriv = 22.1671682967920;
	double secondDeriv = 29.5562243957226;
	FILE *outfile;
	outfile=fopen("Ex_3.13.error.dat","w");
	fprintf(outfile,"#N\teps_F\teps_C\teps_E\teps_F2\teps_C2\n");

	//Calculate values
	for(counter=0;counter<45;counter++){
		forSolution = forwardDiff(2.0,stepSize);
		cenSolution = centralDiff(2.0,stepSize);
		extSolution = extrapolatedDiff(2.0,stepSize);
		for2Solution = forwardSecDiff(2.0,stepSize);
		cen2Solution = centralSecDiff(2.0,stepSize);

		//Calculate error
		forError = fabs((firstDeriv-forSolution) / firstDeriv);
		cenError = fabs((firstDeriv-cenSolution) / firstDeriv);
		extError = fabs((firstDeriv-extSolution) / firstDeriv);
		for2Error = fabs((secondDeriv-for2Solution) / secondDeriv);
		cen2Error = fabs((secondDeriv-cen2Solution) / secondDeriv);
		
		//Print results
		printf("h = %lf\n", stepSize);
		printf("Forward Error = %e\n", forError);
		printf("Central Error = %e\n", cenError);
		printf("Extrapolated Error = %e\n", extError);
		printf("Forward 2nd Deriv Error = %e\n", for2Error);
		printf("Central 2nd Deriv Error = %e\n", cen2Error);
		fprintf(outfile,"%lf\t%.20lf\t%.20lf\t%.20lf\t%.20lf\t%.20lf\n",stepSize,forError,cenError,extError,for2Error,cen2Error);
		stepSize = stepSize - 0.01;
	}

	fclose(outfile);

}
