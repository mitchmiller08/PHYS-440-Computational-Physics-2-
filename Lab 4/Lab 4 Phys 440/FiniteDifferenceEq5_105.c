/*
 *  FiniteDifferenceEq5_105.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 4/1/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/************************************************/
/*												*/
/*		y"(x) -5y'(x) + 10y(x) = 10x			*/
/*												*/
/*		y"(i) = y(i+1) - 2y(i) + y(i-1)			*/
/*			-------------------------			*/
/*					h^2							*/
/*												*/
/*		y'(i) = y(i+1) - y(i-1)					*/
/*				---------------					*/
/*					2*h							*/
/************************************************/

double function(double* yOldArray, double* ynewArray,double xValue, double step,int arrayCounter,double constant1, double constant2, double constant3)
{	
	
	ynewArray[arrayCounter] = (1 / (2-constant3)) * (constant1*yOldArray[arrayCounter+1] + constant2*yOldArray[arrayCounter-1] - constant3*xValue);
	return ynewArray[arrayCounter];
}
double* initializeOldArray(double* yOldArray, int numArrayValues)
{
	int iteration = 0;
	double value = 0;
	while (iteration < numArrayValues-1)
	{
		value = 100*(double)(iteration/10.); // Initial guesses for Y_i values
		yOldArray[iteration] = value;
		iteration++;
	}
	return yOldArray;
}
double* initializeNewArray(double* yNewArray, int numArrayValues)
{
	int iteration = 0;
	while (iteration<numArrayValues)
	{
		yNewArray[iteration] = 0;
		iteration++;
	}
	return yNewArray;
}
int main()
{
	FILE * output;
	output = fopen("Equation_5_105.txt","w");
	double* yNewArray;
	double* yOldArray;
	double relError = 0;
	double upperBoundary = 1;
	double lowerBoundary = 0;
	double step = .1;
	double xValue = 0;
	double tolerance = 2e-4;
	double constant1 = 0;
	double constant2 = 0;
	double constant3 = 0;
	int DONE = 0;
	int iteration = 0;
	int arrayCounter = 0;
	int numArrayValues = 1 + (int)(upperBoundary - lowerBoundary) / step;
	yNewArray = malloc(numArrayValues*sizeof(double));
	yOldArray = malloc(numArrayValues*sizeof(double));
	yNewArray = initializeNewArray(yNewArray, numArrayValues);
	yOldArray = initializeOldArray(yOldArray, numArrayValues);
	// Boundary Conditions
	yOldArray[0] = 0;
	yOldArray[numArrayValues-1] = 100;
	yNewArray[0] = 0;
	yNewArray[numArrayValues-1] = 100;
	
	constant1 = 1 - 2.5*step;
	constant2 = 1 + 2.5*step;
	constant3 = 10*pow(step,2);
	while (DONE == 0)
	{
		
		
			DONE = 1;
			arrayCounter = 1;
		printf("%i\n\n", iteration);
			fprintf(output,"%lf\t%lf\n", 0,yOldArray[0]);
			while (arrayCounter< numArrayValues-1)
			{
				xValue = .1 + .1*(arrayCounter-1);
				fprintf(output,"%lf\t%lf\n", xValue, yOldArray[arrayCounter]);
				yNewArray[arrayCounter] = function(yOldArray,yNewArray, xValue, step, arrayCounter, constant1, constant2, constant3);
				relError = fabs((yNewArray[arrayCounter] - yOldArray[arrayCounter]) / yNewArray[arrayCounter]);
				if (relError > tolerance)
				{
					DONE = 0;
				}
				yOldArray[arrayCounter] = yNewArray[arrayCounter];
				arrayCounter++;
			}
			fprintf(output,"%lf\t%lf\n",upperBoundary,yOldArray[numArrayValues-1]);
			iteration++;
		
	}
	arrayCounter = 0;
	while (arrayCounter<numArrayValues)
	{
		xValue = arrayCounter*step;
		fprintf(output, "%lf\t%lf",xValue, yOldArray[arrayCounter]);
		arrayCounter++;
		
	}
	
	
	
	
	return 47;
	
	
}
