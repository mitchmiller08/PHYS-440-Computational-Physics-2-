/*
 *  Ex_7_1.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 4/3/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
	double numTimeSteps = 20;
	double currentTime = 0;
	double length = 1;
	double massDensity = 1e-3;
	double tension = 10;
	double xStep = 1e-2;
	double timeStep = 1e-6;
	int numXSteps = (int)(length / xStep);
	double maxTime = timeStep*numTimeSteps;
	int counter = 0;
	double* initialDisplacement;
	double* currentDisplacement;
	double* oldDisplacement;
	double* newDisplacement;
	double* dUdt;
	double exponent = 0;
	FILE* stringOutout;
	double epsilon;
	double xValue = 0;
	int i = 0;
	double waveVelocity = sqrt(tension / massDensity);
	stringOutout = fopen("waveStringPosition.txt", "w");
	initialDisplacement = malloc(numXSteps*sizeof(double));
	currentDisplacement = malloc(numXSteps*sizeof(double));
	oldDisplacement = malloc(numXSteps*sizeof(double));
	newDisplacement = malloc(numXSteps*sizeof(double));
	dUdt = malloc(numXSteps*sizeof(double));
	epsilon = pow((timeStep * waveVelocity / xStep),2);
	/********* initial displacement / conditions *****************/
	initialDisplacement[0] = 0;
	initialDisplacement[numXSteps-1] = 0;
	oldDisplacement[0] = 0;
	oldDisplacement[numXSteps-1] = 0;
	currentDisplacement[0] = 0;
	currentDisplacement[numXSteps-1] = 0;
	dUdt[0] = 0;
	dUdt[numXSteps-1] = 0;
	counter = 1;
	while (counter < numXSteps-1)
	{
		xValue = xStep * (double)counter;
		exponent = -100*pow((xValue - .5),2);
		initialDisplacement[counter] = exp(exponent);
		oldDisplacement[counter] = initialDisplacement[counter];
		dUdt[counter] = 0;
		counter++;
	}
	counter = 1;
	while (counter < numXSteps-1)
	{
		currentDisplacement[counter] = .5*epsilon*(initialDisplacement[counter+1] + initialDisplacement[counter-1])
										+ (1 - epsilon) * initialDisplacement[counter] + timeStep*dUdt[counter];
		counter++;
		
	}
	// main time loop
	while (currentTime < maxTime ) 
	{
		currentTime += timeStep;
		newDisplacement[0] = 0;
		newDisplacement[numXSteps-1] = 0;
		for (i = 1; i < numXSteps - 1; i++)
		{
			newDisplacement[i] = epsilon* (currentDisplacement[i+1] + currentDisplacement[i-1])
								+ 2 *(1-epsilon)*currentDisplacement[i] - oldDisplacement[i];
			printf("%e\n",newDisplacement[i]);
		}
		i = 0;
		while (i<numXSteps)
		{
			oldDisplacement[i] = currentDisplacement[i];
			currentDisplacement[i] = newDisplacement[i];
			printf("%e\n",currentDisplacement[i]);
			i++;
		}
		i = 0;
		while (i< numXSteps)
		{
			xValue = xStep * (double)i;
			fprintf( stringOutout, "%lf\t%e\n", xValue, currentDisplacement[i]);
			i++;
		}
	}
}
