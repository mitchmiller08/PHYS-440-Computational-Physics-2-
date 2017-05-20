/*
 *  EX_4_24.c
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "rngs.h"
//#include "rvgs.h"
#include "rvgs.c"
#include "rngs.c"



double euler = 2.718281828459045;
double functionInput(double xValue)
{
		double functionVal = 0;
		functionVal = pow(euler, sqrt(1+2*xValue)-1) / sqrt(1+2*xValue);
		return functionVal;
	
}
double functionInputSquared(double xValue)
{
	double functionVal = 0;
	functionVal = pow(euler, sqrt(1+2*xValue)-1) / sqrt(1+2*xValue);
	functionVal*=functionVal;
	return functionVal;
}
double monteCarloIntegration (double upperLimit, double lowerLimit)
{
	double randomXValue = 0;
	double functionVal = 0;
	randomXValue = Uniform(lowerLimit,upperLimit);
	functionVal = functionInput(randomXValue);
	return functionVal;
}
double monteCarloIntegrationSigma(double upperLimit, double lowerLimit)
{
	double randomXValue = 0;
	double functionVal = 0;
	randomXValue = Uniform(lowerLimit,upperLimit);
	functionVal = functionInputSquared(randomXValue);
	return functionVal;
	
}
int* MCHistogramInitialize(int numBins)
{
	int* MCHistogram;
	int counter=0;
	MCHistogram = (int*)malloc(numBins*sizeof(int));
	for (counter=0;counter<numBins;counter++)
	{
		MCHistogram[counter] = 0;
	}
	return MCHistogram;
	
}
int* histogramInput(double integralVal, double binWidth,int numBins, int* MCHistogram)
{
	int counter=0;
	for (counter = 0; counter < numBins; counter++)
	{
		if (integralVal <  (1.5 + binWidth*counter))
		{
			MCHistogram[counter] += 1;
			break;
		}
	}
	return MCHistogram;

}
int main()
{
	double upperLimit  =		(3 / 2.);
	double lowerLimit =			0;
	double integrationWidth =				0;
	double numPoints =			100;
	double functionAvg =		0;
	double binWidth =			.002;
	int*   MCHistogram;
	int    iterations =			10000;
	int    iterationCounter =	0;
	int    monteCarloCounter =	0;
	int	   numBins =			0;
	double functionAvgSquared = 0;
	double sigma = 0;
	
	FILE *file;
	file = fopen("MCIntegration.txt","w");
	
	integrationWidth = upperLimit - lowerLimit;
	// .5 is the range over which I want to display the histogram, { 1.5, 2 }
	numBins = (int)(.5 / binWidth); 
	// resize the histogram;
	MCHistogram = (int*)malloc(numBins*sizeof(int));
	// intialize the histogram to zero;
	MCHistogram = MCHistogramInitialize(numBins);
	
	for (iterationCounter = 0; iterationCounter < iterations; iterationCounter++)
	{
		functionAvg = 0;
		functionAvgSquared = 0;
		sigma = 0;
		//perform monte-carlo integration
		for (monteCarloCounter = 0; monteCarloCounter < numPoints; monteCarloCounter++)
		{
			functionAvg += monteCarloIntegration(upperLimit,lowerLimit);
			functionAvgSquared += monteCarloIntegrationSigma(upperLimit,lowerLimit);
		}
		
		sigma = sqrt((functionAvgSquared/numPoints) - pow((functionAvg/numPoints),2))/(numPoints-1);
		functionAvg *= integrationWidth / numPoints;
		MCHistogram = histogramInput(functionAvg,binWidth,numBins,MCHistogram);
		
		printf("%lf\n", sigma);
	}
	for (iterationCounter = 0; iterationCounter < numBins; iterationCounter++)
	{
		fprintf(file, "%lf\t%i\n",(1.5+binWidth*iterationCounter), MCHistogram[iterationCounter]);
	}
	fclose(file);
	return 0;
}

