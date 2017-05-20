/*
 *  Ex_7_8.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 4/4/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double **initiateHeatFunction (int numXsteps, int numYSteps, double tempLeft, double tempRight, double tempTop, double tempBottom, double xStep, double yStep)
{
	double** heatFunction;
	double xValue = 0;
	double yValue = 0;
	int xCounter, yCounter;
	xCounter = 0;
	yCounter = 0;
	heatFunction = (double**)malloc(numXsteps*sizeof(double*));
	while (xCounter < numXsteps)
	{
		heatFunction[xCounter] = (double*)malloc(numYSteps*sizeof(double));
		xCounter++;
	}
	xCounter = 0;
	while (xCounter < numXsteps)
	{
		yCounter = 0;
		while (yCounter < numYSteps)
		{
			heatFunction[xCounter][yCounter] = .1;
			yCounter++;
		}
		
		xCounter++;
	}
	xCounter = yCounter = 0;
	while (yCounter < numYSteps)
	{
		heatFunction[0][yCounter] = tempLeft;
		heatFunction[numXsteps-1][yCounter] = tempRight;
		yCounter++;
	}
	while (xCounter < numXsteps) 
	{
		heatFunction[xCounter][0] = tempBottom;
		heatFunction[xCounter][numYSteps-1] = tempTop;
		xCounter++;
	}
	xCounter = yCounter = 0;
	while (xCounter < numXsteps)
	{
		xValue = xStep*(double)xCounter;
		yCounter = 0;
		while (yCounter < numYSteps)
		{
			yValue = yStep*(double)yCounter;
			printf("%lf\t%lf\t%lf\n",xValue,yValue,heatFunction[xCounter][yCounter]);
			yCounter++;
		}
		xCounter++;
	}
	return heatFunction;
}
double successiveOverRelaxation ( double **heatFunction, int xCounter, int yCounter, double xStep, double yStep, int numXsteps, int numYSteps)
{
	double newHeatValue = 0;
	double xForward = heatFunction[xCounter+1][yCounter];
	double xBack = heatFunction[xCounter-1][yCounter];
	double yForward = heatFunction[xCounter][yCounter+1];
	double yBack = heatFunction[xCounter][yCounter-1];
	newHeatValue = ( (xForward +xBack) / (pow(xStep, 2))
					+ ( yForward + yBack) / (pow(yStep , 2))) * pow(xStep*yStep, 2);
					
	newHeatValue /= (2*( pow(xStep, 2) + pow(yStep, 2)));
	return newHeatValue;
}
int main()
{
	double xStep, yStep, tempLeft, tempTop, tempRight, tempBottom, alpha, newHeatValue, xLength, yLength, Pi, tolerance, xValue, yValue;
	double **heatFunction;
	int numXsteps, numYSteps, iteration, xCounter, yCounter, DONE, maxIterations;
	FILE* heatOutPut;
	heatOutPut = fopen("heatData.txt", "w");
	xCounter = 0;
	yCounter = 0;
	heatFunction = (double**)malloc(numXsteps*sizeof(double*));
	while (xCounter < numXsteps)
	{
		heatFunction[xCounter] = (double*)malloc(numYSteps*sizeof(double));
		xCounter++;
	}
	DONE = 0; // 0 ~~> FALSE
	Pi = 3.14159265358979323;
	maxIterations = 100;
	xLength = 1; //  meter
	yLength = 1; //  meter
	numXsteps = 33;
	numYSteps = numXsteps;
	xStep = (xLength / (double)numXsteps);
	yStep = (yLength / (double)numYSteps); 
	tempTop = 0; // celsius
	tempRight = tempTop;
	tempLeft = 100.; // celsius
	tempBottom = tempLeft;
	tolerance = 5e-5;
	alpha = (4 / (2 + sqrt(4 - pow(( cos(Pi / (double)numXsteps) + cos(Pi / (double)numYSteps)), 2)))) - 1;
	heatFunction = initiateHeatFunction( numXsteps, numYSteps, tempLeft, tempRight, tempTop, tempBottom, xStep, yStep);
	while (DONE == 0)
	{
		DONE = 1;
		iteration++;
		if (iteration > maxIterations)
		{
			break;
		}
		xCounter = yCounter = 1;
		while (xCounter < numXsteps - 1)
		{
			yCounter = 1;
			while (yCounter < numYSteps - 1)
			{
				newHeatValue = successiveOverRelaxation( heatFunction, xCounter, yCounter, xStep, yStep, numXsteps, numYSteps);
				if( fabs( (newHeatValue - heatFunction[xCounter][yCounter]) / newHeatValue) > tolerance)
				{
					DONE = 0;
				}
				heatFunction[xCounter][yCounter] = newHeatValue+alpha*(newHeatValue - heatFunction[xCounter][yCounter]);
				yCounter++;
			}
			xCounter++;
		}
	}
	xCounter = 0;
	yCounter = 0;
	while (xCounter < numXsteps)
	{
		xValue = xStep*(double)xCounter;
		yCounter = 0;
		while (yCounter < numYSteps)
		{
			yValue = yStep*(double)yCounter;
			fprintf(heatOutPut, "%lf\t%lf\t%lf\n",xValue, yValue, heatFunction[xCounter][yCounter]);
			yCounter++;
		}
		xCounter++;
	}
	return 0;
}

