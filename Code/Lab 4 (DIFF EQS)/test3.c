/*
 *  EX_5_39.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 4/2/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
int main()
{
	double length = 1;
	double massDensity = .954;
	double tension = 1000;
	double Pi = 3.1415926535897932;
	double* lhAMatrix;
	double* bMatrix;
	double* cMatrix;

	double fundamental = (Pi/length)*sqrt(tension/massDensity);
	double muO = .954; // g / m
	double delta = .5; // g / m^2
	double step = 1/1000; 
	int numPoints = 1000;
	double alpha = .5*pow((Pi/numPoints),2);
	int n1 = numPoints - 2; // number of interior points
	double determ1 = 0;
	double determ2 = 0;
	double determAllNew = 0;
	double determAllOld = 0;
	int counter = 0;
	int determCounter = 0;
	double eigenValue = 0;
	double lambda = 0;
	double xValue = 0;
	lhAMatrix = malloc(numPoints*sizeof(double));
	bMatrix = malloc(numPoints*sizeof(double));
	cMatrix = malloc(numPoints*sizeof(double));
	step = length / numPoints;
	
	while (counter<numPoints)
	{
		xValue = step*(double)counter;
		
		massDensity = muO + (xValue - .5 * length) * delta;
		//massDensity = muO + pow((xValue - .5 * length),2) * pow(delta,2);
		//massDensity = muO;
		lhAMatrix[counter] = tension / (pow(step,2)*massDensity); 
		//scaling to avoid large numbers;
		lhAMatrix[counter] = lhAMatrix[counter] / (2*tension/(muO*pow(step,2)));
		bMatrix[counter] = -2*lhAMatrix[counter];
		cMatrix[counter] = lhAMatrix[counter];
		counter++;
	}
	counter = 1;
	
	while (counter <  n1+1)
	{
		lambda = 50*alpha*counter/n1;
		determAllNew = bMatrix[0] + lambda;
		
		if (numPoints != 1)
		{
			determ1 = determAllNew;
			determAllNew = (bMatrix[1] + lambda)*determ1 - lhAMatrix[1]*cMatrix[0];
		}
		if (numPoints != 2)
		{
			determ2 = determAllNew;
			determCounter = 3;
			while (determCounter < n1)
			{
				determAllNew = (bMatrix[determCounter] + lambda)*determ2 - lhAMatrix[determCounter]*cMatrix[determCounter-1]*determ1;
				determ1 = determ2;
				determ2 = determAllNew;
				determCounter++;
			}
		}
		determAllNew *= 1e300; //this is done to ensure that the multiplication check for zero crossing in the next line does not produce an underflow.
		if (determAllNew*determAllOld < 0)
		{
			eigenValue = sqrt(2*tension*lambda / (muO*pow(step,2)));
			eigenValue /= fundamental;
			printf("%lf\t%lf\n",eigenValue,determAllNew);
		}
		determAllOld = determAllNew;
		counter++;
	}
	
}


