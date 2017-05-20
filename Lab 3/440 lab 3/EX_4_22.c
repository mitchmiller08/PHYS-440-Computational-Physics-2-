/*
 *  EX_4_22.c
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
/****************************/
/*		2D Integration		*/
/****************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
double functionInput(double xValue, double XP, double YP, double yValue, double UpperXIntegrationLimit, double lowerXIntegrationLimit, double UpperYIntegrationLimit, double lowerYIntegrationLimit)
{
	// perform mapping of function into range of {-1,1}
	// if upper and lower limits are already 1,-1, mapping will cause no change in xValue
	double functionVal = 0;
	xValue = ((UpperXIntegrationLimit-lowerXIntegrationLimit) / 2) * xValue 
			+ (UpperXIntegrationLimit+lowerXIntegrationLimit) / 2;
	functionVal = 1 / sqrt(pow(xValue-XP,2) + pow(yValue - YP,2));
	return functionVal;
}
double gaussN5(double UpperXIntegrationLimit, double lowerXIntegrationLimit, double UpperYIntegrationLimit, double lowerYIntegrationLimit)
{
	int nGauss = 5;
	double Gauss[5];
	double weightGauss[5];
	double xP = 0;
	double yP = 0;
	FILE *file;
	int Xcounter=0;
	int Ycounter = 0;
	double integralVal = 0;
	file = fopen("square.txt","w");
	Gauss[0] = - .9061798459386640;
	Gauss[1] = - .5384693101056831;
	Gauss[2] = 0;
	Gauss[3] = .5384693101056831;
	Gauss[4] = .9061798459386640;
	weightGauss[0] = .2369268850561891;
	weightGauss[1] = .4786286704993665;
	weightGauss[2] = .5688888888888889;
	weightGauss[3] = .4786286704993665;
	weightGauss[4] = .2369268850561891;
	for (xP = -20; xP< 20; xP+=.1)
	{
		for (yP = -20; yP < 20; yP+=.1)
		{
			integralVal = 0;
			for (Xcounter = 0; Xcounter < nGauss; Xcounter++)
			{
				for (Ycounter = 0; Ycounter < nGauss; Ycounter++)
				{
					integralVal += weightGauss[Xcounter] * weightGauss[Ycounter] * functionInput(Gauss[Xcounter], xP, yP, Gauss[Ycounter], UpperXIntegrationLimit, lowerXIntegrationLimit, UpperYIntegrationLimit, lowerYIntegrationLimit);
				}
			}
			integralVal *= (UpperXIntegrationLimit - lowerYIntegrationLimit) / 2;
			integralVal *= (UpperYIntegrationLimit - lowerYIntegrationLimit) / 2;
			fprintf(file, "%lf\t%lf\t%lf\n", xP,yP, integralVal);
		}
	}
	// perform last piece of mapping onto {-1,1}. if interval is already correct, no change is calculated.
	integralVal *= (UpperXIntegrationLimit - lowerYIntegrationLimit) / 2;
	integralVal *= (UpperYIntegrationLimit - lowerYIntegrationLimit) / 2;
	return integralVal;
}
int main()
{
	double xValue;
	double yValue;
	double UpperXIntegrationLimit = 1;
	double UpperYIntegrationLimit = 1;
	double lowerXIntegrationLimit = -1;
	double lowerYIntegrationLimit = -1;
	double integralVal = 0;
	integralVal = gaussN5(UpperXIntegrationLimit, lowerXIntegrationLimit, UpperYIntegrationLimit, lowerYIntegrationLimit);

}
