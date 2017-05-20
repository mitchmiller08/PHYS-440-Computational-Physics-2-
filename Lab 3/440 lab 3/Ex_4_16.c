/*
 *  Ex_4_16.c
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
double Euler = 2.718281828459045;
double functionInput(double xValue, double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	// perform mapping of function into range of {-1,1}
	// if upper and lower limits are already 1,-1, mapping will cause no change in xValue
	double functionVal = 0;
	xValue = ((UpperIntegrationLimit-lowerIntegrationLimit) / 2) * xValue 
			+ (UpperIntegrationLimit+lowerIntegrationLimit) / 2;
	functionVal = pow(Euler, -pow(xValue,2));
	return functionVal;
}
double gaussN2(double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	int nGauss = 2;
	double xGauss[2];
	double weightGauss[2];
	int counter=0;
	double integralVal = 0;
	xGauss[0] = -sqrt(pow(3,-1));
	xGauss[1] = sqrt(pow(3,-1));
	weightGauss[0] = 1.;
	weightGauss[1] = 1.;
	for (counter = 0; counter < nGauss; counter++)
	{
		integralVal += weightGauss[counter] * functionInput(xGauss[counter], UpperIntegrationLimit, lowerIntegrationLimit);
	}
	// perform last piece of mapping onto {-1,1}. if interval is already correct, no change is calculated.
	integralVal *= (UpperIntegrationLimit - lowerIntegrationLimit) / 2;
	return integralVal;
}
double gaussN3(double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	int nGauss = 3;
	double xGauss[3];
	double weightGauss[3];
	int counter=0;
	double integralVal = 0;
	xGauss[0] = -(sqrt(15) / 5);
	xGauss[1] = 0;
	xGauss[2] = sqrt(15) / 5;
	weightGauss[0] = 5/9.;
	weightGauss[1] = 8/9.;
	weightGauss[2] = 5/9.;
	for (counter = 0; counter < nGauss; counter++)
	{
		integralVal += weightGauss[counter] * functionInput(xGauss[counter], UpperIntegrationLimit, lowerIntegrationLimit);
	}
	// perform last piece of mapping onto {-1,1}. if interval is already correct, no change is calculated.
	integralVal *= (UpperIntegrationLimit - lowerIntegrationLimit) / 2;
	return integralVal;
}
double gaussN4(double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	int nGauss = 4;
	double xGauss[4];
	double weightGauss[4];
	int counter=0;
	double integralVal = 0;
	xGauss[0] = -.8611363115940526;
	xGauss[1] = -.3399810435848563;
	xGauss[2] = .3399810435848563;
	xGauss[3] = .8611363115940526;
	weightGauss[0] = .3478548451374539;
	weightGauss[1] = .6521451548625461;
	weightGauss[2] = .6521451548625461;
	weightGauss[3] = .3478548451374539;
	for (counter = 0; counter < nGauss; counter++)
	{
		integralVal += weightGauss[counter] * functionInput(xGauss[counter], UpperIntegrationLimit, lowerIntegrationLimit);
	}
	// perform last piece of mapping onto {-1,1}. if interval is already correct, no change is calculated.
	integralVal *= (UpperIntegrationLimit - lowerIntegrationLimit) / 2;
	return integralVal;
}
double gaussN5(double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	int nGauss = 5;
	double xGauss[5];
	double weightGauss[5];
	int counter=0;
	double integralVal = 0;
	xGauss[0] = - .9061798459386640;
	xGauss[1] = - .5384693101056831;
	xGauss[2] = 0;
	xGauss[3] = .5384693101056831;
	xGauss[4] = .9061798459386640;
	weightGauss[0] = .2369268850561891;
	weightGauss[1] = .4786286704993665;
	weightGauss[2] = .5688888888888889;
	weightGauss[3] = .4786286704993665;
	weightGauss[4] = .2369268850561891;
	for (counter = 0; counter < nGauss; counter++)
	{
		integralVal += weightGauss[counter] * functionInput(xGauss[counter], UpperIntegrationLimit, lowerIntegrationLimit);
	}
	// perform last piece of mapping onto {-1,1}. if interval is already correct, no change is calculated.
	integralVal *= (UpperIntegrationLimit - lowerIntegrationLimit) / 2;
	return integralVal;
}
int main()
{
	double UpperIntegrationLimit = 1;
	double lowerIntegrationLimit = 0;
	double n2Gauss = 0;
	double n3Gauss = 0;
	double n4Gauss = 0;
	double n5Gauss = 0;
	double Exact = 0;
	Exact = 0.746824132812427;
	n2Gauss = gaussN2(UpperIntegrationLimit,lowerIntegrationLimit);
	n3Gauss = gaussN3(UpperIntegrationLimit,lowerIntegrationLimit);
	n4Gauss = gaussN4(UpperIntegrationLimit,lowerIntegrationLimit);
	n5Gauss = gaussN5(UpperIntegrationLimit,lowerIntegrationLimit);
	printf("N=2\t						N=3\t						N=4\t						N=5\t						Exact\n");
	printf("%.8f	\t%.8f		\t%.8f		\t%.8f		\t%.8f\n",n2Gauss,n3Gauss,n4Gauss,n5Gauss,Exact);
	return 0;
}

