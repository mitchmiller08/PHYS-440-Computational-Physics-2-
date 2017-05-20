/*
 *  EX_4_18.c
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
double euler = 2.718281828459045;
double functionInput(double xValue, double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	double functionVal = 0;
	functionVal = pow(xValue,3) / (pow(euler,xValue) - 1);
	return functionVal;
}
double gaussLaguerreN2(double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	int nGauss = 2;
	double xGaussLaguerre[2];
	double weightGaussLaguerre[2];
	int counter=0;
	double integralVal = 0;
	xGaussLaguerre[0] = .58578643762690495;
	xGaussLaguerre[1] = 3.4142135623730950;
	weightGaussLaguerre[0] = .85355339059327376;
	weightGaussLaguerre[1] = .14644660940672624;
	for (counter = 0; counter < nGauss; counter++)
	{
		integralVal += weightGaussLaguerre[counter] * pow(euler,xGaussLaguerre[counter])* functionInput(xGaussLaguerre[counter], UpperIntegrationLimit, lowerIntegrationLimit);
	}
	return integralVal;
}
double gaussLaguerreN4(double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	int nGauss = 4;
	double xGaussLaguerre[4];
	double weightGaussLaguerre[4];
	int counter=0;
	double integralVal = 0;
	xGaussLaguerre[0] = .32254768961939231;
	xGaussLaguerre[1] = 1.7457611011583466;
	xGaussLaguerre[2] = 4.5366202969211280;
	xGaussLaguerre[3] = 9.3950709123011331;
	weightGaussLaguerre[0] = .60315410434163360;
	weightGaussLaguerre[1] = .35741869243779969;
	weightGaussLaguerre[2] = .038887908515005384;
	weightGaussLaguerre[3] = .00053929470556132745;
	for (counter = 0; counter < nGauss; counter++)
	{
		integralVal += weightGaussLaguerre[counter] *pow(euler,xGaussLaguerre[counter])* functionInput(xGaussLaguerre[counter], UpperIntegrationLimit, lowerIntegrationLimit);
	}
	return integralVal;
}
double gaussLaguerreN6(double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	int nGauss = 6;
	double xGaussLaguerre[6];
	double weightGaussLaguerre[6];
	int counter=0;
	double integralVal = 0;
	xGaussLaguerre[0] = .22284660417926069;
	xGaussLaguerre[1] = 1.1889321016726230;
	xGaussLaguerre[2] = 2.9927363260593141;
	xGaussLaguerre[3] = 5.7751435691045105;
	xGaussLaguerre[4] = 9.8374674183825899;
	xGaussLaguerre[5] = 15.982873980601702;
	weightGaussLaguerre[0] = .45896467394996359;
	weightGaussLaguerre[1] = .41700083077212099;
	weightGaussLaguerre[2] = .11337338207404498;
	weightGaussLaguerre[3] = 1.0399197453149075e-2;
	weightGaussLaguerre[4] = 2.6101720281493206e-4;
	weightGaussLaguerre[5] = 8.9854790642962124e-7;
	for (counter = 0; counter < nGauss; counter++)
	{
		integralVal += weightGaussLaguerre[counter] *pow(euler,xGaussLaguerre[counter])* functionInput(xGaussLaguerre[counter], UpperIntegrationLimit, lowerIntegrationLimit);
	}
	return integralVal;
}
double gaussLaguerreN8(double UpperIntegrationLimit, double lowerIntegrationLimit)
{
	int nGauss = 8;
	double xGaussLaguerre[8];
	double weightGaussLaguerre[8];
	int counter=0;
	double integralVal = 0;
	xGaussLaguerre[0] =  1.7027963230510100e-1;
	xGaussLaguerre[1] =  9.0370177679937991e-1;
	xGaussLaguerre[2] =  2.2510866298661307;
	xGaussLaguerre[3] =  4.2667001702876588;
	xGaussLaguerre[4] =  7.0459054023934657;
	xGaussLaguerre[5] =  1.0758516010180995e+1;
	xGaussLaguerre[6] =  1.5740678641278005e+1;
	xGaussLaguerre[7] =  2.2863131736889264e+1;
	weightGaussLaguerre[0] = 3.6918858934163753e-1;
	weightGaussLaguerre[1] = 4.1878678081434296e-1;
	weightGaussLaguerre[2] = 1.7579498663717181e-1;
	weightGaussLaguerre[3] = 3.3343492261215652e-2;
	weightGaussLaguerre[4] = 2.7945362352256725e-3;
	weightGaussLaguerre[5] = 9.0765087733582131e-5;
	weightGaussLaguerre[6] = 8.4857467162725315e-7;
	weightGaussLaguerre[7] = 1.0480011748715104e-9;
	for (counter = 0; counter < nGauss; counter++)
	{
		integralVal += weightGaussLaguerre[counter] *pow(euler,xGaussLaguerre[counter])* functionInput(xGaussLaguerre[counter], UpperIntegrationLimit, lowerIntegrationLimit);
	}
	return integralVal;
}
int main()
{
	double UpperIntegrationLimit = 1;
	double lowerIntegrationLimit = 0;
	double n2GaussLaguerre = 0;
	double n4GaussLaguerre= 0;
	double n6GaussLaguerre = 0;
	double n8GaussLaguerre = 0;
	double Exact = 0;
	Exact = 6.49393940226684;
	n2GaussLaguerre = gaussLaguerreN2(UpperIntegrationLimit,lowerIntegrationLimit);
	n4GaussLaguerre = gaussLaguerreN4(UpperIntegrationLimit,lowerIntegrationLimit);
	n6GaussLaguerre = gaussLaguerreN6(UpperIntegrationLimit,lowerIntegrationLimit);
	n8GaussLaguerre = gaussLaguerreN8(UpperIntegrationLimit,lowerIntegrationLimit);
	printf("N=2\t						N=4\t					N=6\t						N=8\t					Exact\n");
	printf("%.8f		\t%.8f		\t%.8f		\t%.8f		\t%.8f\n",n2GaussLaguerre,n4GaussLaguerre,n6GaussLaguerre,n8GaussLaguerre,Exact);
	return 0;
}



