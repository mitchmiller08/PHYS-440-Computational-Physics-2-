/*
 *  RombergIntegration.h
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
 
/********************************************************/
/*				Romberg Integration w/					*/
/*			The Trapezoid Rule Apporximation			*/
/********************************************************/
/* with Romberg integration, the integral is evaluated	*/
/* increasing parameter, m, as in...					*/
/*														*/
/*						Xm - Xo							*/
/*					h = -------							*/
/*							m							*/
/*						  2								*/
/*														*/
/* This creates the first column of the romber table	*/
/* successive rows are created by subtracting elements  */
/* in the first row from one another					*/
/********************************************************/
/*					 k									*/
/*					4 * Tm+k,k-1 - Tm+k-1,k-1			*/
/*		Tm+k,k  = ----------------------------			*/
/*						k								*/
/*					   4	-	1						*/
/********************************************************/

double functionInput(double xValue)
{
	double functionVal = sqrt(1-pow(xValue,2);
	return functionVal;
}

double** RombergtrapezoidRule(double upperIntegralLimit,double lowerIntegralLimit)
{
	/****************************************************/
	/*				TRAPEZOID RULE						*/
	/*				 (composite)						*/
	/*		Xn			Fo						Fn		*/
	/*	INT  F(x) = h ( -- + F1 + F2...+Fn-1 +  -- )	*/
	/*		Xo			2						2		*/
	/****************************************************/
	double **rombergTable[11][11];
	double step=0;
	int numPoints=0;
	int rombergMValue=0;
	double integralValue=0;
	int counter=0;
	double xValue=0;
	//double upperIntegralLimit =  1;
	//double lowerIntegralLimit = -1;
	for (rombergMValue=0;rombergMValue<11;rombergMValue++)
	{
		integralValue=0;
		counter=1;
		numPoints = pow(2,m);
		step = (upperIntegralLimit - lowerIntegralLimit) / numpoints;
		for (xValue = lowerIntegralLimit; xValue < upperIntegralLimit + step; xValue += step)
		{
				if (counter==1 || counter==numPoints)
				{
					integralValue += (functionInput(xValue) / 2);
				}
				else
				{
					integralValue += functionInput(xValue);
				}
		}
		rombergTable[0][rombergMValue] = integralValue;
	}
}
int main()
{
	double** rombergTable[11][11];
	double upperIntegralLimit=1;
	double lowerIntegralLimit=-1
	int i;
	rombergTable = RombergtrapezoidRule(upperIntegralLimit,lowerIntegralLimit);
	for (i=0;i<11;i++)
	{
		printf("%lf",rombergTable[0][i];);
	}
	
}

