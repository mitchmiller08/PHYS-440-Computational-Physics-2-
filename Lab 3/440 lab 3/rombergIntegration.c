/*
 *  rombergIntegration.c
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
 #include <stdio.h>
 #include <stdlib.h>
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

double functionInput(double xValue, int choice)
{
	double functionVal;
	if (choice==1)
	{	//function without variable change
		functionVal = sqrt((1-pow(xValue,2)));
	}
	else
	{	//function with variable change
		// x=cos(0)
		functionVal = pow(sin(xValue),2);
	}
	return functionVal;
}
/*********************************************************************************/
double** rombergTrapezoidRule(double upperIntegralLimit,double lowerIntegralLimit, int choice)
{
	/****************************************************/
	/*				TRAPEZOID RULE						*/
	/*				 (composite)						*/
	/*		Xn			Fo						Fn		*/
	/*	INT  F(x) = h ( -- + F1 + F2...+Fn-1 +  -- )	*/
	/*		Xo			2						2		*/
	/****************************************************/
	double** rombergTable;
	int maxRombergRow = 11;
	double step=0;
	int numPoints=0;
	int rombergMValue=0;
	double integralValue=0;
	int counter=0;
	int i=0;
	int j=0;
	double xValue=0;
	rombergTable = (double**)malloc(maxRombergRow*sizeof(double));
	for (i=0;i<maxRombergRow;i++)
	{
		rombergTable[i] = (double*)malloc(maxRombergRow*sizeof(double));
	}
	for (i=0;i<maxRombergRow;i++)
	{
		for (j=0;j<maxRombergRow;j++)
		{
			rombergTable[i][j]=0;
		}
	}
	for (rombergMValue=0;rombergMValue<11;rombergMValue++)
	{
		integralValue=0;
		counter=0;
		numPoints = pow(2,rombergMValue);
		step = (upperIntegralLimit - lowerIntegralLimit) / numPoints;
		for (xValue = lowerIntegralLimit; xValue < upperIntegralLimit + step; xValue += step)
		{
				
				if (counter==0 || counter==numPoints)
				{
					integralValue += (functionInput(xValue,choice) / 2);
				}
				else
				{
					integralValue += functionInput(xValue,choice);
				}
				counter++;
		}
		integralValue*=step;
		rombergTable[0][rombergMValue] = integralValue;
		
	}
	return rombergTable;
}
/**************************************************************************/
double** rombergTableExtrapolation(double** rombergTable, int maxRombergRow)
{
	/********************************************************/
	/*					 k									*/
	/*					4 * Tm+k,k-1 - Tm+k-1,k-1			*/
	/*		Tm+k,k  = ----------------------------			*/
	/*						k								*/
	/*					   4	-	1						*/
	/********************************************************/
	int k=0;
	int m=0;
	int i=0;
	int row=0;
	int j=0;
	for (k=1;k<maxRombergRow;k++)
	{
		
		for (m=0;m<maxRombergRow;m++)
		{
			row = m+k;
			if (row == maxRombergRow)
			{
				break;
			}
			rombergTable[k][row] = (pow(4,k)*rombergTable[k-1][row] - rombergTable[k-1][row-1]) / (pow(4,k)-1);
			
		}
	}
	for (i=0;i<maxRombergRow;i++)
	{
		printf("\n");
		for (j=0;j<maxRombergRow;j++)
		{
			//rombergTable[i][j]=0;
			//printf("%lf\t", rombergTable[j][i]);
		}
	}
	
	return rombergTable;
}
/****************************************************************************/
double* rombergDiagnostic(double** rombergTable, int maxRombergRow)
{
	/********************************************************************/
	/*						Tm-1,0 - Tm,0								*/
	/*				Rm = ------------------	~~~~~~ > 4					*/
	/*						Tm,0 - Tm+1,0								*/
	/********************************************************************/
	
	double* ratioDiagnostic; // for romberg to work, must be ~~> 4
	int i=0;
	double numerator=0;
	double denominator=0;
	double forward=0;
	double quotient=0;
	double middle=0;
	double back=0;
	ratioDiagnostic = (double*)malloc((maxRombergRow-2)*sizeof(double));
	for (i=1;i<maxRombergRow-1;i++)
	{
		forward = rombergTable[0][i+1];
		middle = rombergTable[0][i];
		back = rombergTable[0][i-1];
		numerator = back-middle;
		denominator = middle-forward;
		quotient = numerator / denominator;
		ratioDiagnostic[i-1] = numerator / denominator;
	}
	return ratioDiagnostic;
}
int main()
{
	double** rombergTable; // in (column, row) form....Return for rombergTrapezoidRule and input for rombergTableExtrapolation()
	double* diagnosticTest; // return array for rombergDiagnostic();
	double diagnosticCheck=0; // Average ratio of diagnostic to 4
	int choice=1;
	double upperIntegralLimit=3.14159; // upper limit of integration
	double lowerIntegralLimit=0; // lower limit of integration
	int maxRombergRow = 11; //maximum extent of Romberg rows and columns
	int i=0; //counter to initiate romberg table
	int j=0; //counter to initiate romberg table
	rombergTable = (double**)malloc(maxRombergRow*sizeof(double));
	diagnosticTest = (double*)malloc((maxRombergRow-2)*sizeof(double));
	//dynamically allocate 2D array...
	for (i=0;i<maxRombergRow;i++)
	{
		rombergTable[i] = (double*)malloc(maxRombergRow*sizeof(double));
	}
	for (i=0;i<maxRombergRow;i++)
	{
		printf("\n");
		for (j=0;j<maxRombergRow;j++)
		{
			rombergTable[j][i]=0;
			//printf("%lf\t", rombergTable[j][i]);
		}
	}
	printf("Select a function to integrate\n\n [1] for Ex 4-5 without variable change \n [2] Ex 4-5 with variable change \n");
	scanf("%i",&choice);
	if (choice==1)
	{
		upperIntegralLimit=1;
		lowerIntegralLimit=-1;
	}
	else
	{
		upperIntegralLimit=3.14159265;
		lowerIntegralLimit=0;
	}
	//perform the compostie trapezoid rule and return the first column of the romberg table.
	rombergTable = rombergTrapezoidRule(upperIntegralLimit,lowerIntegralLimit, choice);
	// run the diagnostic test on the first column.
	diagnosticTest = rombergDiagnostic(rombergTable,maxRombergRow);
	// total numerical value for all diagnotic values.
	for (i=0;i<maxRombergRow-2;i++)
	{
		printf("m=%i\tRatio=%lf\n",i+1,diagnosticTest[i]);
	}
	/*/ Average value of all the diagnostic checks for different Romberg m values
	diagnosticCheck /= (maxRombergRow-2);
	// Ratio of the average against the # 4
	diagnosticCheck /= 4;
	// |(1-ratio)| < tolerance level
	diagnosticCheck -= 1;
	diagnosticCheck = fabs(diagnosticCheck);
	// tolerance level check
	/*if (diagnosticCheck>.1)
	{
		printf("Diagnostic Check = %lf \n Integration failed, perform change of variables",diagnosticCheck);
		return 1;
	}*/
	//else
	//{
		rombergTable = rombergTableExtrapolation(rombergTable,maxRombergRow);
		for (i=0;i<11;i++)
		{
			printf("\n");
			for(j=0;j<11;j++)
			{
				printf("%lf\t",rombergTable[j][i]);
			}
		}
	//}
	return 0;
}




