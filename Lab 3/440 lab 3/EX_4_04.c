/*
 *  EX_4_4.c
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/18/11.
 *  Copyright 2011 __Matthew Beck__. All rights reserved.
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double PI = 3.141592653589793238;
/****************************************************/
/*				Fresnal Equations					*/
/*						      2	          2			*/
/*		    I = .5 Io[(C(v)+.5) + (S(v)+.5)	]		*/
/*													*/
/*													*/
/*			C(v) = INT{0->v} cos(PI*w^2/2) dw		*/
/*													*/
/*			S(v) = INT{0->v} sin(PI*w^2/2) dw		*/
/****************************************************/
double functionInput(double omega, int fresnalCount)
{
	/************************************************************************/
	/*			Fresnal Equation (1)										*/
	/*				    1			   2									*/
	/*			F(w) = ---  cos( PI * w  )									*/
	/*					2	   (  ------ )									*/
	/*						   (    2 	 )									*/
	/*																		*/
	/************************************************************************/
	/************************************************************************/
	/*			Fresnal Equation (2)										*/
	/*					1       		2									*/
	/*			F(w) = ---   sin( PI * w  )									*/
	/*					2    	(  -----  )									*/
	/*					    	(    2	  )									*/
	/*																		*/
	/************************************************************************/
	double fresnalTwo=0;
	double fresnalOne=0;
	fresnalOne = cos(PI*pow(omega,2)/2) / 2;
	fresnalTwo = sin(PI*pow(omega,2)/2) / 2;
	if (fresnalCount==1)
	{
		return fresnalOne;
	}
	else
	{
		return fresnalTwo;
	}
}
/********************************************************************************************/
double** rombergTrapezoidRule(double upperIntegralLimit,double lowerIntegralLimit, int fresnalCount)
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
					integralValue += (functionInput(xValue,fresnalCount) / 2);
				}
				else
				{
					integralValue += functionInput(xValue,fresnalCount);
				}
				counter++;
		}
		integralValue*=step;
		rombergTable[0][rombergMValue] = integralValue;
		
	}
	return rombergTable;
}
/****************************************************************************/
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
	int row=0;
	for (k=1;k<maxRombergRow;k++)
	{
		for (m=0;m<maxRombergRow;m++)
		{
			row = m+k
			if (row == maxRombergRow)
			{
				break;
			}
			rombergTable[k][row] = (pow(4,k)*rombergTable[k-1][row] - rombergTable[k-1][row-1]) / (pow(4,k)-1);
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
	double quotient = 0;
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
/****************************************************************************/
int main()
{
	double** rombergTable; // in (column, row) form....Return for rombergTrapezoidRule and input for rombergTableExtrapolation()
	double* diagnosticTest; // return array for rombergDiagnostic();
	//double diagnosticCheck=0; // Average ratio of diagnostic to 4
	FILE *fresnalFile;
	double fresnalEquation[2]; // number of Fresnal Equations
	double Intensity = 0;
	int fresnalCount=1;
	double upperIntegralLimit=.1; // upper limit of integration
	double lowerIntegralLimit=0; // lower limit of integration
	int maxRombergRow = 11; //maximum extent of Romberg rows and columns
	int i=0; //counter to initiate romberg table
	int j=0; //counter to initiate romberg table
	fresnalFile = fopen("knifesEdge.txt", "w");
	rombergTable = (double**)malloc(maxRombergRow*sizeof(double));
	diagnosticTest = (double*)malloc((maxRombergRow-2)*sizeof(double));
	//dynamically allocate 2D array...
	for (i = 0;i < maxRombergRow; i++)
	{
		rombergTable[i] = (double*)malloc(maxRombergRow*sizeof(double));
	}
	// loop through increasing values of omega to get a handle on how intensity changes as
	// a function of distance.
	for (upperIntegralLimit = 0; upperIntegralLimit < 10; upperIntegralLimit += .01)
	{
		for (fresnalCount=1;fresnalCount<3;fresnalCount++)
		{
			for (i=0;i<maxRombergRow;i++)
			{
				for (j=0;j<maxRombergRow;j++)
				{
					rombergTable[i][j]=0;
				}
			}
		
			//perform the compostie trapezoid rule and return the first column of the romberg table.
			rombergTable = rombergTrapezoidRule(upperIntegralLimit,lowerIntegralLimit, fresnalCount);
			// run the diagnostic test on the first column.
			diagnosticTest = rombergDiagnostic(rombergTable,maxRombergRow);
			// total numerical value for all diagnotic values.
			for (i = 0; i < maxRombergRow-2; i++)
			{
				printf("m=%i\tRatio=%lf\n",i,diagnosticTest[i]);
			}
	
			rombergTable = rombergTableExtrapolation(rombergTable,maxRombergRow);
			for (i = 0; i < 11; i++)
			{
				printf("\n");
				for(j=0;j<11;j++)
				{
					printf("%lf\t",rombergTable[j][i]);
				}
			}
			fresnalEquation[fresnalCount-1] = rombergTable[10][10];
		}
		
		Intensity = (.5)*( (pow( (fresnalEquation[0] + .5), 2)) + (pow( (fresnalEquation[1] + .5), 2)));
		fprintf(fresnalFile, "%lf\t%lf\n",upperIntegralLimit,Intensity);
	
	
	}
	fclose(fresnalFile);
	return 0;

	
}

