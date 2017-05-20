/*
 *  EX_5_47.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 4/4/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "clapack.h"
//#include "cblas.h"
double* alphaSolutionMatrix;
static long systemSolver(long N, long NRHS, double *A, long LDA, long *IPIV, double *B,long LDB)
{
	int counter = 0;
	
	extern void dgesv_(const long *N, const long *NRHSp, double *A, const long *LDA, long *IPIV, double *B, const long *LDB, long *INFO);

	long info;
	dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &info);
	while (counter < (int)N)
	{
		alphaSolutionMatrix[counter] = A[counter];
		counter++;
	}
	return info;
}
double* alphaMatrixInitialize(double* alphaMatrix, int numIntervals)
{
	int counter = 0;
	while (counter < numIntervals)
	{
		alphaMatrix[counter] = 1;
		counter++;
	}
	return alphaMatrix;
}
double* RHSMatrixInitialize(double* RHSMatrix, double* alphaMatrix, int numIntervals, double hStep)
{
	double xValue = 0;
	int counter = 0;
	double value = 0;
	while (counter < numIntervals)
	{
		xValue	= hStep*(double)counter;
		if (counter == 0)
		{
			value = 6*pow(hStep,2)*xValue;
			RHSMatrix[counter] = value;
		}
		if (counter == numIntervals)
		{
			value = 6*pow(hStep,2) - 1;
		}
		else
		{
			value = 6*pow(hStep,2)*xValue;
			RHSMatrix[counter] = value;
		}
		counter++;
	}
	return RHSMatrix;
}
double* LHSMatrixInitiate(double* LHSMatrix, int numIntervals)
{
	double a = -2;
	double b = 1;
	double c = 1;
	int aInput = 0;
	int bInput = 0;
	int cInput = 0;
	int Counter = 0;
	while (Counter < numIntervals)
	{
		LHSMatrix[Counter] = 0;
		Counter++;
	}
	Counter = 0;
	Counter = 0;
	while (Counter < (int)sqrt((double)numIntervals))
	{
		aInput = 10*Counter;
		bInput = 10*Counter+1;
		LHSMatrix[aInput] = a;
		LHSMatrix[bInput] = b;
		Counter++;
	}
	Counter = 1;
	while (Counter < (int)sqrt((double)numIntervals))
	{
		cInput = 10*Counter-1;
		LHSMatrix[cInput] = c;
		Counter++;
	}
	Counter = 0;
	while (Counter < numIntervals) 
	{
		if (Counter > 0 && Counter%9 == 0)
		{
			printf("\n");
		}
		printf("%lf\t",LHSMatrix[Counter]);
		Counter++;
	}
	return LHSMatrix;
}
double basisFunctionCalc(int numIntervals, double xValue, int interval)
{
	double xStep = 1 / (double)numIntervals;
	double xMax = (double)(interval+1)*xStep;
	double xMin = (double)(interval-1)*xStep;
	double xMiddle = (double)(interval)*xStep;
	double maxXValue = 1;
	double basisFunction = 0;
	
	if (xValue < xMin)
	{
		basisFunction = 0;
	}
	if (xValue >= xMin && xValue < xMiddle)
	{
		basisFunction = (xValue - xMin) / (xMiddle - xMin);
	}
	if (xValue >= xMiddle && xValue < xMax)
	{
		basisFunction = (xMax - xValue) / (xMax - xMiddle);
	}
	if (xValue >= xMax)
	{
		basisFunction = 0;
	}
	/*if (interval == numIntervals-1)
	{
		basisFunction = 1;
	}*/
	return basisFunction;
}
int main()
{
	int numIntervals = 10;
	double** basisFunction;
	double xIntervalStep = 1 / (double)numIntervals;
	double xValue = 0;
	double hStep = 1e-3;
	long* pivotMatrix;
	double* alphaMatrix;
	double* RHSMatrix;
	double* LHSMatrix;
	double* functionForm;
	int numBasisValues = 1 / hStep;
	int interval = 0;
	int basisCounter = 0;
	FILE* basisOutPut;
	basisOutPut = fopen("basisFunction.txt","w");
	functionForm = malloc(numBasisValues*sizeof(double));
	basisFunction = (double**)malloc((numIntervals+1)*sizeof(double*));
	while (interval < numIntervals+1)
	{
		basisFunction[interval] = (double*)malloc(numBasisValues*sizeof(double));
		interval++;
	}
	interval = 1;
	while (interval < numIntervals+1)
	{
		xValue = 0;
		basisCounter = 0;
		while (basisCounter < numBasisValues)
		{
			xValue = hStep*(double)basisCounter;
			basisFunction[interval][basisCounter] = basisFunctionCalc(numIntervals, xValue, interval);
			fprintf(basisOutPut,"%lf\t%lf\n", xValue, basisFunction[interval][basisCounter]);
			basisCounter++;
			
		}
		interval++;
	}
	
	/*************************************************************************************/
	/* Solving for the alpha's
	/************************************************************************************/
	int dgesv_INFO = 5;
	pivotMatrix = malloc((numIntervals-1)*sizeof(long));
	alphaMatrix = malloc((numIntervals-1)*sizeof(double));
	RHSMatrix = malloc((numIntervals-1)*sizeof(double));
	LHSMatrix = malloc((numIntervals-1)*(numIntervals-1)*sizeof(double));
	alphaSolutionMatrix = malloc((numIntervals-1)*sizeof(double));
	alphaMatrix = alphaMatrixInitialize(alphaMatrix,numIntervals-1);
	RHSMatrix = RHSMatrixInitialize(RHSMatrix,alphaMatrix,numIntervals-1,hStep);
	LHSMatrix = LHSMatrixInitiate(LHSMatrix,(numIntervals-1)*(numIntervals-1));
	
	dgesv_INFO = systemSolver((numIntervals-1), 1, LHSMatrix,(numIntervals-1), pivotMatrix, RHSMatrix, (numIntervals-1));
	interval = 0;
	while (interval < numIntervals-1 ) 
	{
		printf("%e\n",alphaSolutionMatrix[interval]);
		interval++;
	}
	interval = 0;
	while (interval < numBasisValues)
	{
		basisFunction[0][interval] *= 0;
		basisFunction[10][interval] *=1;
		interval++;
	}
	interval = 1;
	while (interval < numIntervals -1)
	{
		basisCounter =0;
		while (basisCounter < numBasisValues)
		{
			functionForm[basisCounter] = 0;
			basisFunction[interval][basisCounter] *= alphaSolutionMatrix[interval-1];
			basisCounter++;
		}
		interval++;
	}
	interval = 0;
	basisCounter = 0;
	while (interval < numIntervals)
	{
		basisCounter = 0;
		while (basisCounter < numBasisValues)
		{
			functionForm[basisCounter] += basisFunction[interval][basisCounter];
			basisCounter++;
		}
		interval++;
	}
	basisCounter = 0;
	while (basisCounter < numBasisValues)
	{
		xValue = hStep*(double)basisCounter;
		fprintf(basisOutPut,"%lf\t%lf\n", xValue, functionForm[basisCounter]);
		basisCounter++;
	}
	fclose(basisOutPut);
	
	
	
	
	return 0;
	
	
	
}

