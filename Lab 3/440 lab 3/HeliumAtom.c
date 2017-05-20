/*
 *  HeliumAtom.c
 *  440 lab 3
 *
 *  Created by Matthew Beck on 3/4/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <rngs.h>
//#include <rvgs.h>
#include "rvgs.c"
#include "rngs.c"
//#include <rngs.c>
//#include <rvgs.c>

double bohrRadius = .5291772083; //in units of angstroms
 double Pi = 3.141592653589793;
double E = 2.718281828459045;
//double eCharge = −1.602176487;
// e^2 / 4Pi*epsilion ~~> eV*A
double eCharge = 14.39964439;


double one_S_State_Squared(double xValue, double yValue, double zValue)
{
	double functionVal = 0;
	double exponent = 0;
	double Z = 2;
	double radius = sqrt(pow(xValue,2)+pow(yValue,2)+pow(zValue,2));
	
	exponent = -2*Z*radius / bohrRadius;
	functionVal = pow(E,exponent);
	return functionVal;

}
double two_S_State_Squared(double xValue, double yValue, double zValue)
{
	double functionVal = 0;
	double exponent = 0;
	double coefficient = 0;
	double Z = 2;
	double power = 0;
	double radius = sqrt(pow(xValue,2)+pow(yValue,2)+pow(zValue,2));
	coefficient =  pow(1 - (.5)*Z*(radius / bohrRadius),2);
	exponent = -Z*radius / (bohrRadius);
	power = pow(E,exponent);
	functionVal = coefficient*power;
	return functionVal;
}
double two_P_State_Squared(double xValue, double yValue, double zValue)
{
	double functionVal = 0;
	double exponent = 0;
	double Z = 2;
	double power=0;
	double radius = sqrt(pow(xValue,2)+pow(yValue,2)+pow(zValue,2));
	double coefficient = pow(radius,2);
	exponent = -Z*radius / (bohrRadius);
	power = pow(E, exponent);
	functionVal = coefficient*power;
	return functionVal;
	
}
double K1s2s_integrand(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double functionVal = 0;
	double exponent = 0;
	double Z = 2;
	double radius1 = sqrt(pow(x1,2)+pow(y1,2)+pow(z1,2));
	double radius2 = sqrt(pow(x2,2)+pow(y2,2)+pow(z2,2));
	double coefficient = (1-(.5)*Z*radius2/bohrRadius)*(1-(.5)*Z*radius1/bohrRadius);
	exponent = -Z*(radius1+radius2)* (3/ (2*bohrRadius));
	functionVal = pow(E,exponent)*coefficient;
	return functionVal;
}
double K1s2pr1(double x1, double y1, double z1)
{
	double functionVal = 0;
	double exponent = 0;
	double Z = 2;
	double radius1 = sqrt(pow(x1,2)+pow(y1,2)+pow(z1,2));
	//double radius2 = sqrt(pow(x2,2)+pow(y2,2)+pow(z2,2));
	double coefficient = z1;
	exponent = -(3*Z*radius1/(2*bohrRadius));
	functionVal = coefficient*pow(E,exponent);
	return functionVal;
	
}
double K1s2pr2(double x2, double y2, double z2)
{
	double radius2 = sqrt(pow(x2,2)+pow(y2,2)+pow(z2,2));
	double Z=2;
	double coefficient = z2;
	double exponent = -(3*Z*radius2/(2*bohrRadius));
	double functionVal = coefficient*pow(E,exponent);
	return functionVal;
}
double two_S_State(double xValue, double yValue, double zValue)
{
	double functionVal = 0;
	double exponent = 0;
	double coefficient = 0;
	double Z = 2;
	double power = 0;
	double radius = sqrt(pow(xValue,2)+pow(yValue,2)+pow(zValue,2));
	coefficient =  (1 - (.5)*Z*(radius / bohrRadius));
	exponent = -Z*radius / 2*(bohrRadius);
	power = pow(E,exponent);
	functionVal = coefficient*power;
	return functionVal;
}

double repulsion(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double functionVal = 0;
	double radius12 = sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));

	//functionVal = pow(eCharge,2) / fabs(radius1 - radius2);
	functionVal = 1 / fabs(radius12);
	return functionVal; 
}
double J1s1s(double upperlimit, double lowerLimit)
{
	// coordinates of electron 1
	double x1=0;
	double y1=0;
	double z1=0;
	// coordinates of electron 2
	double x2=0;
	double y2=0;
	double z2=0;
	double oneS=0;
	double twoS = 0;
	double repulsionTerm = 0;
	double functionVal = 0;
	double Z = 2;
	double width = upperlimit-lowerLimit;
	int maxIterations = 1000000;
	int iteration = 0;
	double coefficient1S = pow((Z / bohrRadius),3) / Pi;
	double coefficient2S = pow((Z / bohrRadius),3) / (8*Pi);
	for (iteration = 0; iteration<maxIterations; iteration++)
	{
		x1 = Uniform(lowerLimit,upperlimit);
		y1 =  Uniform(lowerLimit, upperlimit);
		z1 =  Uniform(lowerLimit, upperlimit);
		x2 =  Uniform(lowerLimit, upperlimit);
		y2 =  Uniform(lowerLimit, upperlimit);
		z2 =  Uniform(lowerLimit, upperlimit);
		oneS = fabs(one_S_State_Squared(x1,y1,z1));
		twoS = fabs(two_S_State_Squared(x2,y2,z2));
		repulsionTerm = repulsion(x1,y1,z1,x2,y2,z2);
		functionVal += oneS*twoS*repulsionTerm;
	}
	
	functionVal /= maxIterations;
	functionVal *= coefficient1S*coefficient2S*eCharge*pow(width,6);
	return functionVal;
	
}
double K1s2s(double upperlimit, double lowerLimit)
{
	// coordinates of electron 1
	double x1=0;
	double y1=0;
	double z1=0;
	// coordinates of electron 2
	double x2=0;
	double y2=0;
	double z2=0;
	double Z = 2;
	double oneS1=0;
	double oneS2=0;
	double twoS1 = 0;
	double twoS2=0;
	double repulsionTerm = 0;
	double kIntegrand=0;
	double functionVal = 0;
	double width = upperlimit - lowerLimit;
	double zBohr = Z / bohrRadius;
	double coefficient1S = (1 / (8*pow(Pi,2)))*pow(zBohr,6);
	//double coefficient2S = pow((Z / bohrRadius), 3/2) / sqrt(8*Pi);
	int maxIterations = 1000000;
	int iteration = 0;
	for (iteration = 0; iteration<maxIterations; iteration++)
	{
		x1 = Uniform(lowerLimit, upperlimit);
		y1 = Uniform(lowerLimit, upperlimit);
		z1 = Uniform(lowerLimit, upperlimit);
		x2 = Uniform(lowerLimit, upperlimit);
		y2 = Uniform(lowerLimit, upperlimit);
		z2 = Uniform(lowerLimit, upperlimit);
		/*oneS1 = one_S_State(x1,y1,z1);
		oneS2 = one_S_State(x2,y2,z2);
		twoS1 = two_S_State(x1,y1,z1);
		twoS2 = two_S_State(x2,y2,z2);*/
		kIntegrand = K1s2s_integrand(x1,y1,z1,x2,y2,z2);
		repulsionTerm = repulsion(x1,y1,z1,x2,y2,z2);
		functionVal += kIntegrand*repulsionTerm;
	}
	functionVal /= maxIterations;
	functionVal *= eCharge*pow(width,6)*coefficient1S;
	return functionVal;
}
double J1s2P(double upperLimit, double lowerLimit)
{
	// coordinates of electron 1
	double x1=0;
	double y1=0;
	double z1=0;
	// coordinates of electron 2
	double x2=0;
	double y2=0;
	double z2=0;
	double Z = 2;
	double functionVal = 0;
	double repulsionTerm = 0;
	double width = upperLimit - lowerLimit;
	double zBohr = Z / bohrRadius;
	int maxIterations = 1000000;
	int iteration = 0;
	double oneS = 0;
	double twoP = 0;
	double coefficient1S = pow(zBohr,3) / Pi;
	double coefficient2P = (1 / (64*pow(Pi,2)))*pow(zBohr,8);
	for (iteration = 0; iteration<maxIterations; iteration++)
	{
		x1 = Uniform(lowerLimit, upperLimit);
		y1 = Uniform(lowerLimit, upperLimit);
		z1 = Uniform(lowerLimit, upperLimit);
		x2 = Uniform(lowerLimit, upperLimit);
		y2 = Uniform(lowerLimit, upperLimit);
		z2 = Uniform(lowerLimit, upperLimit);
		oneS = fabs(one_S_State_Squared(x1, y1, z1));
		twoP = fabs(two_P_State_Squared(x2, y2, z2));
		repulsionTerm = repulsion(x1, y1, z1, x2, y2, z2);
		functionVal += oneS*twoP*repulsionTerm;
	}
	functionVal /= maxIterations;
	functionVal *= coefficient2P*eCharge*pow(width,6);
	return functionVal;
}
double K1s2p(double upperLimit, double lowerLimit)
{
	
	double x1=0;
	double y1=0;
	double z1=0;
	// coordinates of electron 2
	double x2=0;
	double y2=0;
	double z2=0;
	double Z = 2;
	/*double oneS1=0;
	double oneS2=0;
	double twoS1 = 0;
	double twoS2=0;*/
	double repulsionTerm = 0;
	double rOneTerm = 0;
	double rTwoTerm = 0;
	double functionVal = 0;
	double width = upperLimit - lowerLimit;
	double zBohr = Z / bohrRadius;
	double coefficient2P = (1 / (32*pow(Pi,2)))*pow(zBohr,8);
	//double coefficient2S = pow((Z / bohrRadius), 3/2) / sqrt(8*Pi);
	int maxIterations = 1000000;
	int iteration = 0;
	for (iteration = 0; iteration<maxIterations; iteration++)
	{
		x1 = Uniform(lowerLimit, upperLimit);
		y1 = Uniform(lowerLimit, upperLimit);
		z1 = Uniform(lowerLimit, upperLimit);
		x2 = Uniform(lowerLimit, upperLimit);
		y2 = Uniform(lowerLimit, upperLimit);
		z2 = Uniform(lowerLimit, upperLimit);
		/*oneS1 = one_S_State(x1,y1,z1);
		oneS2 = one_S_State(x2,y2,z2);
		twoS1 = two_S_State(x1,y1,z1);
		twoS2 = two_S_State(x2,y2,z2);*/
		rOneTerm = K1s2pr1(x1,y1,z1);
		rTwoTerm = K1s2pr2(x2,y2,z2);
		repulsionTerm = repulsion(x1,y1,z1,x2,y2,z2);
		functionVal += rOneTerm*repulsionTerm*rTwoTerm;
	}
	functionVal /= maxIterations;
	functionVal *= eCharge*pow(width,6)*coefficient2P;
	return functionVal;

}
int main()
{
	double upperLimit = 5;
	double lowerLimit = -5;
	double j1s2sTerm = 0;
	double k1s2sTerm = 0;
	double j1s2pTerm = 0;
	double k1s2pTerm = 0;
	int j=0;
	int k=0;
	FILE *file;
	double term=0;
	int histogram[2000];
	file = fopen("HeliumAtom.txt", "w");
	for (k=0;k<2000;k++)
	{
		histogram[k] = 0;
	}
	for (j=0;j<1000;j++)
	{
		j1s2sTerm += J1s1s(upperLimit, lowerLimit);
		k1s2sTerm += K1s2s(upperLimit, lowerLimit);
		j1s2pTerm += J1s2P(upperLimit, lowerLimit);
		k1s2pTerm += K1s2p(upperLimit, lowerLimit);
		/*for (k=0;k<2000;k++)
		{
			term = (double)k/100;
			if (jTerm<term)
			{
				histogram[k]++;
				break;
			}
		}*/
	
	
	}
	j1s2sTerm /= 1000;
	k1s2sTerm /= 1000;
	j1s2pTerm /= 1000;
	for (k=0;k<2000;k++)
	{
		term = (double)k/100;
		fprintf(file,"%lf\t%i\n",term, histogram[k]);
	}
	fclose(file);
	return 0;
}

