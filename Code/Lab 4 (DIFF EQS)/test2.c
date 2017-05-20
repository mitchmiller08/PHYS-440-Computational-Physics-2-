/*
 *  QuantumMechanics.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 3/30/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
double *kOrder;
double hBarSquaredOverMass = 7.630517543; //eV*A^2
int choice;
double Potential(double x)
{
double potential = 0;
	if (choice==1)
	{
		potential = fabs(x);
	}
	else
	{
		potential = x*x;
	}
	return potential;
}
double d2PsidX2(double x, double Psi, double dPsidX, double energy)
{
	double forceValue = 0;
	forceValue = -(2/hBarSquaredOverMass)*(energy-Potential(x))*Psi;
	return forceValue;
}
void *rungeKutta(double Psi0, double dPsidX, double x0, double xStep, double energy)
{
		/*double *kOrder;
		kOrder = (double*)malloc(4*sizeof(double));*/
		kOrder[0] = xStep*d2PsidX2(x0,                Psi0,                                                 dPsidX,								 energy);
		kOrder[1] = xStep*d2PsidX2(x0 + .5*xStep,	  Psi0 + .5*xStep*dPsidX,                               dPsidX + .5*kOrder[0],				 energy);
		kOrder[2] = xStep*d2PsidX2(x0 + .5*xStep,     Psi0 + .5*xStep*dPsidX + .25*xStep*kOrder[0],         dPsidX + .5*kOrder[1],				 energy);
		kOrder[3] = xStep*d2PsidX2(x0 + xStep,        Psi0 + xStep*dPsidX + .5*xStep*kOrder[1],             dPsidX + kOrder[2],					 energy);
	
		
	return kOrder;
	
}
double findEnergy(double oldEnergy,double currentEnergy,double x0,double xInitial,double xFinal,double oldPsi,double dPsidX,double xStep)
{
	double midPointEnergy = 0;
	double Psi0 = 0;
	double k1,k2,k3,k4;
	int repetitionCounter = 0;
	double difference=fabs(oldEnergy-currentEnergy);
	while (difference > 5e-16)
	{
		x0 = xInitial;
		Psi0 = 0;
		dPsidX = 1;
		midPointEnergy = (currentEnergy+oldEnergy)/2;
		if (midPointEnergy == oldEnergy)
		{
			midPointEnergy+=1e-16;
		}
		while(x0<-xInitial)
		{
			rungeKutta(Psi0,dPsidX,x0,xStep,midPointEnergy);
			k1 = kOrder[0];
			k2 = kOrder[1];
			k3 = kOrder[2];
			k4 = kOrder[3];
			Psi0 += xStep * (dPsidX + (k1 + k2 + k3)/6);
			dPsidX += (k1 + 2*k2 + 2*k3 + k4)/6;
			x0 += xStep;
		}
		if(Psi0*oldPsi>0)
		{
			oldEnergy = midPointEnergy;
		}
		else
		{
			currentEnergy = midPointEnergy;
		}
		difference = fabs(oldEnergy-currentEnergy);
		oldPsi = Psi0;
		repetitionCounter++;
		if (repetitionCounter>100)
		{
			break;
		}
	}
	return currentEnergy;
}
double outPutWaveFunctions(double energy,double xInitial,double xFinal,int counter,double xStep)
{
	double k1,k2,k3,k4;
	double x0=xInitial;
	double Psi0 = 0;
	double dPsidX = 1;
	double coefficient = 0;
	double integralWeight = 0;
	int integrationCounter = 0;
	int integrationCounterMax = (int)(xFinal+xInitial)/xStep;
	double integralPsiSquared = 0;
	char fileNameWrite[50];
	char fileNameWriteNormal[50];
	char fileNameRead[50];
	sprintf(fileNameWrite, "waveFunctionData%i.txt",counter);
	
	FILE *fileWrite = fopen(fileNameWrite,"w");
	while(x0<xFinal)
	{
		fprintf(fileWrite,"%lf\t%lf\n",x0,Psi0);
		rungeKutta(Psi0,dPsidX,x0,xStep,energy);
		k1 = kOrder[0];
		k2 = kOrder[1];
		k3 = kOrder[2];
		k4 = kOrder[3];
		Psi0 += xStep * (dPsidX + (k1 + k2 + k3)/6);
		if (integrationCounter == 0 || integrationCounter == integrationCounterMax-1)
		{
			integralWeight = .5;
		}
		else
		{
			integralWeight = 1;
		}
		integralPsiSquared += xStep*integralWeight * pow(Psi0,2);
		dPsidX += (k1 + 2*k2 + 2*k3 + k4)/6;
		x0 += xStep;
	}
	xFinal = -xInitial;
	x0=xFinal;
	Psi0 = 0;
	dPsidX = 1;
	/*while (x0 > 0)
	{
		fprintf(fileWrite,"%lf\t%0.32f\n",x0,Psi0);
		rungeKutta(Psi0,dPsidX,x0,xStep,energy);
		k1 = kOrder[0];
		k2 = kOrder[1];
		k3 = kOrder[2];
		k4 = kOrder[3];
		Psi0 += xStep * (dPsidX + (k1 + k2 + k3)/6);
		if (integrationCounter == 0 || integrationCounter == integrationCounterMax-1)
		{
			integralWeight = .5;
		}
		else
		{
			integralWeight = 1;
		}
		integralPsiSquared += xStep*integralWeight * pow(Psi0,2);
		dPsidX += (k1 + 2*k2 + 2*k3 + k4)/6;
		x0 -= xStep;
		if (x0<.040)
		{
			double eat = 1;
		}
	}*/
	coefficient = sqrt((1/integralPsiSquared));
	fclose(fileWrite);
	sprintf(fileNameRead, "waveFunctionData%i.txt",counter);
	sprintf(fileNameWriteNormal, "waveFunctionNormal%i.txt",counter);
	FILE *fileWriteNormal = fopen(fileNameWriteNormal,"w");
	FILE *fileRead = fopen(fileNameRead,"r");
	while ( fscanf(fileRead,"%lf	%lf",&x0,&Psi0) != EOF)
	{
		Psi0*=coefficient;
		fprintf(fileWriteNormal,"%lf\t%0.16f\n",x0,Psi0+energy);
	}	
		return 1;
}
int main()
{
	double xInitial = -10; //in Angstroms
	double xFinal = 10;
	double x0 = 0;
	double psiInitial = 0;
	double dPsiInitial = 1;
	double dPsidX = 0; //guess for slope
	double Psi0 = 0;
	double k1,k2,k3,k4;
	double xStep = 1e-3;
	double energy = 0;
	double oldEnergy = 0;
	double energyStep = .01;
	double midPointEnergy = 0;
	double oldPsi = 0;
	double Pi = 3.141592653589793238;
	int counter = 0;
	printf("Choose a potential: [1]V(x) = |x|	[2]V(x) = 5|x| [-inf,0], V(x) = x^2 [0,inf]");
	scanf("%i",&choice);
	FILE *file = fopen("potentialData.txt","w");
	kOrder = (double*)malloc(4*sizeof(double));
	while (counter < 5)
	{
		x0 = xInitial;
		Psi0 = psiInitial;
		dPsidX = dPsiInitial;
		while(x0<-xInitial)
		{
			rungeKutta(Psi0,dPsidX,x0,xStep,energy);
			k1 = kOrder[0];
			k2 = kOrder[1];
			k3 = kOrder[2];
			k4 = kOrder[3];
			Psi0 += xStep * (dPsidX + (k1 + k2 + k3)/6);
			dPsidX += (k1 + 2*k2 + 2*k3 + k4)/6;
			x0 += xStep;
		}
		//printf("%i\t%lf\n",counter,Psi0);
		if (Psi0*oldPsi<0)
		{
			energy = findEnergy(oldEnergy,energy,xInitial,xInitial,xFinal,oldPsi,dPsiInitial,xStep);
			outPutWaveFunctions(energy,xInitial,xFinal,counter,xStep);
			printf("%i\t%lf\n",counter,Psi0);
			counter++;
		}
		oldPsi = Psi0;
		oldEnergy = energy;
		energy+=energyStep;
	}
	
	for (x0=xInitial;x0<xFinal;x0+=xStep)
	{
		fprintf(file,"%lf\t%lf\n",x0,Potential(x0));
	}
	return 0;

}


