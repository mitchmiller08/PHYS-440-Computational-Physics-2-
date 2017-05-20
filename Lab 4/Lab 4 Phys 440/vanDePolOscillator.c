/*
 *  vanDePolOscillator.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 3/21/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
double force(double time, double x, double velocity, double forceAmplitude)
{
	double eps = 1;
	double forceValue = 0;
	forceValue = -x - eps*(pow(x,2)-1)*velocity;
	return forceValue;
}
double *rungeKuttaFehlberg(double x0, double velocity, double time, double timeStep, double forceAmplitude)
{
		double *kOrder;
		double x,y,a0,a1,a2,a3,a4,a5,b10,b20,b30,b21,b31,b32,c0,c1,c2,c3,c4,c5,d0,d1,d2,d3,d4,d5,yHat;
		a0 = .25;
		
		kOrder = (double*)malloc(4*sizeof(double));
		kOrder[0] = timeStep*force(time,                x0,                                                 velocity,				 forceAmplitude);
		kOrder[1] = timeStep*force(time + .5*timeStep,  x0 + .5*timeStep*velocity,                          velocity ,				 forceAmplitude);
		kOrder[2] = timeStep*force(time + .5*timeStep,  x0 + .5*timeStep*velocity + .25*timeStep*kOrder[0], velocity,				 forceAmplitude);
		kOrder[3] = timeStep*force(time + timeStep,     x0 + timeStep*velocity + .5*timeStep*kOrder[1],     velocity,				 forceAmplitude);
	
		
	return kOrder;
	
}

int main()
{
	double *kOrder;
	double x0 = 3;
	double velocity = 0;
	double k1,k2,k3,k4;
	double time = 0;
	int N = 150;
	double forceAmplitude = .68;
	double Pi = 3.141592653589793238;
	double T = 2*Pi;
	double timeStep = T / N;
	int counter = 0;
	FILE *Attractor;
	FILE *timePlot;
	timePlot = fopen("timePlotVanDePol.txt","w");
	Attractor = fopen("vanDePolAttractor.txt","w");
	kOrder = (double*)malloc(4*sizeof(double));
	while (counter < 1000)
		{
			kOrder = rungeKuttaFehlberg(x0,velocity,time,timeStep, forceAmplitude);
			k1 = kOrder[0];
			k2 = kOrder[1];
			k3 = kOrder[2];
			k4 = kOrder[3];
			x0 += timeStep * (velocity + (k1 + k2 + k3)/6);
			if (x0 < -1*Pi)
			{
				x0 += 2*Pi;
			}
			if (x0 > Pi)
			{
				x0 -= 2*Pi;
			}
			velocity += (k1 + 2*k2 + 2*k3 + k4)/6;
			time += timeStep;
			fprintf(timePlot, "%lf\t%lf\n",x0,velocity);
			//counterOffSet = counter - 75;
				/*if (counter % 150 == 0)
				{
					fprintf(Attractor,"%.15f\t%.15f\t%.15f\n",time,x0,velocity);
				}*/
			counter++;
		}
	
	return 0;

}
