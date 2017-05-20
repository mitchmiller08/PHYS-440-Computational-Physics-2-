/*
 *  RungeKuttaFehlberg.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 3/19/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
double force(double time, double x, double velocity)
{
	double gravity = 9.81;
	double length = 9.81;
	double damping = 1 / 2.;
	double forceAmplitude = 1.2;
	double forceFrequency = (2 / 3.);
	double forceValue = -(gravity / length)*sin(x) - damping*velocity 
						+ forceAmplitude*sin(forceFrequency*time);
	return forceValue;
}
double *rungeKuttaFehlberg(double x0, double velocity, double time, double timeStep)
{
		double *kOrder;
		kOrder = (double*)malloc(4*sizeof(double));
		kOrder[0] = timeStep*force(time, x0, velocity);
		kOrder[1] = timeStep*force(time+.5*timeStep, x0 + .5*timeStep*velocity, velocity + .5*kOrder[0]);
		kOrder[2] = timeStep*force(time + .5*timeStep, x0 + .5*timeStep*velocity + .25*timeStep*kOrder[0], velocity + .5*kOrder[1]);
		kOrder[3] = timeStep*force(time + timeStep, x0 + timeStep*velocity + .5*timeStep*kOrder[1], velocity + kOrder[2]);
	
		
	return kOrder;
	
}
double main()
{
	double x0 = .2;
	double velocity = 0;
	double time = 0;
	double timeStep = .04;
	double *kOrder;
	double k1,k2,k3,k4;
	int counter  = 0;
	FILE *file;
	file = fopen("PendulumData.txt","w");
	kOrder = (double*)malloc(4*sizeof(double));
	while (counter < 1500)
	{
		kOrder = rungeKuttaFehlberg(x0,velocity,time, timeStep);
		k1 = kOrder[0];
		k2 = kOrder[1];
		k3 = kOrder[2];
		k4 = kOrder[3];
		x0 += timeStep*(velocity+(k1 + k2 + k3 + k4)/6);
		/*if (x0 < -3.14159)
		{
			x0 += 3.14159;
		}
		if (x0 > 3.14159)
		{
			x0 -= 3.14159;
		}*/
		velocity += (k1 + 2*k2 + 2*k3 + k4)/6;
		time += timeStep;
		fprintf(file, "%lf\t%lf\n",time,x0);
		counter++;
	}
}

