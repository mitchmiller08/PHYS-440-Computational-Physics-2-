/*
 *  Ex_7_1_F.c
 *  Lab 4 Phys 440
 *
 *  Created by Matthew Beck on 4/3/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
double muValue (double xValue, double length)
{
	double mu0 = 1e-3; // kg / m
	double mu;
	double delta = 7e-3
	; // kg / m^2
	mu = mu0; // + delta*xValue;
	return mu;
}
int main()
{
	double uOld[200], u[200], uNew[200], u_0[200], du_0[200];
	double time, dt, dx, epsilon, T , mu, c, h, x, maxTime, L;
	int i, nx, nt, n;
	FILE* stringOutPut;
	stringOutPut = fopen("stringTrial2.txt", "w");
	n = 100;
	h = 1 / 100.;
	i = 0;
	T = 10; //N (kg * m / s^2
	L = 1; // meters
	mu = 1e-3; // grams
	dt = 1; // seconds
	maxTime = dt * 20;
	
	time = 0;
	while (i < n)
	{
		time = 0;
		x = h*(double)i;
		c = sqrt(T / muValue(x, L));
		u_0[i] = exp(-100.*pow((x  - .5),2));
		i++;
	}
	u_0[0] = 0;
	u_0[100] = 0;
	
	i=0;
	while (i < n)
	{
		time = 0;
		c = sqrt(T / muValue(x, L));
		x = h*(double)i;
		du_0[i] = 0.;
		i++;
	}

	i = 0;
	while ( i <  n )
	{
		uOld[i] = u_0[i];
		i++;
	}
	
	u[0] = 0;
	u[n] = 0;
	i = 1;
	while (i < n-1)
	{
		x = h*(double)i;
		c = sqrt(T / muValue(x, L));
		epsilon = pow((dt * c / h), 2);
		u[i] = .5 * epsilon * (u_0[i+1] + u_0[i-1]) + (1 - epsilon) * u_0[i] + dt * du_0[i];
		i++;
	}
	
	while (time < maxTime)
	{
		time += dt;
		uNew[0] = 0;
		uNew[n] = 0;
		i = 1;
		while (i < n-1 ) 
		{
			x = h*(double)i;
			c = sqrt(T / muValue(x, L));
			epsilon = pow((dt * c / h), 2);
			uNew[i] = epsilon * ( u[i+1] + u[i-1] ) + 2 * (1 - epsilon) * u[i] - uOld[i];
			i++;
		}
		i = 0;
		while (i < n)
		{
			uOld[i] = u[i]; 
			u[i] = uNew[i];
			i++;
		}
		i = 0;
		while (i < n)
		{
			x = h*(double)i;
			fprintf(stringOutPut, "%lf\t%lf\n",x,u[i]);
			i++;
		}
	}
	return 0;
	
}
