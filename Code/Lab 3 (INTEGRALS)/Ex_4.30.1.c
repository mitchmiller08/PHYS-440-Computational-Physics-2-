#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 3
// Exc. 4.30 2/28/12
// Random Walk Simulation

int main(){

	double timeStep;
	int counter=0,steps;
	double x,y,z;
	double vx,vy,vz,v0;
	double theta,phi;
	double distance;

	FILE *outfile;
	outfile = fopen("Ex_4.30.dat","w");

	theta = pi * drand48();
	phi = 2.0 * pi * drand48();
	v0 = 10;

	vx = v0 * sin(theta) * cos(phi);
	vy = v0 * sin(theta) * sin(phi);
	vz = v0 * cos(theta);
	x = 0;
	y = 0;
	z = 0;
	distance = sqrt(x*x+y*y+z*z);

	timeStep = 0.1;
	steps = 10000;

	fprintf(outfile,"%i\t%lf\n",counter,distance);

	for(counter=0;counter<steps;counter++){

		x = x + vx*timeStep;
		y = y + vy*timeStep;
		z = z + vz*timeStep;
		distance = sqrt(x*x+y*y+z*z);

		fprintf(outfile,"%i\t%lf\n",counter,distance);

		theta = pi * drand48();
		phi = 2.0 * pi * drand48();

		vx = v0 * sin(theta) * cos(phi);
		vy = v0 * sin(theta) * sin(phi);
		vz = v0 * cos(theta);

	}

}
