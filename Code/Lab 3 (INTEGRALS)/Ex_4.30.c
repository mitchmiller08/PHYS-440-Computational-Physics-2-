#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 3
// Exc. 4.30 2/28/12
// Random Walk Simulation

typedef struct{
	double a;
	double b;
	double c;
} point;

void step(point * pos, point * vel, double t){

	(*pos).a = (*pos).a + (*vel).a*t;
	(*pos).b = (*pos).b + (*vel).b*t;
	(*pos).c = (*pos).c + (*vel).c*t;

}

void randVeloc(point * vel){

	double theta,phi;
	double v0 = 10;	

	theta = pi * drand48();
	phi = 2.0 * pi * drand48();

	(*vel).a = v0 * sin(theta) * cos(phi);
	(*vel).b = v0 * sin(theta) * sin(phi);
	(*vel).c = v0 * cos(theta);

}

int main(){

	struct point {
		double a;
		double b;
		double c;
	} my_point;

	struct point *position = &my_point;
	struct point *velocity = &my_point;
	double timeStep, distance;
	int counter,steps;

	FILE *outfile;
	outfile = fopen("Ex_4.30.dat","w");

	randVeloc(velocity);

a
	counter=0;
	steps = 10000;
	distance = sqrt((*position).a*(*position).a + (*position).b*(*position).b + (*position).c*(*position).c);

	fprintf(outfile,"%i\t%lf\n",counter,distance);

	for(counter=0;counter<steps;counter++){

		step(position,velocity,timeStep);

		distance = sqrt((*position).a*(*position).a + (*position).b*(*position).b + (*position).c*(*position).c);

		fprintf(outfile,"%i\t%lf\n",counter,distance);

		randVeloc(velocity);

	}

}
