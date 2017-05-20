#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
//#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 4
// Ex. 5.33 4/17/12
// Numerical Differential Equations

int main(){

	double h,c1,c2,c3,c4,x,yy,tol;
	int counter, loops;
	int done=0;
	double *y;
	y = (double*)malloc(11*sizeof(double));
	FILE *outfile;
	outfile = fopen("Ex_5.33.dat","w");

	tol = 5e-4;
	h=0.1;
	c1=1-2.5*h;
	c2=1+2.5*h;
	c3=-10*h*h;
	c4=2+c3;

	for(counter=0;counter<11;counter++){

		y[counter]=100*(counter/10);

	}

	loops=0;
	
	while(done==0){

		loops++;
		if(loops>=100){
			break;
		}

		for(counter=1;counter<10;counter++){

			x = 0.1 + 0.1*(counter-1);
			yy = (c1*y[counter+1]+c2*y[counter-1]+c3*x)/c4;

			if(fabs(yy-y[counter])/yy < tol){
				done=1;
			}
			else{
				done=0;
			}

			y[counter]=yy;

		}

	}

	printf("Loops = %i\n",loops);

	for(counter=0;counter<11;counter++){

		x=0+0.1*counter;
		fprintf(outfile, "%.15lf\t%.15lf\n",x,y[counter]);

	}

}
