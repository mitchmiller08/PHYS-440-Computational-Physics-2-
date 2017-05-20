//  Created by Mitchell Miller on 4/11/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
//#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 4
// Fig. 5.10 4/17/12
// Numerical Differential Equations

int main(){

	double x0,x1,x2,x3;
	double v0,v1,v2,v3;
	double f0,f1,f2,f3;
	double x=3.,v=0,t=0,e=1;
	double h = 2*pi/150;
	int counter,ncount=0;
	FILE *outfile;
	outfile = fopen("Ex_5.19.3.dat","w");
	fprintf(outfile,"%.15lf\t%.15lf\n",x,v);

	for(counter=0;counter<1000;counter++){

		ncount++;
		
		f0 = -1.0*x-e*(x*x-1)*v;

		x1 = x+0.5*h*v;
		v1 = v+0.5*h*f0;
		f1 = -1.*x1-e*(x1*x1-1)*v1;

		x2 = x+0.5*h*v + 0.25*h*f0;
		v2 = v+0.5*h*f1;
		f2 = -1.*x2-e*(x2*x2-1)*v2;

		x3 = x+h*v+0.5*h*f1;
		v3 = v+h*f2;
		f3 = -1.*x3-e*(x3*x3-1)*v3;

		x = x+h*(v+h*(f0+f1+f2)/6);
		v = v+h*(f0+2*f1+2*f2+f3)/6;
		t = t+h;

		//if(ncount==10){

			fprintf(outfile,"%.15lf\t%.15lf\n",x,v);
			ncount = 0;

		//}

	}

	fclose(outfile);

}
