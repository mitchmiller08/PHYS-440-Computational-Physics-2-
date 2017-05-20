#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 4
// Fig. 12.12 4/17/12
// Numerical Differential Equations

double *ders(double x,double *y,double *f){

	f[0] = y[1];
	//f[1] = -1.0*y[0] - 3*(y[0]*y[0]-1) * y[1];
	f[1] = -0.24 * y[1] + y[0] - y[0]*y[0]*y[0] + 0.68*cos(1.7*x);
	//f[1] = -0.1*y[1] - sin(y[0]) + 0.5*sin(2/3*x);

	return(f);

}

int main(){

	// Initialize variables
	double v0,h,hnew,k,m,g,hmax,hmin;
	double f0[2],f1[2],f2[2],f3[2],f4[2],f5[2];
	double y0[2],y[2],yhat[2],x,x0;
	double err,bigerr,maxerr,epsilon;
	double a1,a2,a3,a4,a5;
	double b10,b11,b20,b21;
	double b30,b31,b32,b40,b41,b42,b43;
	double b50,b51,b52,b53,b54;
	double c1,c2,c3,c4,c5,d1,d2,d3,d4,d5;
	int ncount, counter, counter2;
	FILE *outfile;
	outfile = fopen("Fig_12.12.dat","w");

	// Define constants
	a1=1./4;
	a2=3./8;
	a3=12./13;
	a4=1.;
	a5=1./2;

	b10=0.25;
	b20=3./32;
	b21=9./32;
	b30=1932./2197;
	b31=7200./2197;
	b32=7296./2197;
	b40=439./216;
	b41=8.;
	b42=3680./513;
	b43=845./4104;
	b50=8./27;
	b51=2.;
	b52=3544./2565;
	b53=1859./4104;
	b54=1.1/40;
	
	c1=16./135;
	c2=6656/12825;
	c3=28561./56430;
	c4=9./50;
	c5=2./55;

	d1=1./360;
	d2=128./4275;
	d3=2197./75240;
	d4=1./50;
	d5=2./55;

	x0=0.0;
	y0[0]=1.0;
	y0[1]=-1.0;
	hmax=0.00001;
	hmin=0.000001;
	
	h= (2.*pi/1.7)/150;
	epsilon=3e-5;
	ncount=0;
	
	for(counter=0;counter<6000000;counter++){

		ncount++;
		
		ders(x0,y0,f0);
		x = x0 + a1*h;
		for(counter2=0;counter2<2;counter2++){

			y[counter2] = y0[counter2]+b10*h*f0[counter2];
	
		}
		ders(x,y,f1);

		x = x0 + a2*h;
		for(counter2=0;counter2<2;counter2++){

			y[counter2] = y0[counter2] + b20*h*f0[counter2] + b21*h*f1[counter2];
	
		}
		ders(x,y,f2);

		x0 + a3*h;
		for(counter2=0;counter2<2;counter2++){

			y[counter2] = y0[counter2] + b30*h*f0[counter2] - b31*h*f1[counter2] + b32*h*f2[counter2];
	
		}
		ders(x,y,f3);

		x = x0 + a4*h;
		for(counter2=0;counter2<2;counter2++){

			y[counter2] = y0[counter2] + b40*h*f0[counter2] - b41*h*f1[counter2] + b42*h*f2[counter2] - b43*h*f3[counter2];
	
		}
		ders(x,y,f4);

		x = x0 + a5*h;
		for(counter2=0;counter2<2;counter2++){

			y[counter2] = y0[counter2] - b50*h*f0[counter2] + b51*h*f1[counter2] - b52*h*f2[counter2] + b53*h*f3[counter2] - b53*h*f4[counter2];
	
		}
		ders(x,y,f5);

		bigerr = 0.;
		
		for(counter2=0;counter2<2;counter2++){

			yhat[counter2] = y0[counter2] + h*(c1*f1[counter2] + c2*f2[counter2] + c3*f3[counter2] - c4*f4[counter2] + c5*f5[counter2]);
			err = h*fabs(d1*f0[counter2] - d2*f2[counter2] - d3*f3[counter2] + d4*f4[counter2] + d5*f5[counter2]);
			if(err>bigerr){
				bigerr = err;
			}
	
		}

		maxerr = h*epsilon;
		hnew = 0.9*h*sqrt(sqrt(maxerr/err));
		hnew = h;

		if(hnew>4*h){
			hnew = 4*h;
		}
		if(hnew<0.1*h){
			hnew = 0.1*h;
		}
		if(hnew>hmax){
			hnew = hmax;
		}
		if(hnew<hmin){
			exit (1);
		}

		x0 = x0+h;
		for(counter2=0;counter2<2;counter2++){
			y0[counter2] = yhat[counter2];
		}
		/*
		if(y[2]>=pi){
			y[2] = y[2] - 2.*pi;
		}
		if(y[2]<=-pi){
			y[2] = y[2] + 2.*pi;
		}*/
		
		//if(counter>10000){
		//ncount++;
		if(ncount==150){
			fprintf(outfile,"%.15lf\t%.15lf\n",y[0],y[1]);
			ncount = 0;
		//}
		}
		
	}

	fclose(outfile);
		
}






