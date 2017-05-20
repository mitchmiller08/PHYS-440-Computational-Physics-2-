#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 4
// Ex. 7.3 4/17/12
// Numerical Differential Equations

int main(){

	double *uold,*u,*unew,*uo,*duo;
	double time,dt,dx,epsilon,T,mu,c;
	int i,nx,nt,n;
	double h,x,maxTime;
	
	uold = (double*)malloc(100*sizeof(double));
	u = (double*)malloc(100*sizeof(double));
	unew = (double*)malloc(100*sizeof(double));
	uo = (double*)malloc(100*sizeof(double));
	duo = (double*)malloc(100*sizeof(double));

	FILE *outfile;
	outfile = fopen("Ex_7.3.dat","w");

	n=100;	
	h=1./100;
	mu=1e-3;
	dt=1e-4;
	T=10;
	dx=1e-2;
	time=0;
	nx=(int)(1/dx);
	nt=60;
	maxTime=dt*nt;
	epsilon = pow((dt * sqrt(T/mu) / dx),2);
	c = sqrt(T/mu);

	for(i=0;i<n;i++){

		x = h*(i+1);
		uo[i] = pow(e,(-100*(x-0.5)*(x-0.5)));
		duo[i] = -200*c*(x-0.5)*pow(e,-100*pow(x+c*time-0.5,2));

	}
	uo[0] = 0;
	uo[99]= 0;
	duo[0] = 0;
	duo[99] = 0;

	for(i=0;i<n;i++){

		uold[i]=uo[i];

	}

	u[0] = 0;
	u[n-1] = 0;
	
	for(i=1;i<(n-1);i++){

		u[i] = 0.5*epsilon*(uo[i+1]+uo[i-1]) + (1-epsilon)*uo[i] + dt*duo[i];

	}

	while(time<maxTime){

		time = time+dt;
		unew[0] = 0;
		unew[n-1] = 0;
	
		for(i=1;i<(n-1);i++){

			unew[i] = epsilon*(u[i+1]+u[i-1]) + 2*(1-epsilon)*u[i] - uold[i];

		}

		for(i=0;i<n;i++){

			uold[i] = u[i];
			u[i] = unew[i];
		
		}

	}


	for(i=0;i<n;i++){

		x = h*(i+1);
		fprintf(outfile,"%.15lf\t%.15lf\n",x,u[i]);

	}

}
