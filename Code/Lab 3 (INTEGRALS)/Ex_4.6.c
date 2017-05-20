#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

// Mitchell Miller
// PHYS 440 Lab 3
// Exc. 4.6 2/28/12
// Romberg Method

double func(double theta, double x){

	// Period integration
	double funcValue;
	double length=1;
	double gravity=9.81;
	funcValue = 4.*pow(length/(gravity*(1-sin(theta/2)*sin(theta/2)*sin(x)*sin(x))),0.5);
	return(funcValue);

}

double romberg(double a, double b){

	double** t;
	double c = b;
	b = pi/2.0;
	int i,j,m,nMax;
	t = (double**) malloc(10*sizeof(double*));
	for(i=0;i<10;i++){
		t[i] = (double*) malloc(10*sizeof(double));
	}
	double sum,h;
	double a1,a2,a3,a4;
	double final;

	for(i=0;i<10;i++){
		for(j=0;j<10;j++){
			t[i][j] = 0;
		}
	}

	t[0][0]=0.5*(b-a)*(func(c,a)+func(c,b));
	
	for(m=0;m<9;m++){

		nMax = pow(2.,m+1) - 1;
		sum = 0;
		h = (b-a)/(pow(2,m+1));
		
		for(i=1;i<=nMax;i+=2){

			sum = sum + func(c,a+i*h);

		}

		t[m+1][0] = 0.5 * t[m][0] + h*sum;

	}

	for(i=1;i<10;i++){

		for(j=1;j<=i;j++){

			if(j<=3){

				t[i][j] = (pow(4,j) * t[i][j-1] - t[i-1][j-1]) /  (pow(4,j)-1);

			}

		}

	}
/*
	for(i=0;i<10;i++){

		a1=t[i][0];
		a2=t[i][1];
		a3=t[i][2];
		a4=t[i][3];

		printf("%i\t%lf\t%lf\t%lf\t%lf\n",i,a1,a2,a3,a4);

	}*/
	
	final = t[9][3];
	return(final);

}

int main(){

	FILE *outfile;
	outfile = fopen("Ex_4.6.dat","w");

	double theta0;
	double T,degrees;
	double Tactual = 2.00606668071065;
	double error;
	int counter;

	for(counter=0;counter<314;counter++){

		theta0 = counter/100.;
		degrees = theta0 * (180/pi);
		T = romberg(0,theta0);
		error = fabs((T-Tactual)/Tactual);
		fprintf(outfile,"%lf\t%lf\t%lf\n",degrees,T,error);

	}

}
