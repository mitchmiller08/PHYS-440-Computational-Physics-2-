#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654

// Mitchell Miller
// PHYS 440 Lab 3
// Exc. 4.4 2/28/12
// Romberg Method

double func1(double x){

	// C(v) Fresnel Integral
	double funcValue;
	funcValue = cos(pi*x*x/2.);
	return(funcValue);

}

double func2(double x){

	// S(v) Fresnel Integral
	double funcValue;
	funcValue = sin(pi*x*x/2.);
	return(funcValue);

}

double romberg(double a, double b, int flag){

	double** t;
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

	if(flag==0){
		t[0][0]=0.5*(b-a)*(func1(a)+func1(b));
	}
	else if(flag==1){
		t[0][0]=0.5*(b-a)*(func2(a)+func2(b));
	}
	
	for(m=0;m<9;m++){

		nMax = pow(2.,m+1) - 1;
		sum = 0;
		h = (b-a)/(pow(2,m+1));
		
		for(i=1;i<=nMax;i+=2){

			if(flag==0){

				sum = sum + func1(a+i*h);
		
			}
			else if(flag==1){

				sum = sum + func2(a+i*h);

			}

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
	outfile = fopen("Ex_4.4.dat","w");

	double v,I;
	double C,S;
	int counter;

	v=0;

	for(counter=0;counter<1000;counter++){

		v = counter/100.;
		C = romberg(0,v,0);
		S = romberg(0,v,1);
		I = 0.5 * ( (C+0.5)*(C+0.5) + (S+0.5)*(S+0.5) );
		fprintf(outfile,"%lf\t%lf\n",v,I);

	}

}
