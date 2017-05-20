#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.21 2/18/12
// Non-Linear Least Squares Fitting

void *init(double *a, double *h){

	int max=6;
	int m=4;

	if(m>max){

		printf("Dimensional error in init\n");
		exit(1);

	}

	a[0] = 435.84;
	a[1] = 0.03;
	a[2] = 120.;
	a[3] = 2.0;

	h[0] = 0.01;
	h[1] = 0.005;
	h[2] = 4.0;
	h[3] = 2.0;

}
/*
double sum(double *array, int size){

	int i;
	double sumA=0;
	for(i=0;i<size;i++){

		sumA = sumA + array[i];

	}
	return(sumA);

}*/

double sum(double *a){

	double x[21],y[21];
	double tc,lambda,lambdaO,baseline;
	double intensity,gamma;
	double sumA;
	int i;

	x[0] = 435.784;
	x[1] = 435.789;
	x[2] = 435.794;
	x[3] = 435.799;
	x[4] = 435.804;
	x[5] = 435.809;
	x[6] = 435.814;
	x[7] = 435.819;
	x[8] = 435.824;
	x[9] = 435.829;
	x[10] = 435.834;
	x[11] = 435.839;
	x[12] = 435.844;
	x[13] = 435.849;
	x[14] = 435.854;
	x[15] = 435.859;
	x[16] = 435.864;
	x[17] = 435.869;
	x[18] = 435.874;
	x[19] = 435.879;
	x[20] = 435.884;
	
	y[0] = 40.0;
	y[1] = 44.0;
	y[2] = 41.0;
	y[3] = 46.0;
	y[4] = 47.5;
	y[5] = 54.0;
	y[6] = 69.7;
	y[7] = 97.0;
	y[8] = 129.0;
	y[9] = 153.0;
	y[10] = 165.0;
	y[11] = 168.0;
	y[12] = 143.0;
	y[13] = 111.0;
	y[14] = 79.0;
	y[15] = 64.0;
	y[16] = 52.0;
	y[17] = 51.0;
	y[18] = 44.0;
	y[19] = 46.0;
	y[20] = 41.0;

	lambdaO = a[0];
	gamma = a[1];
	intensity = a[2];
	baseline = a[3];

	sumA = 0;
	for(i=0;i<21;i++){

		lambda = x[i];
		tc = baseline;
		tc = tc + intensity / (1.+4.*(lambda-lambdaO)*(lambda-lambdaO) / (gamma*gamma));

		sumA = sumA + (y[i] - tc) * (y[i] - tc);

	}

	return(sumA);

}

void *crude(double *a, double *h){

	int i,k,m=4;
	double sp,sO,sm,test;
	double aPlus[6],aMinus[6];

	for(i=0;i<m;i++){

		for(k=0;k<m;k++){

			if(k==i){

				aPlus[i] = a[i]+h[i];
				aMinus[i] = a[i]-h[i];

			}
			else{

				aPlus[k] = a[k];
				aMinus[k] = a[k];

			}

		}

		sp = sum(aPlus);
		sO = sum(a);
		sm = sum(aMinus);

		test = sp-2.0*sO+sm;

		if(test==0.0){

			test = 0.0000001;

		}

		printf("Test = %lf\n",test);

		a[i] = a[i] - 0.5*h[i] * (sp-sm) / test;

		h[i] = 0.85*h[i];

	}

	printf("S_0 = %lf\n", sO);

}

int main(){

	double *a = malloc(6*sizeof(double));
	double *h = malloc(6*sizeof(double));
	int ik;

	init(a,h);

	for(ik=0;ik<100;ik++){

		crude(a,h);
		printf("Loop %i\n",ik);
		printf("a[1] = %lf\ta[2] = %lf\n",a[0],a[1]);
		printf("a[3] = %lf\ta[4] = %lf\n",a[2],a[3]);

	}

}
