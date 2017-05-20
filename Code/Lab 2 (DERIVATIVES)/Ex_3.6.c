#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.6 2/11/12
// Cubic Spline Fitting

// NOTE THAT N SHOULD BE THE INDEX VALUE, NOT THE NUMBER OF POINTS

double spline(double xValue, double *x, double *f, double *second, int N){

	int counter = 0;
	double sp1,sp2,sp3,sp4,spline;	
	
	// Check bounds
	if(xValue < x[0]){

		printf("X Value is less than x[0] in spline");
		exit(1);
	
	}
	if(xValue > x[N]){

		printf("X Value is greater than x[0] in spline");
		exit(1);

	}


	// Find interval containing xValue
	while(xValue > x[counter+1]){

		counter = counter + 1;

	}

	// Begin interpolation
	sp1 = (xValue - x[counter]) * ((f[counter+1] - f[counter]) / (x[counter+1] - x[counter]));
	sp2 = (xValue - x[counter]) * (x[counter+1] - x[counter]) * (second[counter+1]/6.0 + second[counter]/3.0);
	sp3 = 0.5 * second[counter] * (xValue - x[counter]) * (xValue - x[counter]);
	sp4 = ((second[counter+1] - second[counter]) * (xValue - x[counter]) * (xValue - x[counter]) * (xValue - x[counter]) / (6.0 * (x[counter+1] - x[counter])));
	spline = f[counter] + sp1 - sp2 + sp3 + sp4;

	return spline;

}

void *splineInit(double *x, double *f, double *second, double fp1, double fpn, int N){

	int counter = 0;
	double beta[11] = { 0 };
	double b0,r0,bj,rj,aj,c;
	
	b0 = 2.0 * (x[1] - x[0]);
	beta[0] = b0;
	r0 = 6.0 * ((f[1] - f[0]) / (x[1] - x[0]) - fp1);

	for(counter=1;counter<=N;counter++){

		if(counter == N){

			bj = 2.0 * (x[N] - x[N-1]);
			rj = -6.0 * ((f[N] - f[N-1]) / (x[N] - x[N-1]) - fpn);

		}
		else{
			
			bj = 2.0 * (x[counter+1] - x[counter-1]);
			rj = 6.0 * ((f[counter+1] - f[counter]) / (x[counter+1] - x[counter]) - (f[counter] - f[counter-1]) / (x[counter] - x[counter-1]));

		}

		aj = x[counter] - x[counter-1];
		c = aj;
		beta[counter] = bj - aj * c / beta[counter-1];
		second[counter] = rj - aj * second[counter-1] / beta[counter-1];

		if(beta[counter] == 0.){
	
			printf("Zero element in diagonal of splineinit");
			exit(1);

		}

	}

	second[N] = second[N]/beta[N];

	for(counter=0;counter<N;counter++){

		c = x[N-counter+1] - x[N-counter];
		second[N-counter] = (second[N-counter] - c * second[N-counter+1]) / beta[N-counter];

	}

}

int main(){

	int numberOfPoints = 10;
	int counter = 0;
	double *secondDeriv = malloc(11*sizeof(double));
	double positionMatrix[11] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
	double valueMatrix[11] = {0.,0.4400505857,0.5767248078,0.3390589585,-0.06604332802,-0.3275791376,-0.2766838581,-0.0046828235,0.2346363469,0.2453117866,0.04347274617};
	double valueDeriv1 = 0.5;
	double valueDerivN = -0.250283;
	double interpPoint = 0.0;
	double interpValue;
	FILE *outfile;
	outfile = fopen("Ex_3.6.dat","w");
	
	splineInit(positionMatrix,valueMatrix,secondDeriv,valueDeriv1,valueDerivN,numberOfPoints);

	/*for(counter=0;counter<11;counter++){
		printf("g''(%lf) = %ls\n",counter,secondDeriv[counter]);
	}*/

	for(counter=0;counter<100;counter++){

		interpValue = spline(interpPoint,positionMatrix,valueMatrix,secondDeriv,numberOfPoints);
		printf("g(%lf) = %lf\n", interpPoint, interpValue);
		fprintf(outfile,"%lf\t%lf\n",interpPoint,interpValue);
		interpPoint = interpPoint + 0.1;

	}

}
