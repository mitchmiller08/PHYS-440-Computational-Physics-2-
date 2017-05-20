#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.20 2/18/12
// Least Squares Fitting

void *luSolve(double **a, double *x, double *b, int n){

	int order[15],i,j,k,imax,itemp;
	double scale[15],max,sum,temp,det;

	det = 1.0;

	for(i=0;i<n;i++){

		order[i] = i;
		max = 0.;

		for(j=0;j<n;j++){

			if(fabs(a[i][j])>max){

				max = fabs(a[i][j]);

			}

		}

		scale[i] = 1.0 / max;

	}

	for(k=0;k<n-1;k++){

		if(k!=0){

			for(i=k;i<n;i++){

				sum = a[i][k];

				for(j=0;j<=k-1;j++){

					sum = sum - a[i][j]*a[j][k];

				}

				a[i][k] = sum;

			}

		}

		max = 0;

		for(i=k;i<n;i++){

			if(scale[i]*a[i][k] >= max){

				max = scale[i]*a[i][k];
				imax = i;

			}

		}

		if(imax!=k){

			det = -det;

			for(j=0;j<n;j++){

				temp = a[imax][j];
				a[imax][j] = a[k][j];
				a[k][j] = temp;

			}

			temp = scale[imax];
			scale[imax] = scale[k];
			scale[k] = temp;

			itemp = order[imax];
			order[imax] = order[k];
			order[k] = itemp;

		}

		det = det*a[k][k];

		if(k==0){

			for(j=1;j<n;j++){

				a[0][j] = a[0][j] / a[0][0];

			}

		}
		else{

			for(j=k+1;j<n;j++){

				sum = a[k][j];
		
				for(i=0;i<=k-1;i++){

					sum = sum - a[k][i]*a[i][j];
			
				}

				a[k][j] = sum / a[k][k];

			}

		}

	}

	sum = a[n-1][n-1];
	for(j=0;j<n-1;j++){

		sum = sum - a[n-1][j]*a[j][n-1];

	}
	a[n-1][n-1] = sum;
	det = det*a[n-1][n-1];


	// Begin solution


	for(i=0;i<n;i++){

		x[i] = b[order[i]];

	}

	x[0] = x[0] / a[0][0];

	for(i=1;i<n;i++){

		sum = x[i];

		for(k=0;k<=i-1;k++){

			sum = sum - a[i][k]*x[k];

		}

		x[i] = sum / a[i][i];

	}

	for(i=1;i<n;i++){

		sum = x[n-1-i];

		for(k=n-1-i+1;k<n;k++){

			sum = sum - a[n-1-i][k] * x[k];

		}

		x[n-1-i] = sum;

	}

}

int main(){

	int N=11;
	double** a;
	int i;
	a = (double**) malloc(3*sizeof(double*));
	for(i=0;i<3;i++){
		a[i] = (double*) malloc(3*sizeof(double));
	}
	double *b = malloc(3*sizeof(double));
	double *x = malloc(3*sizeof(double));

	double positionMatrix[11] = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
	double valueMatrix[11] = {1.67203,1.79792,2.37791,2.66408,2.11245,2.43969,1.88843,1.59447,1.79634,1.07810,0.21066};
	FILE *outfile;
	outfile = fopen("Ex_3.20.dat","w");
	double temp=0;

	a[0][0]=11.;
	a[1][0]=5.5;
	a[0][1]=5.5;
	a[0][2]=3.85;
	a[2][0]=3.85;
	a[1][1]=3.85;
	a[1][2]=3.025;
	a[2][1]=3.025;
	a[2][2]=2.5333;

	b[0]=19.6321;
	b[1]=8.38663;
	b[2]=4.99548;

	/*
	// Fill in A matrix
	a[0][0]=11.;
	for(i=0;i<N;i++){

		temp = temp + positionMatrix[i];

	}
	a[0][1]=temp;
	a[1][0]=temp;
	temp=0;
	for(i=0;i<N;i++){

		temp = temp + positionMatrix[i]*positionMatrix[i];

	}
	a[0][2]=temp;
	a[1][1]=temp;
	a[2][0]=temp;
	temp=0;
	for(i=0;i<N;i++){

		temp = temp + positionMatrix[i]*positionMatrix[i]*positionMatrix[i];

	}
	a[1][2]=temp;
	a[2][1]=temp;
	temp=0;
	for(i=0;i<N;i++){

		temp = temp + positionMatrix[i]*positionMatrix[i]*positionMatrix[i]*positionMatrix[i];

	}
	a[2][2]=temp;
	temp=0;

	// Fill in B matrix
	for(i=0;i<N;i++){

		temp = temp + valueMatrix[i];

	}
	b[0] = temp;
	temp=0;
	for(i=0;i<N;i++){

		temp = temp + positionMatrix[i]*valueMatrix[i];

	}
	b[1]=temp;
	temp=0;
	for(i=0;i<N;i++){

		temp = temp + positionMatrix[i]*positionMatrix[i]*valueMatrix[i];

	}
	b[2]=temp;*/

	// Solve normal equations
	luSolve(a,x,b,3);

	for(i=0;i<3;i++){

		printf("C_%i = %lf\n",i,x[i]);

	}

	for(i=0;i<11;i++){

		fprintf(outfile,"%lf\t%lf\n",positionMatrix[i],valueMatrix[i]);

	}

}
