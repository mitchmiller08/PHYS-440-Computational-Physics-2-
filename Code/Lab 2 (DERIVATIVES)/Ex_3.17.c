#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.17 2/18/12
// LU Solve Method

// n is the number of index points, not number of points

double *luSolve(double **a, double *x, double *b, double det, int ndim, int n){

	int order[10],i,j,k,imax,itemp;
	double scale[10],max,sum,temp;

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

			for(j=1;j<n;j++){

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

	double** a;
	int i;
	a = (double**) malloc(5*sizeof(double*));
	for(i=0;i<5;i++){
		a[i] = (double*) malloc(5*sizeof(double));
	}
	double *b = malloc(5*sizeof(double));
	double *x = malloc(5*sizeof(double));
	double det;

	a[0][0] = 2.0;
	a[4][4] = 2.0;
	
	for(i=1;i<5;i++){

		a[i][i] = 2.0;
		a[i-1][i] = -1.0;
		a[i][i-1] = -1.0;

	}

	b[0] = 0.0;
	b[1] = 1.0;
	b[2] = 2.0;
	b[3] = 3.0;
	b[4] = 4.0;

	luSolve(a,x,b,det,5,5.0);

	for(i=0;i<5;i++){

		printf("X[%i] = %lf\n",i,x[i]);

	}

}
