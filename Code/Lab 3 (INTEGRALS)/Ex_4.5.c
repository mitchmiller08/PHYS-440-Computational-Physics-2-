#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 3
// Exc. 4.4 2/28/12
// Romber Method

double func(double x){

	// Specify desired function
	double funcValue;
	funcValue = pow((1-x*x)*(2-x),0.5);
	return(funcValue);

}

double main(){

	double** t;
	int i,j,m,nMax;
	t = (double**) malloc(10*sizeof(double*));
	for(i=0;i<10;i++){
		t[i] = (double*) malloc(10*sizeof(double));
	}
	double a,b,sum,h;
	double rm;
	double a1,a2,a3,a4;

	a = -1.;
	b = 1.;

	for(i=0;i<10;i++){
		for(j=0;j<10;j++){
			t[i][j] = 0;
		}
	}

	t[0][0]=0.5*(b-a)*(func(a)+func(b));
	
	for(m=0;m<9;m++){

		nMax = pow(2.,m+1) - 1;
		sum = 0;
		h = (b-a)/(pow(2,m+1));
		
		for(i=1;i<=nMax;i+=2){

			sum = sum + func(a+i*h);

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

	for(i=0;i<10;i++){

		a1=t[i][0];
		a2=t[i][1];
		a3=t[i][2];
		a4=t[i][3];
		if(i>0&&i<9){

			rm = (t[i-1][0] - t[i][0]) / (t[i][0] - t[i+1][0]);

		}
		
		printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\n",i,a1,a2,a3,a4,rm);

	}

}
