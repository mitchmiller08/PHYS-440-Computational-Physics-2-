#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 4
// Ex. 7.8 4/17/12
// Numerical Differential Equations

int main(){

	double **u,hx=0.1,hy=0.1,uu,alpha=0.49;
	int i,j,count;
	int xSteps=9, ySteps=9,flag=0;
	FILE *outfile;
	outfile = fopen("Ex_7.8.dat","w");

	u = (double**)malloc(xSteps*sizeof(double*));
	for(i=0;i<xSteps;i++){
		u[i] = (double*)malloc(ySteps*sizeof(double));
	}
	
	for(i=0;i<xSteps;i++){
		for(j=0;j<ySteps;j++){

			u[i][j] = 0.01;

		}
	}

	for(j=0;j<ySteps;j++){

		u[0][j] = 10;
		u[8][j] = 0;

	}

	for(i=0;i<xSteps;i++){

		u[i][0] = 0;
		u[i][8] = 10;

	}
	
	count=0;

	while(flag==0){

		count++;
		if(count>100){
			printf("Too many iterations");
			break;
		}

		for(i=1;i<xSteps-1;i++){

			for(j=1;j<ySteps-1;j++){

				uu = ((u[i+1][j] + u[i-1][j])/(hx*hx) + (u[i][j+1] + u[i][j-1])/(hy*hy)) * pow(hx*hy,2);
				uu = uu / (2.*(hx*hx + hy*hy));

				if(fabs((uu-u[i][j])/uu)< 5e-5){
					flag=1;
				}

				u[i][j] = uu + alpha*(uu - u[i][j]);

			}

		}

		printf("Counts = %d\n",count);

	}

	for(i=0;i<xSteps;i++){

		for(j=0;j<ySteps;j++){

			fprintf(outfile,"%d\t%d\t%.15lf\n",i,j,u[i][j]);

		}

	}

}
