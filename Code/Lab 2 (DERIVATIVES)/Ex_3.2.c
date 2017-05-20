#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.2 1/31/12
// Lagrange Interpolation

double lagrangeInterpolation(int numberPoints, double *positionMatrix, double *valueMatrix, double interpolationPoint){

	int counter1,counter2;
	double product,sum=0;

	for(counter1=0;counter1<numberPoints;counter1++){

		product = 1;
		for(counter2=0;counter2<numberPoints;counter2++){

			if(counter1!=counter2){

				//Calculate lambda
				product = product * ((interpolationPoint - positionMatrix[counter2]) / (positionMatrix[counter1] - positionMatrix[counter2]));
				//printf("prod = %lf\n",product);

			}
			
		}

		//Calcualte g(x)
		sum += product*valueMatrix[counter1];

	}

	printf("g(x) = %lf\n",sum);
	return sum;

}

int main(){

	int numberPoints = 21;
	int counter1;
	double positionMatrix[21] = {-5.,-4.5,-4.,-3.5,-3.,-2.5,-2.,-1.5,-1.,-0.5,0,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.};
	double valueMatrix[21] = {0.,0.,0.,0.,0.,0.,0.,0.00000000000000935918,000000000.206115,0.0000453979,0.5,0.999955,1.,1.,1.,1.,1.,1.,1.,1.};
	double pointValue;
	double calcPoint=-1.0;
	FILE *outfile;
	outfile = fopen("Ex_3.2.out.dat","w");

	//Calculate interpolation at points between -1 and 1
	for(counter1=0;counter1<200;counter1++){
		
		pointValue = lagrangeInterpolation(numberPoints,positionMatrix,valueMatrix,calcPoint);
		fprintf(outfile,"%lf\t%lf\n",calcPoint,pointValue);
		calcPoint = calcPoint + 0.01;
		printf("CalcPoint = %lf\n",calcPoint);

	}
	fclose(outfile);

}

