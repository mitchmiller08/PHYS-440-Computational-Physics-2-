#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.2 1/31/12
// Lagrange Interpolation

// Interpolate the J0 Bessel function
// 21 Data points between 0 and 10
// http://keisan.casio.com/has10/SpecExec.cgi

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

	int numberPoints = 3;
	int counter1;
	double positionMatrix[3] = {2.0,2.5,3.0};
	double valueMatrix[3] = {0.2238907791,-0.04838377647,-0.2600519549};
	double pointValue;
	double calcPoint=1.0;
	FILE *outfile;
	outfile = fopen("Ex_3.3.j0.3point.dat","w");

	//Calculate interpolation at points between -1 and 1
	for(counter1=0;counter1<40;counter1++){
		
		pointValue = lagrangeInterpolation(numberPoints,positionMatrix,valueMatrix,calcPoint);
		fprintf(outfile,"%lf\t%lf\n",calcPoint,pointValue);
		calcPoint = calcPoint + 0.1;
		printf("CalcPoint = %lf\n",calcPoint);

	}
	fclose(outfile);

}

