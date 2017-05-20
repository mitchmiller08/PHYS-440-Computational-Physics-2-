#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.2 1/31/12
// Lagrange Interpolation

// I/I0 = (2*J0(x)/x)^2

// Interpolate the Airy pattern
// 10 Data points between 1 and 10
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

	int numberPoints = 11;
	int counter1;
	double positionMatrix[11] = {0.1, 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
	double valueMatrix[11] = {0.997503, 0.774578, 0.332612, 0.0510938, 0.00109043, 0.0171693, 0.008506, 0.00000179011, 0.00344089, 0.00297175, 0.0000755952};
	double pointValue;
	double calcPoint=0.0;
	FILE *outfile;
	FILE *outfileI;
	outfile = fopen("Ex_3.3.airy.dat","w");
	outfileI = fopen("airy.dat","w");

	//Calculate interpolation at points between -1 and 1
	for(counter1=0;counter1<100;counter1++){
		
		pointValue = lagrangeInterpolation(numberPoints,positionMatrix,valueMatrix,calcPoint);
		fprintf(outfile,"%lf\t%lf\n",calcPoint,pointValue);
		calcPoint = calcPoint + 0.1;
		printf("CalcPoint = %lf\n",calcPoint);

	}
	for(counter1=0;counter1<11;counter1++){
		
		fprintf(outfileI,"%lf\t%lf\n",positionMatrix[counter1],valueMatrix[counter1]);

	}
	fclose(outfile);

}

