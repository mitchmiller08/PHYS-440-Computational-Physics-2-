#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.2 2/11/12
// Hermite Interpolation

// I/I0 = (2*J0(x)/x)^2

// Interpolate the Airy pattern
// 10 Data points between 1 and 10
// http://keisan.casio.com/has10/SpecExec.cgi

double hermiteInterpolation(int numberPoints, double *positionMatrix, double *valueMatrix, double *derivativeMatrix, double interpolationPoint){

	int counter1,counter2;
	double product,sum=0,poly=0;
	double hj, hjp;

	for(counter1=0;counter1<numberPoints;counter1++){

		product = 1;
		sum = 0;
		for(counter2=0;counter2<numberPoints;counter2++){

			if(counter1!=counter2){

				//Calculate lambda
				product = product * ((interpolationPoint - positionMatrix[counter2]) / (positionMatrix[counter1] - positionMatrix[counter2]));
				sum = sum + 1./(positionMatrix[counter1] - positionMatrix[counter2]);

			}
			
		}

		//Calcualte g(x)
		hj = (1. - 2.*(interpolationPoint - positionMatrix[counter1]) * sum) * product * product;
		hjp = (interpolationPoint - positionMatrix[counter1]) * product * product;
		poly = poly + hj*valueMatrix[counter1] + hjp*derivativeMatrix[counter1];

	}

	printf("g(x) = %lf\n",poly);
	return poly;

}

int main(){

	int numberPoints = 11;
	int counter1;
	double positionMatrix[11] = {0.1,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
	double valueMatrix[11] = {0.997503, 0.774578, 0.332612, 0.0510938, 0.00109043, 0.0171693, 0.008506, 0.00000179011, 0.00344089, 0.00297175, 0.0000755952};
	double derivativeMatrix[11] = {-0.0498959,-0.404507,-0.406976,-0.146501,0.0120241,0.0048812,-0.0149331,-0.000230446,0.003314,-0.00350941,-0.000885558};
	double pointValue;
	double calcPoint=0.0;
	FILE *outfile;
	outfile = fopen("Ex_3.4.airy.dat","w");

	//Calculate interpolation at points between -1 and 1
	for(counter1=0;counter1<100;counter1++){
		
		printf("CalcPoint = %lf\n",calcPoint);
		pointValue = hermiteInterpolation(numberPoints,positionMatrix,valueMatrix,derivativeMatrix,calcPoint);
		fprintf(outfile,"%lf\t%lf\n",calcPoint,pointValue);
		calcPoint = calcPoint + 0.1;

	}
	fclose(outfile);

}

