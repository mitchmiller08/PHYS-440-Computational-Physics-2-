#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.2 2/11/12
// Hermite Interpolation

// J1' = (J0 - J2) / 2

// Interpolate the J1 Bessel function
// 11 Data points between 0 and 10
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
	double positionMatrix[11] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0};
	double valueMatrix[11] = {0.,0.4400505857,0.5767248078,0.3390589585,-0.06604332802,-0.3275791376,-0.2766838581,-0.0046828235,0.2346363469,0.2453117866,0.04347274617};
	double derivativeMatrix[11] = {0.5, 0.325147, -0.0644716, -0.373072, -0.380639, -0.112081, 0.196759, 0.300748, 0.142321, -0.11759, -0.250283};
	double pointValue;
	double calcPoint=0.0;
	FILE *outfile;
	outfile = fopen("Ex_3.4.j1.dat","w");

	//Calculate interpolation at points between -1 and 1
	for(counter1=0;counter1<100;counter1++){
		
		printf("CalcPoint = %lf\n",calcPoint);
		pointValue = hermiteInterpolation(numberPoints,positionMatrix,valueMatrix,derivativeMatrix,calcPoint);
		fprintf(outfile,"%lf\t%lf\n",calcPoint,pointValue);
		calcPoint = calcPoint + 0.1;

	}
	fclose(outfile);

}

