#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.2 1/31/12
// Lagrange Interpolation

// Interpolate the J1 Bessel function
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

	int numberPoints = 21;
	int counter1;
	double positionMatrix[21] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0};
	double valueMatrix[21] = {0,0.2422684577,0.4400505857,0.5579365079,0.5767248078,0.497094102,0.3390589585,0.1373775274,-0.06604332802,-0.2310604319,-0.3275791376,-0.3414382154,-0.2766838581,-0.1538413014,-0.0046828235,0.1352484276,0.2346363469,0.2731219637,0.2453117866,0.1612644308,0.04347274617};
	double pointValue;
	double calcPoint=0.0;
	FILE *outfile;
	outfile = fopen("Ex_3.3.j1.dat","w");

	//Calculate interpolation at points between -1 and 1
	for(counter1=0;counter1<100;counter1++){
		
		pointValue = lagrangeInterpolation(numberPoints,positionMatrix,valueMatrix,calcPoint);
		fprintf(outfile,"%lf\t%lf\n",calcPoint,pointValue);
		calcPoint = calcPoint + 0.1;
		printf("CalcPoint = %lf\n",calcPoint);

	}
	fclose(outfile);

}

