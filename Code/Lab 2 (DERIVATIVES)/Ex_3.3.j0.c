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

	int numberPoints = 21;
	int counter1;
	double positionMatrix[21] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0};
	double valueMatrix[21] = {1,0.9384698072,0.7651976866,0.5118276717,0.2238907791,-0.04838377647,-0.2600519549,-0.38012774,-0.3971498099,-0.320542509,-0.1775967713,-0.0068438694,0.1506452573,0.260094606,0.3000792705,0.2663396579,0.1716508071,0.0419392518,-0.0903336112,-0.1939287477,-0.2459357645};
	double pointValue;
	double calcPoint=0.0;
	FILE *outfile;
	outfile = fopen("Ex_3.3.j0.dat","w");

	//Calculate interpolation at points between -1 and 1
	for(counter1=0;counter1<100;counter1++){
		
		pointValue = lagrangeInterpolation(numberPoints,positionMatrix,valueMatrix,calcPoint);
		fprintf(outfile,"%lf\t%lf\n",calcPoint,pointValue);
		calcPoint = calcPoint + 0.1;
		printf("CalcPoint = %lf\n",calcPoint);

	}
	fclose(outfile);

}

