#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define e 2.71828182845

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.13 2/17/12
// Richardson Extrapolation

double function(double x){
	//Specify your desired function
	return(x*pow(e,x));
	//return(cos(x));
}

double centralDiff(double xValue, double stepSize){
	double differentiatedValue;

	differentiatedValue = (function(xValue+stepSize)-function(xValue-stepSize)) / (2.0*stepSize);
	return(differentiatedValue);

}

double richardson(int i, double Dih, double Dihh){

	double Dii;

	Dii = (pow(2.0,2.0*i)*Dih - Dihh) / (pow(2.0,2.0*i) - 1);

	return(Dii);

}

int main(){

	double D[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	double stepSize = 0.4;
	int counter = 1;
	int counter2;
	FILE *outfile;
	outfile = fopen("Ex_3.13.richardson.dat","w");
	fprintf(outfile,"h & $D_1$ & $D_2$ & $D_3$ & $D_4$ \\\\ \hline\n");

	// Build first column from central difference method
	for(counter=0;counter<4;counter++){

		stepSize = 0.4 / pow(2,counter);
		D[counter][0] = centralDiff(2.0,stepSize);
	
	}


	// Build table from Richardson extr.
	for(counter2=1;counter2<4;counter2++){

		for(counter=counter2-1;counter<4;counter++){

			D[counter][counter2] = richardson(counter2, D[counter][counter2-1], D[counter-1][counter2-1]);

		}

	}

	for(counter2=0;counter2<4;counter2++){

		fprintf(outfile,"%.3lf",0.4/pow(2,counter2));
		for(counter=0;counter<4;counter++){

			//printf("D_%i(%lf) = %lf\n",counter2+1,0.4/(pow(2,counter)),D[counter][counter2]);
			if(D[counter2][counter]<24 && D[counter2][counter]>21){

				fprintf(outfile," & %lf", D[counter2][counter]);
			
			}
			else{

				fprintf(outfile," & ");

			}

		}

		fprintf(outfile," \\\\ \hline\n");

	}

}
