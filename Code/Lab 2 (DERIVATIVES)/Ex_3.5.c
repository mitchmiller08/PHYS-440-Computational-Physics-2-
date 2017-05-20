#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Mitchell Miller
// PHYS 440 Lab 2
// Exc. 3.5 2/11/12
// Trisolve Method

int main(){

	int counter;
	int diagLength = 5;
	double aMat[5] = {-1.0,-1.0,-1.0,-1.0,-1.0};
	double bMat[5] = {2.,2.,2.,2.,2.};
	double cMat[5] = {-1.0,-1.0,-1.0,-1.0,-1.0};
	double rMat[5] = {0.,1.,2.,3.,4.};
	double rhoMat[5] = {0.,0.,0.,0.,0.};
	double betaMat[5] = {0.,0.,0.,0.,0.};

	betaMat[0] = bMat[0];
	rhoMat[0] = rMat[0];

	for(counter = 1; counter < diagLength; counter++){

		betaMat[counter] = bMat[counter] - (aMat[counter] * cMat[counter-1])/betaMat[counter-1];
		rhoMat[counter] = rMat[counter] - (aMat[counter] * rhoMat[counter-1])/betaMat[counter-1];

	}

	rhoMat[diagLength-1] = rhoMat[diagLength-1] / betaMat[diagLength-1];

	for(counter = 1; counter < diagLength; counter++){

		rhoMat[diagLength-1-counter] = (rhoMat[diagLength-1-counter] - cMat[diagLength-1-counter] * rhoMat[diagLength-counter]) / betaMat[diagLength-1-counter];

	}

	for(counter = 0; counter < diagLength; counter++){

		printf("X_%i = %lf\n", counter, rhoMat[counter]);

	}

}
