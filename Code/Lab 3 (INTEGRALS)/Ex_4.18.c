#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 3
// Exc. 4.18 2/28/12
// Gaussian Quadrature Method

double func(double x){

	double funcValue;
	funcValue = x*x*x/(1.-pow(e,-x));
	return(funcValue);

}

double map(double x){

	double mapValue;
	double a = 0;
	double b = 1;
	mapValue = (b-a)/2.*x + (a+b)/2.;
	return(mapValue);

}

double quad2(){

	double integral=0;
	int counter=0;
	int order=2;
	double x[2] = {5.8578643762690495e-1,3.4142135623730950   };
	double w[2] = {8.5355339059327376e-1,1.4644660940672624e-1};

	for(counter=0;counter<order;counter++){

		integral = integral + func(x[counter])*w[counter];

	}
	
	return(integral);

}

double quad4(){

	double integral=0;
	int counter=0;
	int order=4;
	double x[4] = {3.2254768961939231e-1,1.7457611011583466   ,4.5366202969211280   ,9.3950709123011331   };
	double w[4] = {6.0315410434163360e-1,3.5741869243779969e-1,3.8887908515005384e-2,5.3929470556132745e-4};

	for(counter=0;counter<order;counter++){

		integral = integral + func(x[counter])*w[counter];

	}

	return(integral);

}

double quad6(){

	double integral=0;
	int counter=0;
	int order=6;
	double x[6] = {2.2284660417926069e-1,1.1889321016726230   ,2.9927363260593141   ,5.7751435691045105   ,9.8374674183825899   ,1.5982873980601702e1 };
	double w[6] = {4.5896467394996359e-1,4.1700083077212099e-1,1.1337338207404498e-1,1.0399197453149075e-2,2.6101720281493206e-4,8.9854790642962124e-7};

	for(counter=0;counter<order;counter++){

		integral = integral + func(x[counter])*w[counter];

	}

	return(integral);

}

double quad8(){

	double integral=0;
	int counter=0;
	int order=8;
	double x[8] = {1.7027963230510100e-1,9.037017769937991e-1,2.2510866298661307   ,4.2667001702876588   ,7.0459054023934657   ,1.0758516010180995e1 ,1.5740678641278005e1 ,2.2863131736889296e1 };
	double w[8] = {3.6918858934163753e-1,4.1878678081434296e-1,1.7579498663717181e-1,3.3343492261215652e-2,2.7945362352256725e-3,9.0765087733582131e-5,8.4857467162725315e-7,1.0480011748715104e-9};

	for(counter=0;counter<order;counter++){

		integral = integral + func(x[counter])*w[counter];

	}

	return(integral);

}

int main(){

	double integral2,integral4,integral6,integral8;

	integral2 = quad2();
	integral4 = quad4();
	integral6 = quad6();
	integral8 = quad8();

	printf("N = 2\t%.15lf\n",integral2);
	printf("N = 4\t%.15lf\n",integral4);
	printf("N = 6\t%.15lf\n",integral6);
	printf("N = 8\t%.15lf\n",integral8);

}
