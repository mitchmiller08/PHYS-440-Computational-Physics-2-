#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define pi 3.141592654
//#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 4
// Fig. 5.10 5/4/12
// Lab 5

/* Define complex double data type */
typedef complex<double> dcomp;

int main(){

	dcomp i, a[600],u,w,t;
	int i,k,inv,n;
	i = dcomp (0., 1.);
	double r1,ai,spec,time;
	double ang,pi;
	double t[1001],x[1001],v[1001];
	int m,n1,nd2,l,le,le1,ip;

	FILE *infile;
	infile = fopen("van.dat","r");
	FILE *outfile2;
	outfile2 = fopen("Fig_6.18.spec.dat","w");

	k=8;
	n=pow(2,k);

	for(i=0;i<n;i++){

		fscanf(infile,"%lf\t%lf\t%lf", t[i],x[i],v[i]);
		printf("%lf\t%lf\n",t[i],x[i]

	}




}
