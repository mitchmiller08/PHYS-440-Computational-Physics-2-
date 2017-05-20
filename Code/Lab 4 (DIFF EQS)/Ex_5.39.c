#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 4
// Ex. 5.39 4/17/12
// Numerical Differential Equations

int main(){

	double x,d1,d2,det,detOld=0,eig;
	double T,muo,mu,L,lambda,omgf,h,delta,alpha;
	double *a,*b,*c;
	int i,j,k,n,nl;
	n = 1000;
	a = (double*)malloc(n*sizeof(double));
	b = (double*)malloc(n*sizeof(double));
	c = (double*)malloc(n*sizeof(double));

	muo = 0.954;
	delta = 0.5;
	T = 1000;
	L = 1;
	nl = n-2;
	h = 1e-3;
	omgf = (pi/L)*sqrt(T/muo);

	alpha = 0.5*pow((pi/n),2);
	printf("%.15lf\t%.15lf\t%.15lf\n",h,omgf,alpha);
	
	for(i=0;i<n;i++){

		x = h*((double)i+1);
		//printf("%lf\n",x);
		mu = muo + (x - 0.5*L)*delta;
		a[i] = T/(h*h*mu);

		a[i] = a[i]/(2*T/(muo*h*h));
		b[i] = -2*a[i];
		c[i] = a[i];
		//printf("%lf\t%lf\t%lf\n",a[i],b[i],c[i]);
	
	}
	
	for(j=1;j<nl+1;j++){

		lambda = 50*alpha*(j)/nl;
		
		det = b[0] + lambda;
		//printf("%lf\n",det);

		if(n != 1){

			d1 = det;
			det = (b[1]+lambda)*d1-a[1]*c[0];
			//printf("%lf\n",det);

			if(n != 2){

				d2 = det;

				for(k=2;k<n;k++){

					det = (b[k]+lambda)*d2 - a[k]*c[k-1]*d1;
					d1 = d2;
					d2 = det;
					//printf("%lf\n",det);

				}

			}

		}


		if(det*detOld<0){

			eig = sqrt(2*T*lambda/(muo*h*h));
			eig = eig / omgf;
			printf("%.15lf\n",eig);

		}

		detOld = det;

	}

}
