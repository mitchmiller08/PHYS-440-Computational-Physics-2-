//
//  Fig_3.6.c
//  Lab 4
//
//  Created by Mitchell Miller on 4/11/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define pi 3.141592654
#define e  2.71828182845905

// Mitchell Miller
// PHYS 440 Lab 4
// Fig. 3.6 4/17/12
// Numerical Differential Equations

double func(double theta,double angVel,double q,double fD,double omegaD,double t){
    
    double funcValue;
    
    funcValue = -1.0 * sin(theta) - q*angVel + fD*sin(omegaD*t);
    
    return(funcValue);
}

int main(){
    
	// Initialize variables
	int ncount=0,counter;
	double theta=0.2;
	double angVel=0;
	double q=0.5;
	double fD=1.465;
	double omegaD=2./3;
	double dt=0.01;
	double t=0;
	double k1v,k2v,k3v,k4v;
    	double k1x,k2x,k3x,k4x;
	FILE *outfile;
	outfile = fopen("Fig_3.10.1465.dat","w");
	
	for(counter=0;counter<10000;counter++){
        
		ncount++;
        
		// Begin RK 4 block

		k1v = func(theta,angVel,q,fD,omegaD,t)*dt;
   		k1x = angVel*dt;
    		t = t+0.5*dt;
    
	    	k2v = func(theta+0.5*k1x,angVel+0.5*k1v,q,fD,omegaD,t)*dt;
	    	k2x = (angVel+0.5*k1v)*dt;
	    	k3v = func(theta+0.5*k2x,angVel+0.5*k2v,q,fD,omegaD,t)*dt;
	    	k3x = (angVel+0.5*k2v)*dt;
	    	t = t+0.5*dt;
	    
	    	k4v = func(theta+k3x,angVel+k3v,q,fD,omegaD,t)*dt;
	    	k4x = (angVel+k3v)*dt;
	    
	    	angVel = angVel + (k1v + 2.*k2v + 2.*k3v + k4v)/6.;
	    	theta = theta + (k1x + 2.*k2x + 2.*k3x + k4x)/6.;
	    	
	    	if(theta >= pi){
			theta = theta - 2.0*pi;
	    	}
	    
	    	if(theta <= -1.*pi){
			theta = theta + 2.0*pi;
	    	}

		// End RK 4 block
		
		fprintf(outfile,"%.15lf\t%.15lf\t%.15lf\n",t,theta,angVel);
        
	}
    
}
