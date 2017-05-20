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

double func(theta1,angVel1,q1,fD1,omegaD1,t1){
    
    double funcValue;
    
    funcValue = -sin(theta1) - q1*angVel1 + fD1*sin(omegaD1*t1);
    
    return(funcValue);
}

int main(){
    
	// Initialize variables
	int ncount=0,counter;
	double theta=0.2;
	double angVel=0.;
	double q=0.5;
	double fD=0.5;
	double omegaD=2./3;
	double dt=0.04;
	double t=0;
	double k1v,k2v,k3v,k4v;
    double k1x,k2x,k3x,k4x;
	FILE *outfile;
	outfile = fopen("/Users/Mitch/Dropbox/PHYS 440/Lab 4/Lab 4/Lab 4/Fig_3.6.dat","w");
	
	for(counter=0;counter<1000;counter++){
        
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
	    
        k4v = func(theta+0.5*k3x,angVel+0.5*k3v,q,fD,omegaD,t)*dt;
        k4x = (angVel+0.5*k3v)*dt;
	    
        angVel = angVel + (k1v + 2.*k2v + 2.*k3v + k4v)/6.;
        theta = theta + (k1x + 2.*k2x + 2.*k3x + k4x)/6.;
	    
        if(theta >= pi){
			theta = theta - 2.0*pi;
        }
	    
        if(theta <= -pi){
			theta = theta + 2.0*pi;
        }
        
		// End RK 4 block
		
		fprintf(outfile,"%lf\t%lf\t%lf\n",t,theta,angVel);
        
	}
    
}