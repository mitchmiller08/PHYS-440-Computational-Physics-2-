//
//  main.c
//  HeAtom
//
//  Created by Mitchell Miller on 3/23/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double bohrRadius = .5291772083; //angstroms
double Pi = 3.141592653589793;
double E = 2.718281828459045;
double eCharge = 14.39964439;
double Z = 2;
int numberOfLoops = 10000000;

double oneSState(double x,double y,double z){
    
    double funcValue;
    double exp = 0;
    double radius = sqrt(x*x+y*y+z*z);
    
    exp = -2*Z*radius/bohrRadius;
    funcValue = pow(E,exp);
    
    return(funcValue);
    
}

double twoSState(double x,double y,double z){
    
    double funcValue;
    double exp;
    double coef;
    double expTerm;
    
    double radius = sqrt(x*x+y*y+z*z);
    
    coef = pow(1. - 0.5*Z*radius/bohrRadius,2);
    exp = -Z*radius/bohrRadius;
    expTerm = pow(E,exp);
    funcValue = coef*expTerm;
    
    return(funcValue);
}

double twoPState(double x,double y,double z){
    
    double funcValue;
    double exp;
    double coef;
    double expTerm;
    
    double radius = sqrt(x*x+y*y+z*z);
    
    exp = -2*Z*radius/bohrRadius;
    coef = radius*radius;
    expTerm = pow(E,exp);
    funcValue = coef*expTerm;
    
    return(funcValue);
}

double K1s2sIntegral(double x1,double y1,double z1,double x2,double y2,double z2){
    
    double funcValue;
    double coef;
    double exp;
    double radius1 = sqrt(x1*x1+y1*y1+z1*z1);
    double radius2 = sqrt(x2*x2+y2*y2+z2*z2);
    
    coef = (1-(.5)*Z*radius2/bohrRadius)*(1-(.5)*Z*radius1/bohrRadius);
    exp = -Z*(radius1+radius2)* (3/ (2*bohrRadius));
    funcValue = coef*pow(E,exp);
    
    return(funcValue);
}

double K1s2pR(double x,double y,double z){
    
    double funcValue;
    double exp;
    
    double radius = sqrt(x*x+y*y+z*z);
    exp = -(3*Z*radius/(2*bohrRadius));
    funcValue = z*pow(E,exp);
    
    return(funcValue);
}

double repulsion(double x1,double y1,double z1,double x2,double y2,double z2){
    
    double funcValue;
    double radius = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    
    funcValue = 1 / fabs(radius);
    
    return(funcValue);
}

double J1S2S(double low, double high){
    
    double funcValue=0;
    int counter;
    
    double x1=0;
    double x2=0;
    double y1=0;
    double y2=0;
    double z1=0;
    double z2=0;
    
    double oneS=0;
    double twoS=0;
    double repul=0;
    
    double width = high - low;
    
    double coef1S = pow((Z / bohrRadius),3) / Pi;
    double coef2S = pow((Z / bohrRadius),3) / (8*Pi);
    
    for(counter=0;counter<numberOfLoops;counter++){
        
        x1 = (drand48() - 0.5) * width;
        y1 = (drand48() - 0.5) * width;
        z1 = (drand48() - 0.5) * width;
        x2 = (drand48() - 0.5) * width;
        y2 = (drand48() - 0.5) * width;
        z2 = (drand48() - 0.5) * width;
        
        oneS = fabs(oneSState(x1, y1, z1));
        twoS = fabs(twoSState(x2, y2, z2));
        repul = repulsion(x1, y1, z1, x2, y2, z2);
        
        funcValue += oneS * twoS * repul;
        
    }
    
    funcValue /= numberOfLoops;
    funcValue *= coef1S*coef2S*eCharge*pow(width,6);
    return(funcValue);
    
}

double K1S2S(double low, double high){
    
    double funcValue=0;
    int counter;
    
    double x1=0;
    double x2=0;
    double y1=0;
    double y2=0;
    double z1=0;
    double z2=0;
    
    double repul;
    double kInteg;
    double zBohr = Z / bohrRadius;
    double coef1S = (1 / (8*pow(Pi,2)))*pow(zBohr,6);
    
    double width = high - low;
    
    
    for(counter=0;counter<numberOfLoops;counter++){
        
        x1 = (drand48() - 0.5) * width;
        y1 = (drand48() - 0.5) * width;
        z1 = (drand48() - 0.5) * width;
        x2 = (drand48() - 0.5) * width;
        y2 = (drand48() - 0.5) * width;
        z2 = (drand48() - 0.5) * width;
        
        kInteg = K1s2sIntegral(x1,y1,z1,x2,y2,z2);
        repul = repulsion(x1, y1, z1, x2, y2, z2);
        funcValue += kInteg*repul;
        
    }
    
    funcValue = funcValue / numberOfLoops;
    funcValue = funcValue * eCharge*pow(width,6)*coef1S;
    return(funcValue);
    
}

double J1S2P(double low, double high){
    
    double funcValue=0;
    int counter;
    double width = high - low;
    
    double x1=0;
    double x2=0;
    double y1=0;
    double y2=0;
    double z1=0;
    double z2=0;
    
    double zBohr = Z / bohrRadius;
    double repul=0;
    double oneS=0;
    double twoP=0;
    
    double coef1S = pow(zBohr,3) / Pi;
    double coef2P = (1 / (64*pow(Pi,2)))*pow(zBohr,8);
    
    for(counter=0;counter<numberOfLoops;counter++){
        
        x1 = (drand48() - 0.5) * width;
        y1 = (drand48() - 0.5) * width;
        z1 = (drand48() - 0.5) * width;
        x2 = (drand48() - 0.5) * width;
        y2 = (drand48() - 0.5) * width;
        z2 = (drand48() - 0.5) * width;
     
        oneS = oneSState(x1, y1, z1);
        twoP = twoPState(x2, y2, z2);
        repul = repulsion(x1, y1, z1, x2, y2, z2);
        
        funcValue += oneS*twoP*repul;
        
    }
    
    funcValue = funcValue / numberOfLoops;
    funcValue = funcValue * coef1S*coef2P*eCharge*pow(width,6);    // May need to add 'coef1S' here
    
    return(funcValue);
}

double K1S2P(double low, double high){
    
    double funcValue=0;
    int counter;
    double width = high - low;
    
    double x1=0;
    double x2=0;
    double y1=0;
    double y2=0;
    double z1=0;
    double z2=0;
    
    double zBohr = Z / bohrRadius;
    double repul;
    double r1;
    double r2;
    double coef2P = (1 / (32*pow(Pi,2)))*pow(zBohr,8);

    for(counter=0;counter<numberOfLoops;counter++){
        
        x1 = (drand48() - 0.5) * width;
        y1 = (drand48() - 0.5) * width;
        z1 = (drand48() - 0.5) * width;
        x2 = (drand48() - 0.5) * width;
        y2 = (drand48() - 0.5) * width;
        z2 = (drand48() - 0.5) * width;
        
        r1 = K1s2pR(x1, y1, z1);
        r2 = K1s2pR(x2, y2, z2);
        repul = repulsion(x1, y1, z1, x2, y2, z2);
        
        funcValue += r1*r2*repul;
        
    }
    
    funcValue = funcValue / numberOfLoops;
    funcValue = funcValue * eCharge*pow(width,6)*coef2P;
    
    return(funcValue);
    
}

int main ()
{

    double lower = -5;
    double upper = 5;
    double j1s2s=0;
    double k1s2s=0;
    double j1s2p=0;
    double k1s2p=0;
    
    int counter;
    int loops=1000;
    
    FILE *outfile;
    outfile = fopen("HeAtom.dat","w");
    
    srand48( time(NULL));
    
    for(counter=0;counter<loops;counter++){
        
        printf("Begin Loop\t");
        j1s2s += J1S2S(lower, upper);
        j1s2p += J1S2P(lower, upper);
        
        k1s2s += K1S2S(lower, upper);
        k1s2p += K1S2P(lower, upper);
        
        printf("Loop %i Complete\n",counter);
        
    }
    
    j1s2s = j1s2s / loops;
    j1s2p = j1s2p / loops;
    k1s2s = k1s2s / loops;
    k1s2p = k1s2p / loops;
    
    printf("J1S2S = %lf\n",j1s2s);
    printf("K1S2S = %lf\n",k1s2s);
    printf("J1S2P = %lf\n",j1s2p);
    printf("K1S2P = %lf\n",k1s2p);
    
    fclose(outfile);
    
    return 0;
}

