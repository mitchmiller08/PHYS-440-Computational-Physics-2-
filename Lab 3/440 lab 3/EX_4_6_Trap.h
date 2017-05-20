/*
 *  EX_4_6_Trap.h
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/23/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
 #include <stdio.h>
 #include <stdlib.h>
 
static double Pi = 3.141592653589793;
double functionInput(double xValue, double kVal)
{
	double functionVal = 0;
	functionVal = (1 / sqrt(1 - pow( kVal, 2 )*pow( sin(xValue), 2 )));
	return functionVal;
}