/*
 *  UniformDistribution.c
 *  440 lab 3
 *
 *  Created by Matthew Beck on 2/25/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
#include "rngs.h"
#include "rvgs.h"
   double Uniform(double a, double b)
/* =========================================================== 
 * Returns a uniformly distributed real number between a and b. 
 * NOTE: use a < b
 * ===========================================================
 */
{ 
  return (a + (b - a) * Random());
}

