/******************************************************************************************/
/*                                                                                        */
/*    GenericPotential.c                                                                  */
/*    Anthony Ruth                                                                        */
/*    Copyright 2011                                                                      */
/*    This program reads in a file named ForGenericPotentialFinder.txt                    */
/*    This file consists of locations in angstroms and potentials at those locations in eV*/
/*    The program then determines the energy eigenvalues that fit the potential file      */
/*    Once the program has determined an energy eigenvalue it creates the corresponding   */
/*    normalized wavefunction eigenvector. The energy is displayed and the wavefunction   */
/*    is printed to a file. The program is currently set to only print the first 4        */ 
/*    wavefunctions.																	  */
/*    Since the program uses Simpson's rule integration to normalize the wavefunction     */
/*    the number of points should be odd, or the wavefunction will not be properly        */
/*    nomalized.                                                                          */
/*    The program outputs both the energy and the calculated value of psi at the right    */
/*    endpoint. There is currently a bug in the program, so if the value of psi at the    */
/*    right endpoint times the width of the well is greater than 5, the printed           */
/*    wavefunctions are garbage, but the energies are still good.                         */ 
/******************************************************************************************/
#include <stdio.h>
#include <math.h>
//#include <io.h>
#include <stdlib.h>
static double hbarsquaredoverm = 7.61996348636199; // ev/angstrom squared 

void readinpotentials(double *potential, double *location, int numberoflines)
{
	/******************************************************/
	/* Forms a potential array and a location array       */
	/* for use in other parts of the program              */
	/* potential file is formated with locations in the   */
	/* first column in units of angstroms and potentials  */
	/* in the second column in units of eV                */
	/******************************************************/
	FILE *infile;
	if((infile=fopen("ForGenericPotentialFinder.txt", "r")) == NULL)
		 {
		printf("Cannot open file.\n");
		exit(1);
		}
		int index;
		for (index = 0;index<numberoflines ;index++ ) // read data from file into arrays
		{
		fscanf(infile,"%lg %lg",&location[index],&potential[index]);
		}
		fclose(infile);
}

int countnumberoflines()
{
	/******************************************************************/
	/* counts the number of lines in the file and returns it as an    */
	/* integer so the rest of the program can use it                  */
	/******************************************************************/
	FILE *infile;
	if((infile=fopen("ForGenericPotentialFinder.txt", "r")) == NULL)
		 {
		printf("Cannot open file.\n");
		exit(1);
		}
	double z1=-1.011111,z2=1.054e-12;// these numbers are random to ensure that values found in the file are not equal to them
	double V1=0.0000045,V2=1.424234;
	int numberoflinescounted = -1; //starts at -1 because last line of file is counted twice

	while(z1!=z2 || V1!=V2)  // count lines to determine size of array
	{
		numberoflinescounted++;
		z2 = z1;
		V2 = V1;
		fscanf(infile,"%lg %lg",&z1,&V1);
	}
	fclose(infile);
	return numberoflinescounted;
}

double calculatepsiinfinity(double *potentialforinffunction,double Eforinffunction, int numberoflinesforinffunction,double *location)
{
	/*************************************************************/
	/* calculates the wavefunction at the right endpoint and     */
	/* returns it as a double. This information is used to       */
	/* locate the energy eigenvalues.                            */
	/* Pi+1 = Pi + Pi' * dz                                      */
    /* Pi+1' = Pi' + Pi'' * dz                                   */
    /* -(h^2)/(2m) * Pi'' + V(zi)*Pi = En *Pi                    */
    /* Pi+1' = Pi' + (V(zi)-En)*2m/(h^2) * Pi * dz               */
    /* P0 = 0                                                    */
    /* P0' = 1e-4	this value seems to produce a smaller value  */
	/* of psi at the right endpoint than if P0' = 1				 */
	/*************************************************************/
	    double psi0 = 0.0, psi1;  
		double psi0prime = 1.0e-4,psi1prime;
		double potential = 0;
		double step = 0;
		int index;
		for (index = 1;index<numberoflinesforinffunction ;index++ )
		{
			step = location[index]-location[index-1];
			potential = potentialforinffunction[index-1];
			
			psi1 = psi0	+ psi0prime * step;
			psi1prime = psi0prime + 2.0/hbarsquaredoverm * (potential-Eforinffunction) * psi1 * step;
			//printf("step size %g  potential  %g  psi   %g    psi'    %g\n",step,potential,psi0,psi0prime);
			psi0 = psi1;
			psi0prime = psi1prime;
		}
		return psi1;
}

double rootfinder(double potentialforrootfinder[],double Eforrootfinder, int numberoflinesforrootfinder,double sigmaeforrootfinder,double *location)
{
	/********************************************************************/
	/* This function forms an interval over which an energy eigenvalue  */
	/* may reside. It guesses a number in that interval and uses the    */
	/* calculatepsiinfinity function to determine if the guessed number */
	/* should become the new left bound or new right bound of the       */
	/* interval. The interval becomes smaller and smaller until         */
	/* adequate accuracy has been reached or it becomes apparent        */
	/* that the interval has shrunk down to a single point, whichever   */
	/* comes first. It then returns the most accurate energy eigenvalue */
	/* found as a double                                                */
	/********************************************************************/
	double Elinearestimate1 = Eforrootfinder-sigmaeforrootfinder; 
	double Elinearestimate2 = Eforrootfinder,Elinearestimate3;
	double psiinf1forrootfinder=calculatepsiinfinity(potentialforrootfinder,Elinearestimate1,numberoflinesforrootfinder,location);
	double psiinf2forrootfinder=calculatepsiinfinity(potentialforrootfinder,Elinearestimate2,numberoflinesforrootfinder,location);
	double psiinf3forrootfinder;
	while ((fabs((Elinearestimate1- Elinearestimate2)/(Elinearestimate1 + Elinearestimate2))>1.12e-16))
	{	
		Elinearestimate3 = (Elinearestimate1 + Elinearestimate2)/2.0;
		//printf("%.17g %.17g %.17g %.17g %.17g %.17g\n",Elinearestimate1,Elinearestimate2,Elinearestimate3,psiinf1forrootfinder,psiinf2forrootfinder,psiinf3forrootfinder);
		psiinf3forrootfinder = calculatepsiinfinity(potentialforrootfinder,Elinearestimate3,numberoflinesforrootfinder,location);
		if (psiinf3forrootfinder * psiinf2forrootfinder > 0.0)
		{
			Elinearestimate2 = Elinearestimate3;
			psiinf2forrootfinder = psiinf3forrootfinder;
		}
		else
		{
			Elinearestimate1 = Elinearestimate3;
			psiinf1forrootfinder = psiinf3forrootfinder;
		}

	}
	return Elinearestimate3;
}
void fillpsiarray(double Eforpsiarray, double *potentialforpsiarray, int numberoflinesforpsiarray, double *psi, double *location)
{
	/*************************************************************/
	/*fills array with psi values                                */
	/* Pi+1 = Pi + Pi' * dz                                      */
    /* Pi+1' = Pi' + Pi'' * dz                                   */
    /* -(h^2)/(2m) * Pi'' + V(zi)*Pi = En *Pi                    */
    /* Pi+1' = Pi' + (V(zi)-En)*2m/(h^2) * Pi * dz               */
    /* P0 = 0                                                    */
    /* P0' = 1e-4	this value seems to produce a smaller value  */
	/* of psi at the right endpoint than if P0' = 1				 */
	/*************************************************************/
		double psi0 = 0.0,psi1;
		double psi0prime = 5,psi1prime;
		psi[0]=psi0;
		int index;
		double step = 0;
		double potential = 0;
		for (index = 1;index<numberoflinesforpsiarray ;index++ ) 
		{
			step = location[index]-location[index-1];
			potential = potentialforpsiarray[index-1];
			
			psi1 = psi0	+ psi0prime * step;
			psi1prime = psi0prime + 2.0/hbarsquaredoverm * (potential-Eforpsiarray) * psi1 * step;
			//printf("step size %g  potential  %g  psi   %g    psi'    %g\n",step,potential,psi0,psi0prime);
			psi0 = psi1;
			psi0prime = psi1prime;
			psi[index]=psi1;
		}
}
double simpsonsruleintegration(int numberofpoints, double *psi, double *location)
{
	/****************************************************************************************************/
	/*				Simpson's rule integration										                    */
	/*																				                    */
	/*			b	   (n-1)/2															                */
	/*			?f(x) ~  S  (f(2i*(b-a)/n+a)+4*f((2i+1)*(b-a)/n) + f((2i+2)*(b-a)/n))/3 * (b-a)/n       */
	/*			a		i=0															                    */
	/*         which is equal to													                  	*/
	/*         w =  (b-a)/n{1/3,4/3,2/3,4/3,2/3,.... ,2/3,4/3,1/3}                                      */
	/*		   f =  f{a,a+(b-a)/n,...,a + (n-1)(b-a)/n,b}                                               */
	/*																					                */
	/*			b						 n                                                              */
	/*			?f(x)	 =				 S  w[i]*f[i]												    */
	/*			a						 i = 0	                                                        */
	/****************************************************************************************************/
	double sum = 0.0,increment;
	int index;
	for (index = 0; index < numberofpoints; index++)
	{
		
		if (index % 2 == 0)
		{
			if (index == 0)
			{
				increment = location[index+1]-location[index];
				sum += increment * psi[index] * psi[index]/3.0;//left and right endpoints are at even indeces and weighted 1/3
			}
			else
			{
				if (index == numberofpoints-1)
				{
					increment = location[index]-location[index-1]; // location[index+1] does not exist here
					sum += increment * psi[index] * psi[index]/3.0;
				}
				else
				{
					increment = location[index+1]-location[index];
					sum += 2.0 * increment * psi[index] * psi[index]/3.0;//Even indeces that are not endpoints are weighted 2/3
				}
			}
		}
		else
		{
			{
				increment = location[index+1]-location[index];
				sum += 4.0 * increment * psi[index] * psi[index]/3.0;//Odd indeces are weighted 4/3
			}
		}
	}
	return sum;
}

void normalizewavefunctions(double *psi, int numberoflinesfornormalizing, double *location)
{
	/****************************************************************/
	/*  Uses Simpson's rule integration to determine the integral of*/
	/*  psi^2. Then divides the entire array to satisfy:            */
	/*	      b                                                     */
	/*        ?psi^2 = 1                                            */
	/*        a                                                     */
	/****************************************************************/
	double sum = simpsonsruleintegration(numberoflinesfornormalizing, psi, location);
	sum = pow(sum,0.5);
	int index;
	for (index = 0;index < numberoflinesfornormalizing ;index++ )
	{
		psi[index] = psi[index]/sum;
	}
}

void printtofile(double *psi, double *location, int eigennumber, int numberoflinesforprinting)
{
	/*********************************************************************/
	/* Prints the first four wavefunctions to files. The other           */
	/* wavefunctions are thrown in the garbage.txt file                  */
	/*********************************************************************/
			FILE *outfile;
			int index;
			if (eigennumber == 1)
			{
				outfile = fopen("GenericPsi1.txt", "w");
			}
			if (eigennumber == 2)
			{
				outfile = fopen("GenericPsi2.txt", "w");
			}
			if (eigennumber == 3)
			{
				outfile = fopen("GenericPsi3.txt", "w");
			}
			if (eigennumber == 4)
			{
				outfile = fopen("GenericPsi4.txt", "w");
			}
			if (eigennumber >= 5)
			{
				outfile = fopen("Garbage.txt","w");
			}
			for (index = 0;index < numberoflinesforprinting ;index++ ) //prints psi and location to a file so they can be graphed
			{
				fprintf(outfile,"%g   %g\n",location[index],psi[index]);
			}
			fclose(outfile);
}

void main()
{
		int numberoflines = countnumberoflines();
		printf("The number of lines in the file is %d \n",numberoflines);
		printf("The program outputs both the energy and the calculated value of psi at the right \n"); 
		printf("endpoint. There is currently a bug in the program, so if the value of psi at the    \n"); 
		printf("right endpoint times the width of the well is greater than 5, the printed    \n"); 
        printf("wavefunctions are garbage, but the energies are still good. \n");
		// prints the message in the top of the file which describes and determines when the
		// bug in the program has mutilated the wavefunctions                      
		double position[numberoflines];
		double potential[numberoflines]; //initialize arrays and pointers
		double wavefunctions[numberoflines];
		double *psi,*location, *potentialenergy;
		location = position;
		psi = wavefunctions;
		potentialenergy = potential;
		readinpotentials(potentialenergy,location,numberoflines);
		double psiinf1,psiinf2;
		double E,dE,Erootfinder;
		double sum;
		int numberoffoundvalues = 0;
		FILE *outfile;
		dE = 1.0e-4; //runs into problems if the difference between the energies of 2 adjacent eigenfunctions is <1e-4 eV, but otherwise makes program run faster
		psiinf1 = 0.0;
		psiinf2 = 0.0;
		for (E = dE;E < 1000; E+=dE )
	{

		psiinf2 = calculatepsiinfinity(potential,E,numberoflines,location);
		//printf("E= %g ev psi(-infinity)= %g\n",E,psiinf2);
		
		if (psiinf1*psiinf2<0) //detects when the correct energy has been crossed
		{
			Erootfinder = rootfinder(potential,E,numberoflines,dE,location); //finds a close approximation to the correct energy
			numberoffoundvalues++;
			fillpsiarray(Erootfinder,potential,numberoflines,psi,location); // fills wavefunction array at the correct eigenvalue of energy
			//normalizewavefunctions(psi,numberoflines,location);             // normalizes said array
			double plusmachine = calculatepsiinfinity(potentialenergy,Erootfinder+ Erootfinder*1.11e-16,numberoflines,location);
			double minusmachine = calculatepsiinfinity(potentialenergy,Erootfinder - Erootfinder*1.11e-16,numberoflines,location);
			printf("E= %g eV psi(right endpoint)= %g  width of well = %g\n",Erootfinder,psi[numberoflines-1],location[numberoflines-1]-location[0]); //displays energy eigenvalue and information that determines if wavefunctions are usable or not
			printf("E= %g eV psi(right endpoint)= %g  width of well = %g\n",Erootfinder + Erootfinder*1.11e-16,plusmachine ,location[numberoflines-1]-location[0]);
			printf("E= %g eV psi(right endpoint)= %g  width of well = %g\n",Erootfinder - Erootfinder*1.11e-16,minusmachine,location[numberoflines-1]-location[0]);
			printtofile(psi,location,numberoffoundvalues,numberoflines); //prints wavefunction and location arrays to file for graphing
		}
		psiinf1 = psiinf2;
		
	}	
}
