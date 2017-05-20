#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Lambda(int i, int Degree, double XMin, double XMax, double XIncrement, double* DataPointSubIPointer, double* InterpolationPointer, int NumberOfIntervals);

int main()
{
  printf("\nXMax-Xmin must be a multiple of XIncrement!\n");

  int NumberOfDataPoints, Degree, NumberOfIntervals;

  double XMin, XMax, XIncrement;

  FILE *Data = fopen("DataToBeFitted.txt", "r");
  FILE *FittedData = fopen("FittedData.txt", "a+");

  
  if(Data == NULL)
  {
     printf("error\n");
     return 0;
  }


  printf("\nEnter number of data points to be fitted: ");

  scanf("%d", &NumberOfDataPoints);
  printf("\nEnter the degree of the polynomial you wish to fit your data to: ");

  scanf("%d", &Degree);

  Degree=Degree+1;



  printf("\nEnter X Min: ");

  scanf("%lf", &XMin);
  
  printf("\nEnter X Max: ");
  
  scanf("%lf", &XMax);

  printf("\nEnter X Increment: ");

  scanf("%lf", &XIncrement);

  printf("\nEnter (Xmin-XMax)/XIncrement (Integer): ");

  scanf("%d", &NumberOfIntervals);

  int InterpolationSize;

  int i,j,a,b;

  InterpolationSize = Degree*NumberOfIntervals;

  double* InterpolationPointer;

  double** InterpolationMatrix;

  InterpolationMatrix = (double **) malloc(Degree * sizeof(double *));

  for(i = 0; i < Degree; i++)
  {
    InterpolationMatrix[i] = (double *) malloc(NumberOfIntervals * sizeof(double));
  }





  double PolynomialOfX;

  

  double* FunctionOfDataPoints=(double*)malloc(NumberOfDataPoints*sizeof(double));

  double* DataPointSubIPointer;

  InterpolationPointer = (double*)malloc(InterpolationSize*sizeof(double));

  DataPointSubIPointer = (double*)malloc(NumberOfDataPoints*sizeof(double));

  
  

  
  
  for(i=0;i<=NumberOfDataPoints-1;i++)
  {
    
    fscanf(Data, "%le %le\n", &DataPointSubIPointer[i], &FunctionOfDataPoints[i]); 
    
  }

  fclose(Data);

  

  for(i=0;i<Degree;i++)
  {
    Lambda(i, Degree, XMin, XMax, XIncrement, DataPointSubIPointer, InterpolationPointer, NumberOfIntervals);
  }

  for (a=0;a<NumberOfIntervals;a++)
  {
    for (b=0;b<Degree;b++)
    {
      InterpolationPointer[b+Degree*a] = InterpolationMatrix[a][b];
    }
  }

  for(j=0;j<NumberOfIntervals;j++)
  {
    PolynomialOfX=0;
    for(i=0;i<Degree;i++)
    {
      PolynomialOfX = PolynomialOfX + FunctionOfDataPoints[i]*InterpolationMatrix[i][j];
    }

    fprintf(FittedData,"%lf %lf\n" , XMin+j*XIncrement, PolynomialOfX);
  }
  
  fclose(FittedData);
  free(FunctionOfDataPoints);
  free(DataPointSubIPointer);
  free(InterpolationPointer);

  /* Deallocate memory of 2D array. */
  for(i = 0; i < Degree; i++)
  {
    free(InterpolationMatrix[i]);
  }
  free(InterpolationMatrix);

  return 0;

}//end main




void Lambda(int i, int Degree, double XMin, double XMax, double XIncrement, double* DataPointSubIPointer, double* InterpolationPointer, int NumberOfIntervals)
{
  int j,a,b,X;
  int n=0;

  double** InterpolationMatrixForm;

  InterpolationMatrixForm = (double **) malloc(Degree * sizeof(double *));

  for(a = 0; a < Degree; a++)
  {
    InterpolationMatrixForm[a] = (double *) malloc(NumberOfIntervals * sizeof(double));
  }


  //put into a 2d array of Degree x (Xmax-Xmin)/XIncrement for convenience 
  for (a=0;a<(XMax-XMin)/XIncrement;a++)
  {
    for (b=0;b<Degree;b++)
    {
      InterpolationPointer[b+Degree*a] = InterpolationMatrixForm[a][b];
    }
  }

  
  for(X=XMin;X<XMax;X+=XIncrement)
  {
    for(j=0;j<Degree;j++)
    {
      if(j!=i)
      {
        InterpolationMatrixForm[i][n]=0;
      }
    }
    n=n+1;
  }

  n=0;
  for(X=XMin;X<XMax;X+=XIncrement)
  {
    for(j=0;j<Degree;j++)
    {
      if(j!=i)
      {
        InterpolationMatrixForm[i][n]=InterpolationMatrixForm[i][n]*((X-DataPointSubIPointer[j])/(DataPointSubIPointer[i]-DataPointSubIPointer[j]));
      }
    }
    n=n+1;
  }

  for (a=0;a<(XMax-XMin)/XIncrement;a++)
  {
    for (b=0;b<Degree;b++)
    {
      InterpolationMatrixForm[a][b] = InterpolationPointer[b+Degree*a];
    }
  }

  /* Deallocate memory of 2D array. */
  for(a = 0; a < Degree; a++)
  {
    free(InterpolationMatrixForm[a]);
  }
  free(InterpolationMatrixForm);
}
