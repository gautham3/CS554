#ifndef HELPERS_C
#define HELPERS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>


float *create1dZeroVec(int n){
    float *x;
    x = malloc(sizeof(float*) * n);
    for(int i = 0; i < n; i++){
        x[i] = 0;
    }
    return x;
}

/*Rand Vector creation*/
float *create1dRandRHS(int n){
    float *x;
    time_t t;
    x = malloc(sizeof(float*) * n);
    /* Intializes random number generator */
    srand((unsigned) time(&t));

    for(int i = 0; i < n; i++){
        x[i] = rand()%10;
    }
    return x;
}

void print_vec(int n, float* v)
{
  for (int i=0; i<n; i++)
    printf("%3f \n", v[i]);
  printf("\n");
}

void print_intvec(int n, int* v)
{
  for (int i=0; i<n; i++)
    printf("%3d \n", v[i]);
  printf("\n");
}

float VecErrNorm(int n, float* a, float* b)
{
  // Compare the L-infinity norm of sol a w.r.t. exact b
  float tmp, errnorm; 
  float max = 0; 
  float maxb = 0;
  for (int i=0; i<n; i++)
  {
    tmp = fabs(a[i] - b[i]); 
    max = (tmp>max) ? tmp : max;
    maxb = (fabs(b[i])>maxb) ? fabs(b[i]) : maxb;
  }
  maxb = (maxb==0)? 1: maxb;
  errnorm = max/maxb;
  //printf("%3f  , %3f , %3f", max,maxb,errnorm);printf("\n  ");
  return errnorm;      
}

float* VecAdd( float* a, float* b, float k, int n)
{
    float* c = create1dZeroVec(n);           
    for (int i=0; i<n; i++) {
            c[i] = a[i] + k*b[i];
    }       
    return c;
}

int* row_load_allot(int n, int ptot)
{
    int nwrks, offset, avrow, rows, lrow, roweq, rowrem;
    int* offsv = malloc(sizeof(int*) * ptot);    

    nwrks = ptot-1; // no. of workers
    avrow = floor((float)n/(float)nwrks); // average no. of rows of A each worker deals with      
    lrow = 0;
    roweq = avrow*nwrks;
    if (roweq<n){lrow=1; rowrem = n-roweq;} 
    offset = 0; offsv[0]= 0;

    for (int k=1; k<ptot; k++)
    {
        if(k>rowrem){lrow=0;}
        rows = avrow + lrow; 
        offset = offset + rows; 
        offsv[k] = offset; 
    }

    return offsv;
}

void communicate_1d_BCs(float a,float b, int rank, int p)
{
    if (rank==p){ MPI_Bcast(&b, 1, MPI_FLOAT, rank-1, MPI_COMM_WORLD);}
    if (rank==0){ MPI_Bcast(&a, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);}
    return;
}

float get1dLapMaxErr(int n, int rank, float* a, float* b)
{
    float err =  VecErrNorm(n, a, b); //local vector error
    float errmax;
    MPI_Reduce(&err, &errmax, 1, MPI_FLOAT, MPI_MAX,0, MPI_COMM_WORLD);
    return errmax;
}

#endif