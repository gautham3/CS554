#ifndef HELPERS_C
#define HELPERS_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


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

#endif