#include <stdio.h>
#include <stdlib.h>
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


int main()
{
    float* u_b = create1dRandRHS(10);
    float* u0 = create1dRandRHS(10);
    float* u = create1dRandRHS(10);

    printf("Hello World\n");
    
    int itermax = 100000;
    int n = 10;
    
    for (int k=0; k<itermax; k++){
        for (int i=1; i<n-1; i++){
            u_b[i] = (u[i-1]+u[i+1])/2;
        }
	for (int i=0; i<n; i++){
		u[i] = u_b[i];
	}
    }
    print_vec(n,u_b);  print_vec(n,u);
    free(u_b); free(u); free(u0);
    return 0;
}










