#ifndef MPI_LAP1D_SOLVERS_C
#define MPI_LAP1D_SOLVERS_C

#include "helpers.c"


struct par_dat {
    int* offsetv_d ; /* used to determine start index of element for each worker */
    int rows_d ; /* no. of elements for each worker */
    int rank_d ;
    int n_d ;
    int nwrks_d ;
};

struct par_dat par_struct_assign(int *offsvec, int rows, int rank, int N, int nwrks)
{
    struct par_dat pmdat;
    pmdat.offsetv_d = offsvec;
    pmdat.rows_d = rows;
    pmdat.rank_d = rank;
    pmdat.n_d = N;
    pmdat.nwrks_d = nwrks; 
    return pmdat;
} 

struct solret {
    int iter;
    float* x;
};

void create_lap1d_sol(float* usol,float ul,float ur,int N,int n, int offs)
{
    float dsol = (ur-ul)/(N-1);
    usol[0] = ul+offs*dsol;
    for (int i=1; i<n; i++){
        usol[i] = usol[i-1]+dsol;
    }
    return;
}

struct solret mpiLap1Dsolve(struct par_dat pmd, float* u, float* u_b, int itermax)
{
    float ul, ur, ul1, ur1, ulr, url;
    int rank = pmd.rank_d;
    int numprocs = pmd.nwrks_d;
    int n = pmd.rows_d;
    MPI_Status status; 

    struct solret sol;
    for (int k=0; k<itermax; k++){

            ul = u[0]; 
            ur = u[n-1];
            ulr = 0; url = 0;
            ul1 = 2*u[0]-u[1]; // ensure Dirichlet BCs, for end points where not overwritten
            ur1 = 2*u[n-1]-u[n-2];
            if (rank>0){ MPI_Send(&ul, 1, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD);}
            if (rank<numprocs-1){ MPI_Send(&ur, 1, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD);}
            if (rank<numprocs-1){ MPI_Recv(&ur1, 1, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);}
            if (rank>0){ MPI_Recv(&ul1, 1, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);}

            for (int i=1; i<n-1; i++){
                u_b[i] = (u[i-1]+u[i+1])/2;
            }
            if (n>1){ u_b[0] = (ul1+u[1])/2; u_b[n-1] = (ur1+u[n-2])/2;}
            else {
                if ((rank==0)||(rank==numprocs-1)){u_b[0] = u[0];}
                else {u_b[0] = (ul1+ur1)/2;}
            }
            
            for (int i=0; i<n; i++){
                u[i] = u_b[i];
            }
        }
    sol.x = u_b;
    sol.iter = itermax;
    return sol;
}

void lap1D_MPI_timer_output(int ntimer, struct par_dat pmd, float* u, float* u_b, int itermax)
{
    
    int rank = pmd.rank_d;
    int numprocs = pmd.nwrks_d;
    int N = pmd.n_d;
    int nwrks = pmd.nwrks_d;
    int n = pmd.rows_d;
    int* offsvec = pmd.offsetv_d;
    clock_t beg, end; 
    struct solret sol;

    float ul = u[0];
    float ur = u[n-1];
    if (rank==0) {printf("rank:%d, left BC:%f\n",rank, ul); }
    if (rank==numprocs-1) {printf("rank:%d, right BC:%f\n",rank, ur); }
    float* usol = create1dZeroVec(n);


    // Start Timer
    beg = clock();

    for(int j=0; j<ntimer; j++){sol = mpiLap1Dsolve(pmd,u,u_b,itermax);}
    u_b = sol.x;
    int iters = sol.iter;

    // End timer
    end = clock();
    float t_tot = 1.0*(end-beg)/CLOCKS_PER_SEC;

    //check error in solution
    communicate_1d_BCs(ul,ur,rank,numprocs);
    create_lap1d_sol(usol,ul,ur,N,n,offsvec[rank]);
    //float err =  VecErrNorm(n, u_b, usol); //local vector error
    float errmax = get1dLapMaxErr(n, rank, u_b, usol); // print global vector error

    //Output
    if(rank==0){  
        printf("1D Laplacian solve:- \n");        
        printf("Size : %d on %d ranks\n",N,nwrks);
        printf("iters  %d \n", iters);  
        printf("Total time taken: %f \n", t_tot/ntimer);
        printf("Solution Error Norm: %f\n",errmax);
        printf("-------------------------------------------------------\n\n");
    }
    free (usol);
    return;
}






#endif