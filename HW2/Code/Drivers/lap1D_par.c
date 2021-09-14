
#include "../Helpers/helpers.c"
#include "../Helpers/mpi_lap1d_solvers.c"


int main(int argc, char *argv[])
{
    int numprocs, rank, rc;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

    int N = atoi(argv[1]);  // size of full 1-D problem
    int ntimer = atoi(argv[2]); // No. of times lap1D solve done- time will be averaged
    int itermax = 10000; // No. of iterations for each solve
    int n = 0;
    int* offsvec =  row_load_allot(N,numprocs+1); // vector of offsets for each worker (here master also works)
    if (rank>=0) { n = offsvec[rank+1]-offsvec[rank]; } //Calc. no. of elements of u this worker rank deals with

    float* u_b = create1dRandRHS(n);
    float* u = create1dRandRHS(n);

    // create struct with data for parallel routines
    struct par_dat pmdat = par_struct_assign(offsvec,n,rank,N,numprocs); 

    //time 1D Laplacian solves
    lap1D_MPI_timer_output(ntimer,pmdat,u,u_b,itermax);

    //print_vec(n,usol);print_vec(n,u_b);
    free(u_b); free(u);
    MPI_Finalize();
    return 0;
}










