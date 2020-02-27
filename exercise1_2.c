#include <mpi.h>
#include <stdio.h>
#include <math.h>

// Run with: mpicc exercise2_2.c && mpiexec -n 2 ./a.out &&rm a.out
// or: mpicc exercise2_2.c && mpiexec --oversubscribe  -n 10 ./a.out &&rm a.out

int main(int argc, char** argv) {
    int n = 1000, i;
    double pi = 0.0, partial_sum_pi = 0.0, h = 1.0/n;
    MPI_Init(NULL, NULL);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (n % world_size != 0){
      printf("Problems with dividing the processes!\n");
      return -1;
    }

    int n_i = n / world_size;
    int lower_limit = n_i * world_rank, upper_limit = n_i * (world_rank + 1);
    printf("Rank %d: [%d, %d]: n_i = %d\n", world_rank, lower_limit, upper_limit, n_i);

    if (world_rank == 0){
      for (i = 1; i < world_size; ++i){
        MPI_Send(&n_i, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Recv(&partial_sum_pi, 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        pi += partial_sum_pi;
      }
      for (int i = lower_limit; i < upper_limit; ++i){
        pi += 4.0 / (1 + h*h*(i - 0.5)*(i - 0.5));
      }
      pi *= h;
      printf("Rank %d out of %d processors has pi = %lf\n", world_rank, world_size, pi);

    } else {
      MPI_Recv(&n_i, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      double partial_sum_pi_here = 0.0;
      for (int i = lower_limit; i < upper_limit; ++i){
        partial_sum_pi_here += 4.0 / (1 + h*h*(i - 0.5)*(i - 0.5));
      }

      MPI_Send(&partial_sum_pi_here, 1, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD);
      printf("Rank %d out of %d processors is computing..\n", world_rank, world_size);
    }

    MPI_Finalize();
    return 0;
}
