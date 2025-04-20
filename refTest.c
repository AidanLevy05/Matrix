#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "MatrixC.h"

void divisor(const char *msg)
{
  printf("\n---- %s ----\n\n", msg);
}

// Global matrices
Matrix A_parallel, A_seq;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  srand(time(NULL) + rank);

  int N = 1000;
  initSize(&A_parallel, N, N);
  initSize(&A_seq, N, N);

  // Rank 0 generates the random matrix
  if (rank == 0)
  {
    setRandom(&A_parallel, 100);
    copyMatrix(&A_seq, &A_parallel);
  }

  // Broadcast A_parallel to all ranks
  MPI_Bcast(&(A_parallel.matrix[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ----------- Parallel ref() -----------
  if (rank == 0)
    printf("Running parallel ref()...\n");

  double start = MPI_Wtime();
  ref(&A_parallel);
  double end = MPI_Wtime();

  if (rank == 0)
  {
    printf("ref() done in %.3f seconds\n", end - start);
    divisor("ref");
  }

  // ----------- Sequential Seq_ref() -----------
  if (rank == 0)
  {
    printf("Running sequential Seq_ref()...\n");

    start = MPI_Wtime();
    Seq_ref(&A_seq);
    end = MPI_Wtime();

    printf("Seq_ref() done in %.3f seconds\n", end - start);
    divisor("Seq_ref");
  }

  MPI_Finalize();
  return 0;
}