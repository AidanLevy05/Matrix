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

  // Rank 0 sets data
  if (rank == 0)
  {
    setRandom(&A_parallel, 100);
    copyMatrix(&A_seq, &A_parallel);
  }

  // Broadcast to all processes
  MPI_Bcast(&(A_parallel.matrix[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ----------- Parallel rref() -----------
  if (rank == 0)
    printf("Running parallel rref()...\n");

  double start = MPI_Wtime();
  rref(&A_parallel);
  double end = MPI_Wtime();

  if (rank == 0)
  {
    printf("rref() done in %.3f seconds\n", end - start);
    divisor("rref");
  }

  // ----------- Sequential Seq_rref() -----------
  if (rank == 0)
  {
    printf("Running sequential Seq_rref()...\n");

    start = MPI_Wtime();
    Seq_rref(&A_seq);
    end = MPI_Wtime();

    printf("Seq_rref() done in %.3f seconds\n", end - start);
    divisor("Seq_rref");
  }

  MPI_Finalize();
  return 0;
}
