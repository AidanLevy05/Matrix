#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "MatrixC.h"

void divisor(const char *msg)
{
  printf("\n---- %s ----\n\n", msg);
}

Matrix A, B, C_parallel, C_seq;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  srand(time(NULL) + rank);

  int N = 1000;
  initSize(&A, N, N);
  initSize(&B, N, N);
  initSize(&C_parallel, N, N);
  initSize(&C_seq, N, N);

  // Rank 0 generates the random data
  if (rank == 0)
  {
    setRandom(&A, 100);
    setRandom(&B, 100);
  }

  // Broadcast A and B
  MPI_Bcast(&(A.matrix[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(B.matrix[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ------------------ Parallel ------------------
  if (rank == 0)
    printf("Running parallel multiplyMatrix()...\n");

  double start = MPI_Wtime();
  multiplyMatrix(&A, &B, &C_parallel);
  double end = MPI_Wtime();

  if (rank == 0)
  {
    printf("multiplyMatrix() done in %.3f seconds\n", end - start);
    divisor("multiplyMatrix");
  }

  // ------------------ Sequential ------------------
  if (rank == 0)
  {
    printf("Running sequential Seq_multiplyMatrix()...\n");

    start = MPI_Wtime();
    Seq_multiplyMatrix(&A, &B, &C_seq);
    end = MPI_Wtime();

    printf("Seq_multiplyMatrix() done in %.3f seconds\n", end - start);
    divisor("Seq_multiplyMatrix");
  }

  MPI_Finalize();
  return 0;
}