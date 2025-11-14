#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "MatrixC.h"

void divisor(const char *msg);

Matrix A_original, L_parallel, U_parallel;
Matrix L_seq, U_seq;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  srand(time(NULL) + rank);

  int N = 10, use_sequential = 0;

  if (rank == 0)
  {
    printf("Enter matrix size (N for NxN): ");
    scanf("%d", &N);

    printf("Run sequential also? (1 = yes, 0 = no): ");
    scanf("%d", &use_sequential);
  }

  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&use_sequential, 1, MPI_INT, 0, MPI_COMM_WORLD);

  initSize(&A_original, N, N);
  initSize(&L_parallel, N, N);
  initSize(&U_parallel, N, N);

  if (rank == 0)
  {
    initSize(&L_seq, N, N);
    initSize(&U_seq, N, N);
    setRandom(&A_original, 100);
  }

  MPI_Bcast(&(A_original.matrix[0][0]), N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ------------------ Parallel LU ------------------
  if (rank == 0)
    printf("Running parallel LU()...\n");

  double start = MPI_Wtime();
  LU(&A_original, &L_parallel, &U_parallel);
  double end = MPI_Wtime();

  if (rank == 0)
  {
    printf("Parallel LU() done in %.3f seconds\n", end - start);
    divisor("Parallel LU");

    if (use_sequential == 1)
    {
      Matrix A_copy;
      initSize(&A_copy, N, N);
      copyMatrix(&A_copy, &A_original);

      printf("Running sequential Seq_LU()...\n");

      start = MPI_Wtime();
      Seq_LU(&A_copy, &L_seq, &U_seq);
      end = MPI_Wtime();

      printf("Sequential Seq_LU() done in %.3f seconds\n", end - start);
      divisor("Sequential LU");
    }
  }

  MPI_Finalize();
  return 0;
}

void divisor(const char *msg)
{
  printf("\n---- %s ----\n\n", msg);
}
