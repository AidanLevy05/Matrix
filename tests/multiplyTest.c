#include "../src/MatrixC.h"
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void divisor(const char *msg) { printf("\n---- %s ----\n\n", msg); }

Matrix A = {0}, B = {0}, C_parallel = {0}, C_seq = {0};

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  srand(time(NULL) + rank);

  int N = 100, use_parallel = 0;

  if (rank == 0) {
    printf("Enter matrix size: ");
    scanf("%d", &N);

    printf("Run sequential too? (1 = yes, 0 = no): ");
    scanf("%d", &use_parallel);
  }

  // Broadcast settings
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&use_parallel, 1, MPI_INT, 0, MPI_COMM_WORLD);

  const size_t element_count = (size_t)N * (size_t)N;
  if (element_count > INT_MAX) {
    if (rank == 0)
      fprintf(stderr, "Error: matrix size is too large for MPI counts\n");
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  initSize(&A, N, N);
  initSize(&B, N, N);
  initSize(&C_parallel, N, N);
  initSize(&C_seq, N, N);

  // Rank 0 generates the random data
  if (rank == 0) {
    setRandom(&A, 100);
    setRandom(&B, 100);
  }

  // Broadcast A and B
  MPI_Bcast(A.matrix, (int)element_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(B.matrix, (int)element_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ------------------ Parallel ------------------
  if (rank == 0)
    printf("Running parallel multiplyMatrix()...\n");

  double start = MPI_Wtime();
  multiplyMatrix(&A, &B, &C_parallel);
  double end = MPI_Wtime();

  if (rank == 0) {
    printf("multiplyMatrix() done in %.3f seconds\n", end - start);
    divisor("multiplyMatrix");
  }

  // ------------------ Sequential ------------------
  if (rank == 0 && use_parallel == 1) {
    printf("Running sequential Seq_multiplyMatrix()...\n");

    start = MPI_Wtime();
    Seq_multiplyMatrix(&A, &B, &C_seq);
    end = MPI_Wtime();

    printf("Seq_multiplyMatrix() done in %.3f seconds\n", end - start);
    divisor("Seq_multiplyMatrix");
  }

  destroyMatrix(&A);
  destroyMatrix(&B);
  destroyMatrix(&C_parallel);
  destroyMatrix(&C_seq);

  MPI_Finalize();
  return 0;
}
