#include "../src/MatrixC.h"
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void divisor(const char *msg) { printf("\n---- %s ----\n\n", msg); }

// Global matrices
Matrix A_parallel = {0}, A_seq = {0};

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  srand(time(NULL) + rank);

  int N = 10, use_parallel = 1;

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

  initSize(&A_parallel, N, N);
  initSize(&A_seq, N, N);

  // Rank 0 generates the random matrix
  if (rank == 0) {
    setRandom(&A_parallel, 100);
    copyMatrix(&A_seq, &A_parallel);
  }

  // Broadcast A_parallel to all ranks
  MPI_Bcast(A_parallel.matrix, (int)element_count, MPI_DOUBLE, 0,
            MPI_COMM_WORLD);

  // ----------- Parallel ref() -----------
  if (rank == 0)
    printf("Running parallel ref()...\n");

  double start = MPI_Wtime();
  ref(&A_parallel);
  double end = MPI_Wtime();

  if (rank == 0) {
    printf("ref() done in %.3f seconds\n", end - start);
    divisor("ref");
  }

  // ----------- Sequential Seq_ref() -----------
  if (rank == 0 && use_parallel == 1) {
    printf("Running sequential Seq_ref()...\n");

    start = MPI_Wtime();
    Seq_ref(&A_seq);
    end = MPI_Wtime();

    printf("Seq_ref() done in %.3f seconds\n", end - start);
    divisor("Seq_ref");
  }

  destroyMatrix(&A_parallel);
  destroyMatrix(&A_seq);

  MPI_Finalize();
  return 0;
}
