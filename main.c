#include "MatrixC.h"

void divisor(const char *);

int main(int argc, char **argv)
{
  // Initialization
  srand(time(0));
  const int RANDOM_MAX = 100;
  MPI_Init(&argc, &argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  Matrix A, B, C, D;

  initValue(&A, 2, 3, 1);
  initValue(&B, 3, 2, 2);
  multiplyMatrix(&A, &B, &C);

  if (rank == 0)
  {
    divisor("A");
    display(&A);
    divisor("B");
    display(&B);
    divisor("C = A * B");
    display(&C);

    divisor("REF on C");
    ref(&C);
    display(&C);
    divisor("RREF on C");
    rref(&C);
    display(&C);

    divisor("Randomize D");
  }

  initSize(&D, 3, 3);
  setRandom(&D, RANDOM_MAX);

  // Synchronize processes before printing (only root process will print)
  MPI_Barrier(MPI_COMM_WORLD);

  if (rank == 0)
  {
    display(&D);
    divisor("RREF D");
    rref(&D);
    display(&D);
  }

  MPI_Finalize();
  return 0;
}

void divisor(const char *message)
{
  printf("\n");
  printf("%s ", message);
  for (int i = 0; i < 25; i++)
    printf("-");
  printf("\n\n");
}
