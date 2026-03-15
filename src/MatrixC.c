#include "MatrixC.h"

#include <limits.h>
#include <stdint.h>

#define MAT(m, r, c)                                                          \
  ((m)->matrix[(size_t)(r) * (size_t)((m)->cols) + (size_t)(c)])
#define ROW_PTR(m, r) ((m)->matrix + (size_t)(r) * (size_t)((m)->cols))

static size_t checkedProduct(const size_t a, const size_t b,
                             const char *message) {
  if (a != 0 && b > SIZE_MAX / a)
    error(message);
  return a * b;
}

static size_t checkedBytesForDoubles(const size_t count, const char *message) {
  return checkedProduct(count, sizeof(double), message);
}

static size_t checkedMatrixCount(const int rows, const int cols,
                                 const char *message) {
  if (rows <= 0 || cols <= 0)
    error(message);

  return checkedProduct((size_t)rows, (size_t)cols, message);
}

static size_t matrixCount(const Matrix *m) {
  return checkedProduct((size_t)m->rows, (size_t)m->cols,
                        "Matrix size is too large");
}

static int toMpiCount(const size_t count, const char *message) {
  if (count > INT_MAX)
    error(message);
  return (int)count;
}

static double *allocDoubles(const size_t count, const bool zero_fill,
                            const char *message) {
  if (count == 0)
    return NULL;

  const size_t bytes = checkedBytesForDoubles(count, message);
  double *buffer = zero_fill ? calloc(count, sizeof(double)) : malloc(bytes);
  if (buffer == NULL)
    error(message);

  return buffer;
}

static void ensureMatrixShape(Matrix *m, const int rows, const int cols,
                              const bool zero_fill,
                              const char *message) {
  if (m == NULL)
    error("Matrix pointer is NULL");

  const size_t count = checkedMatrixCount(rows, cols, message);
  const size_t bytes = checkedBytesForDoubles(count, message);

  if (m->matrix != NULL && m->rows == rows && m->cols == cols) {
    if (zero_fill)
      memset(m->matrix, 0, bytes);
  } else {
    free(m->matrix);
    m->matrix = allocDoubles(count, zero_fill, message);
  }

  m->rows = rows;
  m->cols = cols;
}

/*
Name: error()
Parameters: const char*
Return: void
Description: Prints error message and exits program
*/
void error(const char *message) {
  fprintf(stderr, "Error: %s\n", message);
  exit(EXIT_FAILURE);
}

/*
Name: isValid()
Parameters: const Matrix*
Return: bool
Description: Returns true if matrix has positive dimensions and storage
*/
bool isValid(const Matrix *m) {
  return (m != NULL && m->rows > 0 && m->cols > 0 && m->matrix != NULL);
}

/*
Name: getRows()
Parameters: const Matrix*
Return: int
Description: Returns the number of rows in a matrix
*/
int getRows(const Matrix *m) { return m->rows; }

/*
Name: getCols()
Parameters: const Matrix*
Return: int
Description: Returns the number of cols in a matrix
*/
int getCols(const Matrix *m) { return m->cols; }

/*
Name: init()
Parameters: Matrix*
Return: void
Description: Default constructor
*/
void init(Matrix *m) {
  if (m == NULL)
    error("Matrix pointer is NULL");

  m->rows = 0;
  m->cols = 0;
  m->matrix = NULL;
}

/*
Name: destroyMatrix()
Parameters: Matrix*
Return: void
Description: Releases matrix storage
*/
void destroyMatrix(Matrix *m) {
  if (m == NULL)
    return;

  free(m->matrix);
  m->matrix = NULL;
  m->rows = 0;
  m->cols = 0;
}

/*
Name: initSize()
Parameters: Matrix*, const int r, const int c
Return: void
Description: Constructor with set size parameters
*/
void initSize(Matrix *m, const int r, const int c) {
  ensureMatrixShape(m, r, c, true, "Unable to allocate matrix");
}

/*
Name: initValue()
Parameters: Matrix*, const int r, const int c, const double defval
Return: void
Description: Constructor with set size parameters and default value
*/
void initValue(Matrix *m, const int r, const int c, const double defval) {
  ensureMatrixShape(m, r, c, false, "Unable to allocate matrix");

  const size_t count = matrixCount(m);
  for (size_t i = 0; i < count; ++i)
    m->matrix[i] = defval;
}

/*
Name: copyMatrix()
Parameters: Matrix *dest, const Matrix *src
Return: void
Description: Copies src matrix into dest matrix
*/
void copyMatrix(Matrix *dest, const Matrix *src) {
  if (dest == src)
    return;

  if (!isValid(src))
    error("Cannot copy invalid matrix");

  initSize(dest, src->rows, src->cols);
  memcpy(dest->matrix, src->matrix,
         checkedBytesForDoubles(matrixCount(src), "Matrix copy is too large"));
}

/*
Name: display()
Parameters: Matrix
Return: void
Description: Displays the matrix
*/
void display(const Matrix *m) {
  if (!isValid(m))
    error("Invalid matrix");

  const int maxRows = 10;
  const int maxCols = 10;
  const int showRows = (m->rows <= maxRows) ? m->rows : maxRows;
  const int showCols = (m->cols <= maxCols) ? m->cols : maxCols;

  for (int i = 0; i < showRows; i++) {
    printf("[");
    for (int j = 0; j < showCols; j++)
      printf(" %7.2f", MAT(m, i, j));

    if (showCols < m->cols)
      printf("   ...");
    printf(" ]\n");
  }

  if (showRows < m->rows) {
    const int width = showCols * 8 + ((showCols < m->cols) ? 6 : 2);
    for (int i = 0; i < width; i++)
      printf(" ");
    printf("...\n");
  }
}

/*
Name: rowScale()
Parameters: Matrix*, int row, double scalar
Return: void
Description: Multiplies a row by a scalar
*/
void rowScale(Matrix *m, const int row, const double scalar) {
  if (!isValid(m) || row >= m->rows || row < 0)
    error("Invalid matrix size for row scaling");

  double *row_ptr = ROW_PTR(m, row);
  for (int i = 0; i < m->cols; i++)
    row_ptr[i] *= scalar;
}

/*
Name: rowSwap()
Parameters: Matrix*, int row1, int row2
Return: void
Description: Swaps two valid rows of the matrix
*/
void rowSwap(Matrix *m, const int row1, const int row2) {
  if (!isValid(m) || row1 >= m->rows || row2 >= m->rows || row1 < 0 ||
      row2 < 0) {
    error("Invalid row access for row swap");
  }

  double *row_a = ROW_PTR(m, row1);
  double *row_b = ROW_PTR(m, row2);
  for (int i = 0; i < m->cols; i++) {
    const double temp = row_a[i];
    row_a[i] = row_b[i];
    row_b[i] = temp;
  }
}

/*
Name: rowReplace()
Parameters: Matrix*, int targetRow, int sourceRow, double scalar
Return: void
Description: Replaces target row with itself plus scalar multiplied by source
row
*/
void rowReplace(Matrix *m, const int targetRow, const int sourceRow,
                const double scalar) {
  if (!isValid(m) || targetRow >= m->rows || sourceRow >= m->rows ||
      targetRow < 0 || sourceRow < 0) {
    error("Invalid row access for row replace");
  }

  double *target = ROW_PTR(m, targetRow);
  const double *source = ROW_PTR(m, sourceRow);
  for (int i = 0; i < m->cols; i++)
    target[i] += scalar * source[i];
}

/*
Name: setValue()
Parameters: Matrix*, int r, int c, double value
Return: void
Description: Sets an element of the row to a value
*/
void setValue(Matrix *m, const int r, const int c, const double value) {
  if (!isValid(m) || r >= m->rows || c >= m->cols || r < 0 || c < 0)
    error("Invalid row or column access for set value");

  MAT(m, r, c) = value;
}

/*
Name: setRandom()
Parameters: Matrix*, const int maxRand
Return: void
Description: Sets matrix to random values
*/
void setRandom(Matrix *m, const int maxRand) {
  if (!isValid(m))
    error("Invalid matrix in set random");
  if (maxRand <= 0)
    error("Random bound must be positive");

  const size_t count = matrixCount(m);
  for (size_t i = 0; i < count; ++i)
    m->matrix[i] = rand() % maxRand;
}

/*
Name: getValue()
Parameters: Matrix*, const int r, const int c
Return: double
Description: Returns matrix element at a slot
*/
double getValue(Matrix *m, const int r, const int c) {
  if (!isValid(m) || r >= m->rows || c >= m->cols || r < 0 || c < 0)
    error("Invalid row or column access in get value");

  return MAT(m, r, c);
}

/*
Name: displayDimensions()
Parameters: Matrix*
Return: void
Description: Displays dimensions of the matrix
*/
void displayDimensions(const Matrix *m) {
  if (!isValid(m))
    error("Invalid matrix");

  printf("[%d rows, %d cols]\n", m->rows, m->cols);
}

/*
Name: isSquare()
Parameters: const Matrix* A
Return: bool
Description: Returns true if the matrix is square
*/
bool isSquare(const Matrix *A) { return isValid(A) && A->rows == A->cols; }

/*
Name: transpose()
Parameters: const Matrix *input, Matrix *output
Return: void
Description: Transposes the input matrix, stores result in output matrix
*/
void transpose(const Matrix *input, Matrix *output) {
  if (output == NULL)
    error("Invalid output matrix in transpose");
  if (input == output)
    error("In-place transpose is not supported");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int input_rows = 0;
  int input_cols = 0;
  if (rank == 0) {
    if (!isValid(input))
      error("Invalid matrix in transpose");
    input_rows = input->rows;
    input_cols = input->cols;
  }

  MPI_Bcast(&input_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

  const size_t input_count =
      checkedMatrixCount(input_rows, input_cols, "Transpose input is too large");
  double *input_buffer =
      (rank == 0) ? input->matrix
                  : allocDoubles(input_count, false,
                                 "Unable to allocate transpose input buffer");

  MPI_Bcast(input_buffer,
            toMpiCount(input_count, "Transpose input is too large for MPI"),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);

  initSize(output, input_cols, input_rows);

  const int rows_per_proc = input_rows / size;
  const int remaining_rows = input_rows % size;
  const int start =
      rank * rows_per_proc + (rank < remaining_rows ? rank : remaining_rows);
  const int end = start + rows_per_proc + (rank < remaining_rows ? 1 : 0);

  for (int i = start; i < end; i++) {
    for (int j = 0; j < input_cols; j++) {
      MAT(output, j, i) = input_buffer[(size_t)i * (size_t)input_cols + j];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, output->matrix,
                toMpiCount(input_count,
                           "Transpose output is too large for MPI reduction"),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (rank != 0)
    free(input_buffer);
}

/*
Name: writeToFile()
Parameters: const Matrix*, const char* filename
Return: void
Description: Writes a matrix to a file
*/
void writeToFile(const Matrix *m, const char *filename) {
  if (!isValid(m))
    error("Invalid matrix");

  FILE *file = fopen(filename, "w");
  if (file == NULL)
    error("Unable to open file for writing");

  fprintf(file, "Rows: %d, Cols: %d\n", m->rows, m->cols);
  for (int i = 0; i < 25; i++)
    fprintf(file, "-");
  fprintf(file, "\n");

  for (int i = 0; i < m->rows; i++) {
    for (int j = 0; j < m->cols; j++)
      fprintf(file, "%lf ", MAT(m, i, j));
    fprintf(file, "\n");
  }

  fclose(file);
}

/*
Name: solve()
Parameters: const Matrix*, const double[], const int size
Return: bool
Description: Returns true if the solution set works with the matrix
*/
bool solve(const Matrix *m, const double solutions[], const int size) {
  if (!isValid(m) || solutions == NULL || size < 1)
    error("Cannot check solution set on empty matrix or empty solution set");
  if (size < m->cols - 1)
    error("Solution vector is too small for the matrix");

  const double EPSILON = 1e-5;

  for (int i = 0; i < m->rows; i++) {
    double value = 0.0;
    for (int j = 0; j < m->cols - 1; j++)
      value += MAT(m, i, j) * solutions[j];

    if (fabs(value - MAT(m, i, m->cols - 1)) > EPSILON)
      return false;
  }

  return true;
}

/*
Name: multiplyMatrix()
Parameters: const Matrix *A, const Matrix *B, Matrix *C
Return: void
Description: Multiplies matrix A and B using OpenMPI and stores result in matrix
C.
*/
void multiplyMatrix(const Matrix *A, const Matrix *B, Matrix *C) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int a_rows = 0;
  int a_cols = 0;
  int b_rows = 0;
  int b_cols = 0;

  if (rank == 0) {
    if (!isValid(A) || !isValid(B))
      error("Invalid input matrices");

    a_rows = A->rows;
    a_cols = A->cols;
    b_rows = B->rows;
    b_cols = B->cols;

    if (a_cols != b_rows)
      error("Incompatible dimensions for multiplication");
  }

  MPI_Bcast(&a_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (!isValid(B) || B->rows != b_rows || B->cols != b_cols)
    error("Matrix B must be initialized on all ranks before multiplication");

  const size_t b_count =
      checkedMatrixCount(b_rows, b_cols, "Matrix B is too large to broadcast");
  MPI_Bcast(B->matrix,
            toMpiCount(b_count, "Matrix B is too large for MPI broadcast"),
            MPI_DOUBLE, 0, MPI_COMM_WORLD);

  const int rows_per_proc = a_rows / size;
  const int remainder = a_rows % size;
  const int local_rows = (rank < remainder) ? rows_per_proc + 1 : rows_per_proc;
  printf("Rank %d received %d rows, computing...\n", rank, local_rows);

  const size_t local_a_count =
      checkedProduct((size_t)local_rows, (size_t)a_cols,
                     "Local multiplication input is too large");
  const size_t local_c_count =
      checkedProduct((size_t)local_rows, (size_t)b_cols,
                     "Local multiplication output is too large");
  double *local_A =
      allocDoubles(local_a_count, false,
                   "Unable to allocate local multiplication input");
  double *local_C =
      allocDoubles(local_c_count, true,
                   "Unable to allocate local multiplication output");

  int sendcounts[size];
  int displs[size];
  int offset = 0;
  for (int i = 0; i < size; i++) {
    const int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) *
                      a_cols;
    sendcounts[i] = count;
    displs[i] = offset;
    offset += count;
  }

  MPI_Scatterv(rank == 0 ? A->matrix : NULL, sendcounts, displs, MPI_DOUBLE,
               local_A,
               toMpiCount(local_a_count,
                          "Local multiplication input is too large for MPI"),
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int i = 0; i < local_rows; i++) {
    const double *a_row = local_A + (size_t)i * (size_t)a_cols;
    double *c_row = local_C + (size_t)i * (size_t)b_cols;
    for (int j = 0; j < b_cols; j++) {
      double sum = 0.0;
      for (int k = 0; k < a_cols; k++)
        sum += a_row[k] * B->matrix[(size_t)k * (size_t)b_cols + (size_t)j];
      c_row[j] = sum;
    }
  }

  int recvcounts[size];
  int recvdispls[size];
  offset = 0;
  for (int i = 0; i < size; i++) {
    const int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) *
                      b_cols;
    recvcounts[i] = count;
    recvdispls[i] = offset;
    offset += count;
  }

  if (rank == 0)
    initSize(C, a_rows, b_cols);

  MPI_Gatherv(local_C,
              toMpiCount(local_c_count,
                         "Local multiplication output is too large for MPI"),
              MPI_DOUBLE, rank == 0 ? C->matrix : NULL, recvcounts, recvdispls,
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(local_A);
  free(local_C);
}

/*
Name: ref()
Parameters: Matrix *m
Return: void
Description: Computes the Row Echelon Form (REF) of matrix m using OpenMPI.
*/
void ref(Matrix *m) {
  if (!isValid(m))
    error("Invalid matrix");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int rows = m->rows;
  const int cols = m->cols;
  const int rows_per_proc = rows / size;
  const int remainder = rows % size;
  const int local_rows = (rank < remainder) ? rows_per_proc + 1 : rows_per_proc;
  printf("Rank %d preparing to allocate %d rows × %d cols = %.2f MB\n", rank,
         local_rows, cols,
         (local_rows * cols * sizeof(double)) / (1024.0 * 1024.0));

  const int start_row = (rank < remainder) ? rank * (rows_per_proc + 1)
                                           : rank * rows_per_proc + remainder;
  const size_t local_count =
      checkedProduct((size_t)local_rows, (size_t)cols,
                     "Local REF buffer is too large");
  double *local_matrix =
      allocDoubles(local_count, false, "Unable to allocate local REF buffer");
  double *pivot_row =
      allocDoubles((size_t)cols, false, "Unable to allocate REF pivot row");
  printf("Rank %d received %d rows, computing...\n", rank, local_rows);

  int sendcounts[size];
  int displs[size];
  int offset = 0;
  for (int i = 0; i < size; i++) {
    const int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) *
                      cols;
    sendcounts[i] = count;
    displs[i] = offset;
    offset += count;
  }

  MPI_Scatterv(rank == 0 ? m->matrix : NULL, sendcounts, displs, MPI_DOUBLE,
               local_matrix,
               toMpiCount(local_count, "Local REF buffer is too large for MPI"),
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int r = 0; r < rows; r++) {
    int owner = -1;
    int local_r = -1;
    for (int i = 0; i < size; i++) {
      const int start = (i < remainder) ? i * (rows_per_proc + 1)
                                        : i * rows_per_proc + remainder;
      const int end =
          start + ((i < remainder) ? rows_per_proc + 1 : rows_per_proc);
      if (r >= start && r < end) {
        owner = i;
        local_r = r - start;
        break;
      }
    }

    if (rank == owner)
      memcpy(pivot_row, local_matrix + (size_t)local_r * (size_t)cols,
             checkedBytesForDoubles((size_t)cols, "REF pivot row is too large"));

    MPI_Bcast(pivot_row,
              toMpiCount((size_t)cols,
                         "REF pivot row is too large for MPI broadcast"),
              MPI_DOUBLE, owner, MPI_COMM_WORLD);

    if (fabs(pivot_row[r]) < 1e-12)
      error("REF encountered a zero pivot");

    for (int i = 0; i < local_rows; i++) {
      const int global_row = start_row + i;
      if (global_row <= r)
        continue;

      double *local_row = local_matrix + (size_t)i * (size_t)cols;
      const double factor = local_row[r] / pivot_row[r];
      for (int j = r; j < cols; j++)
        local_row[j] -= factor * pivot_row[j];
    }
  }

  MPI_Gatherv(local_matrix,
              toMpiCount(local_count, "Local REF buffer is too large for MPI"),
              MPI_DOUBLE, rank == 0 ? m->matrix : NULL, sendcounts, displs,
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(local_matrix);
  free(pivot_row);
}

/*
Name: rref()
Parameters: Matrix *m
Return: void
Description: Computes Reduced Row Echelon Form (RREF) using OpenMPI.
*/
void rref(Matrix *m) {
  if (!isValid(m))
    error("Invalid matrix");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int rows = m->rows;
  const int cols = m->cols;
  const int rows_per_proc = rows / size;
  const int remainder = rows % size;
  const int local_rows = (rank < remainder) ? rows_per_proc + 1 : rows_per_proc;
  const int start_row = (rank < remainder) ? rank * (rows_per_proc + 1)
                                           : rank * rows_per_proc + remainder;

  const size_t local_count =
      checkedProduct((size_t)local_rows, (size_t)cols,
                     "Local RREF buffer is too large");
  double *local_matrix =
      allocDoubles(local_count, false, "Unable to allocate local RREF buffer");
  double *pivot_row =
      allocDoubles((size_t)cols, false, "Unable to allocate RREF pivot row");

  printf("Rank %d received %d rows, computing...\n", rank, local_rows);

  int sendcounts[size];
  int displs[size];
  int offset = 0;
  for (int i = 0; i < size; i++) {
    const int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) *
                      cols;
    sendcounts[i] = count;
    displs[i] = offset;
    offset += count;
  }

  MPI_Scatterv(rank == 0 ? m->matrix : NULL, sendcounts, displs, MPI_DOUBLE,
               local_matrix,
               toMpiCount(local_count, "Local RREF buffer is too large for MPI"),
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  for (int r = 0; r < rows; r++) {
    int owner = -1;
    int local_r = -1;
    for (int i = 0; i < size; i++) {
      const int start = (i < remainder) ? i * (rows_per_proc + 1)
                                        : i * rows_per_proc + remainder;
      const int end =
          start + ((i < remainder) ? rows_per_proc + 1 : rows_per_proc);
      if (r >= start && r < end) {
        owner = i;
        local_r = r - start;
        break;
      }
    }

    if (rank == owner) {
      double *owner_row = local_matrix + (size_t)local_r * (size_t)cols;
      const double pivot = owner_row[r];
      if (pivot != 0.0) {
        for (int j = r; j < cols; j++)
          owner_row[j] /= pivot;
      }
      memcpy(pivot_row, owner_row,
             checkedBytesForDoubles((size_t)cols,
                                    "RREF pivot row is too large"));
    }

    MPI_Bcast(pivot_row,
              toMpiCount((size_t)cols,
                         "RREF pivot row is too large for MPI broadcast"),
              MPI_DOUBLE, owner, MPI_COMM_WORLD);

    for (int i = 0; i < local_rows; i++) {
      const int global_i = start_row + i;
      if (global_i == r)
        continue;

      double *local_row = local_matrix + (size_t)i * (size_t)cols;
      const double factor = local_row[r];
      for (int j = r; j < cols; j++)
        local_row[j] -= factor * pivot_row[j];
    }
  }

  MPI_Gatherv(local_matrix,
              toMpiCount(local_count, "Local RREF buffer is too large for MPI"),
              MPI_DOUBLE, rank == 0 ? m->matrix : NULL, sendcounts, displs,
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(local_matrix);
  free(pivot_row);
}

/*
Name: addMatrix()
Parameters: const Matrix *A, const Matrix *B, Matrix *result
Return: Void
Description: Stores result of addition of A and B
*/
void addMatrix(const Matrix *A, const Matrix *B, Matrix *result) {
  if (!isValid(A) || !isValid(B) || A->cols != B->cols || A->rows != B->rows)
    error("Invalid matrix size");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  initSize(result, A->rows, A->cols);

  const int rows_per_proc = A->rows / size;
  const int remaining_rows = A->rows % size;
  const int start =
      rank * rows_per_proc + (rank < remaining_rows ? rank : remaining_rows);
  const int end = start + rows_per_proc + (rank < remaining_rows ? 1 : 0);

  for (int i = start; i < end; i++) {
    for (int j = 0; j < A->cols; j++)
      MAT(result, i, j) = MAT(A, i, j) + MAT(B, i, j);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, result->matrix,
                toMpiCount(matrixCount(result),
                           "Matrix addition result is too large for MPI"),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

/*
Name: subtractMatrix()
Parameters: const Matrix *A, const Matrix *B, Matrix *result
Return: Void
Description: Stores result of subtraction of A and B
*/
void subtractMatrix(const Matrix *A, const Matrix *B, Matrix *result) {
  if (!isValid(A) || !isValid(B) || A->cols != B->cols || A->rows != B->rows)
    error("Invalid matrix size");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  initSize(result, A->rows, A->cols);

  const int rows_per_proc = A->rows / size;
  const int remaining_rows = A->rows % size;
  const int start =
      rank * rows_per_proc + (rank < remaining_rows ? rank : remaining_rows);
  const int end = start + rows_per_proc + (rank < remaining_rows ? 1 : 0);

  for (int i = start; i < end; i++) {
    for (int j = 0; j < A->cols; j++)
      MAT(result, i, j) = MAT(A, i, j) - MAT(B, i, j);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, result->matrix,
                toMpiCount(matrixCount(result),
                           "Matrix subtraction result is too large for MPI"),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

static void luDecompose(const Matrix *A, Matrix *L, Matrix *U) {
  if (!isValid(A) || !isSquare(A))
    error("LU decomposition requires a valid square matrix");

  const int n = A->rows;
  initSize(L, n, n);
  initSize(U, n, n);

  for (int i = 0; i < n; ++i)
    memcpy(ROW_PTR(U, i), ROW_PTR(A, i),
           checkedBytesForDoubles((size_t)n, "LU row copy is too large"));

  for (int k = 0; k < n; ++k) {
    double *pivot_row = ROW_PTR(U, k);
    const double pivot = pivot_row[k];

    if (fabs(pivot) < 1e-12)
      error("LU decomposition encountered a zero pivot");

    const double inv_pivot = 1.0 / pivot;
    MAT(L, k, k) = 1.0;

    for (int i = k + 1; i < n; ++i) {
      double *u_row = ROW_PTR(U, i);
      const double multiplier = u_row[k] * inv_pivot;
      MAT(L, i, k) = multiplier;
      u_row[k] = 0.0;

      for (int j = k + 1; j < n; ++j)
        u_row[j] -= multiplier * pivot_row[j];
    }
  }
}

/*
Name: LU
Parameters: Matrix* A, Matrix* L, Matrix* U
Return: void
Description: Performs LU decomposition
*/
void LU(Matrix *A, Matrix *L, Matrix *U) { luDecompose(A, L, U); }

/*
Name: Seq_multiplyMatrix()
Parameters: const Matrix *A, const Matrix *B, Matrix *C
Return: void
Description: Multiplies matrix A and matrix B sequentially and stores result in
matrix C.
*/
void Seq_multiplyMatrix(const Matrix *A, const Matrix *B, Matrix *C) {
  if (!isValid(A) || !isValid(B))
    error("Invalid input matrices");
  if (A->cols != B->rows)
    error("Incompatible dimensions for multiplication");

  initSize(C, A->rows, B->cols);

  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < B->cols; j++) {
      double sum = 0.0;
      for (int k = 0; k < A->cols; k++)
        sum += MAT(A, i, k) * MAT(B, k, j);
      MAT(C, i, j) = sum;
    }
  }
}

/*
Name: Seq_ref()
Parameters: Matrix *m
Return: void
Description: Puts the matrix into row echelon form sequentially.
*/
void Seq_ref(Matrix *m) {
  if (!isValid(m))
    error("Invalid matrix");

  int lead = 0;
  const int rowCount = m->rows;
  const int colCount = m->cols;

  for (int r = 0; r < rowCount; r++) {
    if (lead >= colCount)
      return;

    int i = r;
    while (MAT(m, i, lead) == 0.0) {
      i++;
      if (i == rowCount) {
        i = r;
        lead++;
        if (lead == colCount)
          return;
      }
    }

    if (i != r)
      rowSwap(m, r, i);

    const double pivot = MAT(m, r, lead);
    if (pivot != 0.0) {
      for (int j = 0; j < colCount; j++)
        MAT(m, r, j) /= pivot;
    }

    for (int row = r + 1; row < rowCount; row++) {
      const double factor = MAT(m, row, lead);
      for (int j = 0; j < colCount; j++)
        MAT(m, row, j) -= factor * MAT(m, r, j);
    }

    lead++;
  }
}

/*
Name: Seq_rref()
Parameters: Matrix *m
Return: void
Description: Converts matrix to reduced row echelon form without MPI.
*/
void Seq_rref(Matrix *m) {
  if (!isValid(m))
    error("Invalid matrix");

  Seq_ref(m);

  const int rowCount = m->rows;
  const int colCount = m->cols;

  for (int r = rowCount - 1; r >= 0; r--) {
    int leadCol = -1;
    for (int j = 0; j < colCount; j++) {
      if (fabs(MAT(m, r, j) - 1.0) < 1e-6) {
        leadCol = j;
        break;
      }
    }

    if (leadCol == -1)
      continue;

    for (int i = 0; i < r; i++) {
      const double factor = MAT(m, i, leadCol);
      for (int j = 0; j < colCount; j++)
        MAT(m, i, j) -= factor * MAT(m, r, j);
    }
  }
}

/*
Name: Seq_LU()
Parameters: const Matrix *A, Matrix *L, Matrix *U
Return: void
Description: Performs LU decomposition sequentially.
*/
void Seq_LU(const Matrix *A, Matrix *L, Matrix *U) { luDecompose(A, L, U); }
