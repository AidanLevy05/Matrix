#include "MatrixC.h"

/*
Name: error()
Parameters: const char*
Return: void
Description: Prints error message and exits program
*/
void error(const char *message)
{
  fprintf(stderr, "Error: %s\n", message);
  exit(EXIT_FAILURE);
}

/*
Name: isValid()
Parameters: const Matrix*
Return: bool
Description: Returns true if cols and rows are greater than 0
*/
bool isValid(const Matrix *m)
{
  return (m->rows > 0 && m->cols > 0);
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
void init(Matrix *m)
{
  m->rows = 0;
  m->cols = 0;
}

/*
Name: initSize()
Parameters: Matrix*, const int r, const int c
Return: void
Description: Constructor with set size parameters
*/
void initSize(Matrix *m, const int r, const int c)
{
  m->rows = r;
  m->cols = c;

  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < c; j++)
    {
      m->matrix[i][j] = 0;
    }
  }
}

/*
Name: initValue();
Parameters: Matrix*, const int r, const int c, const double defval
Return: void
Description: Constructor with set size parameters and default value
*/
void initValue(Matrix *m, const int r, const int c, const double defval)
{
  m->rows = r;
  m->cols = c;

  for (int i = 0; i < r; i++)
  {
    for (int j = 0; j < c; j++)
    {
      m->matrix[i][j] = defval;
    }
  }
}

/*
Name: copyMatrix()
Parameters: Matrix *dest, const Matrix *src
Return: void
Description: Copies src matrix into dest matrix
*/
void copyMatrix(Matrix *dest, const Matrix *src)
{
  if (dest == src)
    return;

  if (!isValid(dest) || !isValid(src))
  {
    error("Cannot copy invalid matrix");
  }

  dest->rows = src->rows;
  dest->cols = src->cols;

  for (int i = 0; i < src->rows; i++)
  {
    for (int j = 0; j < src->cols; j++)
    {
      dest->matrix[i][j] = src->matrix[i][j];
    }
  }
}

/*
Name: display()
Parameters: Matrix
Return: void
Description: Displays the matrix
*/
void display(const Matrix *m)
{
  int maxRows = 10, maxCols = 10;
  int showAllRows = m->rows <= maxRows;
  int showAllCols = m->cols <= maxCols;

  for (int i = 0; i < (showAllRows ? m->rows : maxRows); i++)
  {
    printf("[");
    for (int j = 0; j < (showAllCols ? m->cols : maxCols); j++)
    {
      printf(" %7.2f", m->matrix[i][j]);
    }
    if (!showAllCols)
      printf("   ...");
    printf(" ]\n");
  }
  if (!showAllRows)
  {
    int width = (showAllCols ? m->cols : maxCols) * 8 + (showAllCols ? 2 : 6);
    for (int i = 0; i < width; i++)
      printf(" ");
    printf("...\n");
  }
}

/*
Name: rowScale()
Parameters: Matrix*, int row, double scalar
Return: void
Description: Multiplies a row by a non-zero scalar
*/
void rowScale(Matrix *m, const int row, const double scalar)
{
  if (!isValid(m) || row >= m->rows || row < 0)
  {
    error("Invalid matrix size for row scaling");
  }

  for (int i = 0; i < m->cols; i++)
  {
    m->matrix[row][i] *= scalar;
  }
}

/*
Name: rowSwap()
Parameters: Matrix*, int row1, int row2
Return: void
Description: Swaps two valid rows of the matrix
*/
void rowSwap(Matrix *m, const int row1, const int row2)
{
  if (!isValid(m) || row1 >= m->rows || row2 >= m->rows || row1 < 0 || row2 < 0)
  {
    error("Invalid row access for row swap");
  }

  for (int i = 0; i < m->cols; i++)
  {
    double temp = m->matrix[row1][i];
    m->matrix[row1][i] = m->matrix[row2][i];
    m->matrix[row2][i] = temp;
  }
}

/*
Name: rowReplace()
Parameters: Matrix*, int targetRow, int sourceRow, double scalar
Return: void
Description: Replaces target row with itself plus scalar multiplied by source row
*/
void rowReplace(Matrix *m, const int targetRow, const int sourceRow, const double scalar)
{
  if (!isValid(m) || targetRow >= m->rows || sourceRow >= m->rows || targetRow < 0 || sourceRow < 0)
  {
    error("Invalid row access for row replace");
  }

  for (int i = 0; i < m->cols; i++)
  {
    m->matrix[targetRow][i] += scalar * m->matrix[sourceRow][i];
  }
}

/*
Name: setValue()
Parameters: Matrix*, int r, int c, double value
Return: void
Description: Sets an element of the row to a value
*/
void setValue(Matrix *m, const int r, const int c, const double value)
{
  if (!isValid(m) || r >= m->rows || c >= m->cols || r < 0 || c < 0)
  {
    error("Invalid row or column access for set value");
  }

  m->matrix[r][c] = value;
}

/*
Name: setRandom()
Parameters: Matrix*, const int maxRand
Return: void
Description: Sets matrix to random values
*/
void setRandom(Matrix *m, const int maxRand)
{
  if (!isValid(m))
  {
    error("Invalid matrix in set random");
  }

  for (int i = 0; i < m->rows; i++)
  {
    for (int j = 0; j < m->cols; j++)
    {
      m->matrix[i][j] = rand() % maxRand;
    }
  }
}

/*
Name: getValue()
Parameters: Matrix*, const int r, const int c
Return: double
Description: Returns matrix element at a slot
*/
double getValue(Matrix *m, const int r, const int c)
{
  if (!isValid(m) || r >= m->rows || c >= m->cols || r < 0 || c < 0)
  {
    error("Invalid row or column access in get value");
  }

  return m->matrix[r][c];
}

/*
Name: displayDimensions()
Parameters: Matrix*
Return: void
Description: Displays dimensions of the matrix
*/
void displayDimensions(const Matrix *m)
{
  if (!isValid(m))
  {
    error("Invalid matrix");
  }

  printf("[%d rows, %d cols]\n", m->rows, m->cols);
}

/*
Name: isSquare()
Parameters: const Matrix* A
Return: bool
Description: Returns true if the matrix is square
*/
bool isSquare(const Matrix *A)
{
  return A->rows == A->cols;
}

/*
Name: transpose()
Parameters: const Matrix *input, Matrix *output
Return: void
Description: Transposes the input matrix, stores result in output matrix
*/
void transpose(const Matrix *input, Matrix *output)
{

  if (!isValid(input))
  {
    error("Invalid matrix in transpose");
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Bcast((void *)input, sizeof(Matrix), MPI_BYTE, 0, MPI_COMM_WORLD);

  output->rows = input->cols;
  output->cols = input->rows;

  int rows_per_proc = input->rows / size;
  int remaining_rows = input->rows % size;
  int start = rank * rows_per_proc + (rank < remaining_rows ? rank : remaining_rows);
  int end = start + rows_per_proc + (rank < remaining_rows ? 1 : 0);

  for (int i = start; i < end; i++)
  {
    for (int j = 0; j < input->cols; j++)
    {
      output->matrix[j][i] = input->matrix[i][j];
    }
  }

  // back to all processes
  for (int proc = 0; proc < size; proc++)
  {
    int proc_start = proc * rows_per_proc + (proc < remaining_rows ? proc : remaining_rows);
    int proc_end = proc_start + rows_per_proc + (proc < remaining_rows ? 1 : 0);

    for (int row = 0; row < input->cols; row++)
    {
      for (int col = proc_start; col < proc_end; col++)
      {
        MPI_Bcast(&output->matrix[row][col], 1, MPI_DOUBLE, proc, MPI_COMM_WORLD);
      }
    }
  }
}

/*
Name: writeToFile()
Parameters: const Matrix*, const char* filename
Return: void
Description: Writes a matrix to a file
*/
void writeToFile(const Matrix *m, const char *filename)
{
  if (!isValid(m))
  {
    error("Invalid matrix");
  }

  FILE *file = fopen(filename, "w");
  if (file == NULL)
  {
    error("Unable to open file for writing");
  }

  fprintf(file, "Rows: %d, Cols: %d\n", m->rows, m->cols);
  for (int i = 0; i < 25; i++)
    fprintf(file, "-");
  fprintf(file, "\n");

  for (int i = 0; i < m->rows; i++)
  {
    for (int j = 0; j < m->cols; j++)
    {
      fprintf(file, "%lf ", m->matrix[i][j]);
    }
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
bool solve(const Matrix *m, const double solutions[], const int size)
{
  if (!isValid(m) || size < 1)
  {
    error("Cannot check solution set on empty matrix or empty solution set");
  }

  const double EPSILON = 1e-5;
  double value = 0;

  for (int i = 0; i < m->rows; i++)
  {
    value = 0;
    for (int j = 0; j < m->cols - 1; j++)
    {
      value += m->matrix[i][j] * solutions[j];
    }
    if (fabs(value - m->matrix[i][m->cols - 1] > EPSILON))
      return false;
  }
  return true;
}

/*
Name: multiplyMatrix()
Parameters: const Matrix *A, const Matrix *B, Matrix *C
Return: void
Description: Multiplies matrix A and B using OpenMPI and stores result in matrix C.
*/
void multiplyMatrix(const Matrix *A, const Matrix *B, Matrix *C)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int a_rows = A->rows;
  int a_cols = A->cols;
  int b_cols = B->cols;

  if (rank == 0)
  {
    if (!isValid(A) || !isValid(B))
      error("Invalid input matrices");
    if (a_cols != B->rows)
      error("Incompatible dimensions for multiplication");
  }

  // Broadcast matrix dimensions
  MPI_Bcast(&a_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Broadcast B (whole matrix)
  MPI_Bcast((void *)B->matrix, MAX_ROWS * MAX_COLS, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Calculate local rows
  int rows_per_proc = a_rows / size;
  int remainder = a_rows % size;
  int local_rows = (rank < remainder) ? rows_per_proc + 1 : rows_per_proc;
  printf("Rank %d received %d rows, computing...\n", rank, local_rows);

  // Allocate only as much as needed
  double local_A[local_rows][a_cols];
  double local_C[local_rows][b_cols];
  memset(local_C, 0, sizeof(local_C));

  // Prepare scatter metadata
  int sendcounts[size];
  int displs[size];
  int offset = 0;
  for (int i = 0; i < size; i++)
  {
    int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) * a_cols;
    sendcounts[i] = count;
    displs[i] = offset;
    offset += count;
  }

  // Scatter A rows to processes
  MPI_Scatterv(&(A->matrix[0][0]), sendcounts, displs, MPI_DOUBLE,
               &(local_A[0][0]), local_rows * a_cols, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

  // Compute local result
  for (int i = 0; i < local_rows; i++)
  {
    for (int j = 0; j < b_cols; j++)
    {
      for (int k = 0; k < a_cols; k++)
      {
        local_C[i][j] += local_A[i][k] * B->matrix[k][j];
      }
    }
  }

  // Prepare gather metadata
  int recvcounts[size];
  int recvdispls[size];
  offset = 0;
  for (int i = 0; i < size; i++)
  {
    int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) * b_cols;
    recvcounts[i] = count;
    recvdispls[i] = offset;
    offset += count;
  }

  // Gather results into final matrix C
  MPI_Gatherv(&(local_C[0][0]), local_rows * b_cols, MPI_DOUBLE,
              &(C->matrix[0][0]), recvcounts, recvdispls, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  // Set dimensions on root
  if (rank == 0)
  {
    C->rows = a_rows;
    C->cols = b_cols;
  }
}

/*
Name: ref()
Parameters: Matrix *m
Return: void
Description: Computes the Row Echelon Form (REF) of matrix m using OpenMPI.
*/
void ref(Matrix *m)
{
  if (!isValid(m))
    error("Invalid matrix");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rows = m->rows;
  int cols = m->cols;

  int rows_per_proc = rows / size;
  int remainder = rows % size;
  int local_rows = (rank < remainder) ? rows_per_proc + 1 : rows_per_proc;
  printf("Rank %d preparing to allocate %d rows × %d cols = %.2f MB\n", rank, local_rows, cols, (local_rows * cols * sizeof(double)) / (1024.0 * 1024.0));

  int start_row = (rank < remainder)
                      ? rank * (rows_per_proc + 1)
                      : rank * rows_per_proc + remainder;

  double local_matrix[local_rows][cols];
  printf("Rank %d received %d rows, computing...\n", rank, local_rows);

  // Scatter rows of m into local_matrix
  int sendcounts[size], displs[size];
  int offset = 0;
  for (int i = 0; i < size; i++)
  {
    int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) * cols;
    sendcounts[i] = count;
    displs[i] = offset;
    offset += count;
  }

  MPI_Scatterv(&(m->matrix[0][0]), sendcounts, displs, MPI_DOUBLE,
               &(local_matrix[0][0]), local_rows * cols, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

  // Temporary pivot row buffer
  double pivot_row[cols];

  for (int r = 0; r < rows; r++)
  {
    int owner = -1, local_r = -1;
    for (int i = 0; i < size; i++)
    {
      int start = (i < remainder) ? i * (rows_per_proc + 1) : i * rows_per_proc + remainder;
      int end = start + ((i < remainder) ? rows_per_proc + 1 : rows_per_proc);
      if (r >= start && r < end)
      {
        owner = i;
        local_r = r - start;
        break;
      }
    }

    // Broadcast pivot row
    if (rank == owner)
      memcpy(pivot_row, local_matrix[local_r], sizeof(double) * cols);

    MPI_Bcast(pivot_row, cols, MPI_DOUBLE, owner, MPI_COMM_WORLD);

    // Eliminate below pivot
    for (int i = 0; i < local_rows; i++)
    {
      int global_row = start_row + i;
      if (global_row <= r)
        continue;

      double factor = local_matrix[i][r] / pivot_row[r];
      for (int j = r; j < cols; j++)
      {
        local_matrix[i][j] -= factor * pivot_row[j];
      }
    }
  }

  // Gather results back to m
  MPI_Gatherv(&(local_matrix[0][0]), local_rows * cols, MPI_DOUBLE,
              &(m->matrix[0][0]), sendcounts, displs, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
}

/*
Name: rref()
Parameters: Matrix *m
Return: void
Description: Computes Reduced Row Echelon Form (RREF) using OpenMPI.
*/
void rref(Matrix *m)
{
  if (!isValid(m))
    error("Invalid matrix");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rows = m->rows;
  int cols = m->cols;

  int rows_per_proc = rows / size;
  int remainder = rows % size;
  int local_rows = (rank < remainder) ? rows_per_proc + 1 : rows_per_proc;

  int start_row = (rank < remainder)
                      ? rank * (rows_per_proc + 1)
                      : rank * rows_per_proc + remainder;

  double local_matrix[local_rows][cols];

  printf("Rank %d received %d rows, computing...\n", rank, local_rows);

  // Setup for scatter/gather
  int sendcounts[size], displs[size];
  int offset = 0;
  for (int i = 0; i < size; i++)
  {
    int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) * cols;
    sendcounts[i] = count;
    displs[i] = offset;
    offset += count;
  }

  // Scatter matrix
  MPI_Scatterv(&(m->matrix[0][0]), sendcounts, displs, MPI_DOUBLE,
               &(local_matrix[0][0]), local_rows * cols, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

  double pivot_row[cols];

  for (int r = 0; r < rows; r++)
  {
    int owner = -1, local_r = -1;
    for (int i = 0; i < size; i++)
    {
      int s = (i < remainder) ? i * (rows_per_proc + 1) : i * rows_per_proc + remainder;
      int e = s + ((i < remainder) ? rows_per_proc + 1 : rows_per_proc);
      if (r >= s && r < e)
      {
        owner = i;
        local_r = r - s;
        break;
      }
    }

    // Normalize pivot row on owner
    if (rank == owner)
    {
      double pivot = local_matrix[local_r][r];
      if (pivot != 0)
      {
        for (int j = r; j < cols; j++)
          local_matrix[local_r][j] /= pivot;
      }
      memcpy(pivot_row, local_matrix[local_r], sizeof(double) * cols);
    }

    // Broadcast normalized pivot row
    MPI_Bcast(pivot_row, cols, MPI_DOUBLE, owner, MPI_COMM_WORLD);

    // Eliminate all other rows (above and below)
    for (int i = 0; i < local_rows; i++)
    {
      int global_i = start_row + i;
      if (global_i == r)
        continue;

      double factor = local_matrix[i][r];
      for (int j = r; j < cols; j++)
        local_matrix[i][j] -= factor * pivot_row[j];
    }
  }

  // Gather result to root
  MPI_Gatherv(&(local_matrix[0][0]), local_rows * cols, MPI_DOUBLE,
              &(m->matrix[0][0]), sendcounts, displs, MPI_DOUBLE,
              0, MPI_COMM_WORLD);
}

/*
Name: addMatrix()
Parameters: const Matrix *A, const Matrix *B, Matrix *result
Return: Void
Description: Stores result of addition of A and B
*/
void addMatrix(const Matrix *A, const Matrix *B, Matrix *result)
{
  if (!isValid(A) || !isValid(B) || A->cols != B->cols || A->rows != B->rows)
  {
    error("Invalid matrix size");
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rows_per_proc = A->rows / size;
  int remaining_rows = A->rows % size;
  int start = rank * rows_per_proc + (rank < remaining_rows ? rank : remaining_rows);
  int end = start + rows_per_proc + (rank < remaining_rows ? 1 : 0);

  for (int i = start; i < end; i++)
  {
    for (int j = 0; j < A->cols; j++)
    {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, result->matrix, MAX_ROWS * MAX_COLS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

/*
Name: subtractMatrix()
Parameters: const Matrix *A, const Matrix *B, Matrix *result
Return: Void
Description: Stores result of subtraction of A and B
*/
void subtractMatrix(const Matrix *A, const Matrix *B, Matrix *result)
{
  if (!isValid(A) || !isValid(B) || A->cols != B->cols || A->rows != B->rows)
  {
    error("Invalid matrix size");
  }

  if (!isValid(A) || !isValid(B) || A->cols != B->cols || A->rows != B->rows)
  {
    error("Invalid matrix size");
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int rows_per_proc = A->rows / size;
  int remaining_rows = A->rows % size;
  int start = rank * rows_per_proc + (rank < remaining_rows ? rank : remaining_rows);
  int end = start + rows_per_proc + (rank < remaining_rows ? 1 : 0);

  for (int i = start; i < end; i++)
  {
    for (int j = 0; j < A->cols; j++)
    {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, result->matrix, MAX_ROWS * MAX_COLS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

/*
Name: LU()
Parameters: Matrix *A, Matrix *L, Matrix *U
Return: void
Description: Computes LU decomposition of matrix A using OpenMPI. A = L × U
*/
void LU(Matrix *A, Matrix *L, Matrix *U)
{
  if (!isValid(A))
    error("Invalid matrix");

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int n = A->rows;

  if (A->rows != A->cols)
    error("LU decomposition requires a square matrix");

  int rows_per_proc = n / size;
  int remainder = n % size;
  int local_rows = (rank < remainder) ? rows_per_proc + 1 : rows_per_proc;

  int start_row = (rank < remainder)
                      ? rank * (rows_per_proc + 1)
                      : rank * rows_per_proc + remainder;

  double local_matrix[local_rows][n];
  double pivot_row[n];

  // Setup scatter metadata
  int sendcounts[size], displs[size];
  int offset = 0;
  for (int i = 0; i < size; i++)
  {
    int count = ((i < remainder) ? rows_per_proc + 1 : rows_per_proc) * n;
    sendcounts[i] = count;
    displs[i] = offset;
    offset += count;
  }

  // Scatter A into local_matrix
  MPI_Scatterv(&(A->matrix[0][0]), sendcounts, displs, MPI_DOUBLE,
               &(local_matrix[0][0]), local_rows * n, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    initSize(L, n, n);
    initSize(U, n, n);
  }

  for (int k = 0; k < n; k++)
  {
    // Determine pivot owner and local row index
    int owner = -1, local_k = -1;
    for (int i = 0; i < size; i++)
    {
      int s = (i < remainder) ? i * (rows_per_proc + 1) : i * rows_per_proc + remainder;
      int e = s + ((i < remainder) ? rows_per_proc + 1 : rows_per_proc);
      if (k >= s && k < e)
      {
        owner = i;
        local_k = k - s;
        break;
      }
    }

    if (rank == owner)
      memcpy(pivot_row, local_matrix[local_k], sizeof(double) * n);

    MPI_Bcast(pivot_row, n, MPI_DOUBLE, owner, MPI_COMM_WORLD);

    for (int i = 0; i < local_rows; i++)
    {
      int global_i = start_row + i;
      if (global_i <= k)
        continue;

      double factor = local_matrix[i][k] / pivot_row[k];
      local_matrix[i][k] = factor;

      for (int j = k + 1; j < n; j++)
        local_matrix[i][j] -= factor * pivot_row[j];
    }
  }

  // Gather modified matrix back to A
  MPI_Gatherv(&(local_matrix[0][0]), local_rows * n, MPI_DOUBLE,
              &(A->matrix[0][0]), sendcounts, displs, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  // Extract L and U on rank 0
  if (rank == 0)
  {
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (i > j)
        {
          L->matrix[i][j] = A->matrix[i][j];
          U->matrix[i][j] = 0;
        }
        else if (i == j)
        {
          L->matrix[i][j] = 1;
          U->matrix[i][j] = A->matrix[i][j];
        }
        else
        {
          L->matrix[i][j] = 0;
          U->matrix[i][j] = A->matrix[i][j];
        }
      }
    }
  }
}

/*





Sequentual functions below





*/

/*
Name: Seq_multiplyMatrix()
Parameters: const Matrix *A, const Matrix *B, Matrix *C
Return: void
Description: Multiplies matrix A and matrix B sequentially and stores result in matrix C.
*/
void Seq_multiplyMatrix(const Matrix *A, const Matrix *B, Matrix *C)
{
  if (!isValid(A) || !isValid(B))
    error("Invalid input matrices");
  if (A->cols != B->rows)
    error("Incompatible dimensions for multiplication");

  for (int i = 0; i < A->rows; i++)
  {
    for (int j = 0; j < B->cols; j++)
    {
      C->matrix[i][j] = 0;
      for (int k = 0; k < A->cols; k++)
      {
        C->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
      }
    }
  }
}

/*
Name: Seq_ref()
Parameters: Matrix *m
Return: void
Description: Puts the matrix into row echelon form sequentially.
*/
void Seq_ref(Matrix *m)
{
  if (!isValid(m))
    error("Invalid matrix");

  int lead = 0;
  int rowCount = m->rows;
  int colCount = m->cols;

  for (int r = 0; r < rowCount; r++)
  {
    if (lead >= colCount)
      return;

    int i = r;
    while (m->matrix[i][lead] == 0)
    {
      i++;
      if (i == rowCount)
      {
        i = r;
        lead++;
        if (lead == colCount)
          return;
      }
    }

    if (i != r)
    {
      for (int j = 0; j < colCount; j++)
      {
        double temp = m->matrix[r][j];
        m->matrix[r][j] = m->matrix[i][j];
        m->matrix[i][j] = temp;
      }
    }

    double lv = m->matrix[r][lead];
    if (lv != 0)
    {
      for (int j = 0; j < colCount; j++)
        m->matrix[r][j] /= lv;
    }

    for (int i = r + 1; i < rowCount; i++)
    {
      double lv2 = m->matrix[i][lead];
      for (int j = 0; j < colCount; j++)
        m->matrix[i][j] -= lv2 * m->matrix[r][j];
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
void Seq_rref(Matrix *m)
{
  if (!isValid(m))
    error("Invalid matrix");

  Seq_ref(m);

  int rowCount = m->rows;
  int colCount = m->cols;

  for (int r = rowCount - 1; r >= 0; r--)
  {
    int leadCol = -1;
    for (int j = 0; j < colCount; j++)
    {
      if (fabs(m->matrix[r][j] - 1.0) < 1e-6)
      {
        leadCol = j;
        break;
      }
    }

    if (leadCol == -1)
      continue;

    for (int i = 0; i < r; i++)
    {
      double factor = m->matrix[i][leadCol];
      for (int j = 0; j < colCount; j++)
        m->matrix[i][j] -= factor * m->matrix[r][j];
    }
  }
}

/*
Name: Seq_LU()
Parameters: const Matrix *A, Matrix *L, Matrix *U
Return: void
Description: Performs LU decomposition sequentially.
*/
void Seq_LU(const Matrix *A, Matrix *L, Matrix *U)
{
  if (!isValid(A))
    error("Invalid matrix");

  int n = A->rows;

  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      L->matrix[i][j] = 0;
      U->matrix[i][j] = 0;
    }
  }

  for (int k = 0; k < n; ++k)
  {
    for (int j = k; j < n; ++j)
    {
      double sum = 0;
      for (int s = 0; s < k; ++s)
        sum += L->matrix[k][s] * U->matrix[s][j];
      U->matrix[k][j] = A->matrix[k][j] - sum;
    }

    L->matrix[k][k] = 1.0;

    for (int i = k + 1; i < n; ++i)
    {
      double sum = 0;
      for (int s = 0; s < k; ++s)
        sum += L->matrix[i][s] * U->matrix[s][k];
      L->matrix[i][k] = (A->matrix[i][k] - sum) / U->matrix[k][k];
    }
  }
}
