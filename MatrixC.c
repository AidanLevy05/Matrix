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
  for (int i = 0; i < m->rows; i++)
  {
    printf("[");
    for (int j = 0; j < m->cols; j++)
    {
      printf(" %6.2f ", m->matrix[i][j]);
    }
    printf("]\n");
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
Name: isSquareMatrix()
Parameters: const Matrix* A
Return: bool
Description: Returns true if the matrix is square
*/
bool isSquare(const Matrix *A)
{
  return A->rows == A->cols;
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
Parameters: const Matrix* A, const Matrix* B, Matrix* result
Return: void
Description: Multiplies A * B and saves the result using MPI parallelism.
*/
void multiplyMatrix(const Matrix *A, const Matrix *B, Matrix *result)
{

  if (!isValid(A) || !isValid(B) || A->cols != B->rows)
  {
    error("Invalid matrix in matrix multiplication");
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
    for (int j = 0; j < B->cols; j++)
    {
      result->matrix[i][j] = 0;
      for (int k = 0; k < A->cols; k++)
      {
        result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, result->matrix, MAX_ROWS * MAX_COLS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
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
Name: ref()
Parameters: Matrix *
Return: void
Description: Puts the matrix into row echelon form
*/
void ref(Matrix *m)
{

  if (!isValid(m))
  {
    error("Invalid matrix");
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int lead = 0;
  int rowCount = m->rows;
  int colCount = m->cols;

  for (int r = 0; r < rowCount; r++)
  {
    if (lead >= colCount)
      return;

    // find pivot row for column lead
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

    // swap current row with pivot row
    if (i != r)
    {
      rowSwap(m, r, i);
    }

    // scale pivot row to make pivot 1
    double divisor = m->matrix[r][lead];
    for (int j = 0; j < colCount; j++)
    {
      m->matrix[r][j] /= divisor;
    }

    // OpenMPI starts here
    MPI_Bcast(m->matrix[r], colCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int rows_per_proc = (rowCount - (r + 1)) / size;
    int remaining_rows = (rowCount - (r + 1)) % size;

    // eliminate all entries below the pivot
    for (int i = r + 1; i < rowCount; i++)
    {
      if (m->matrix[i][lead] != 0)
      {
        double factor = m->matrix[i][lead];
        for (int j = 0; j < colCount; j++)
        {
          m->matrix[i][j] -= factor * m->matrix[r][j];
        }
      }
    }

    // gather everything back
    for (int proc = 0; proc < size; proc++)
    {
      int proc_start = r + 1 + proc * rows_per_proc + (proc < remaining_rows ? proc : remaining_rows);
      int proc_end = proc_start + rows_per_proc + (proc < remaining_rows ? 1 : 0);

      for (int row = proc_start; row < proc_end; row++)
      {
        MPI_Bcast(m->matrix[row], colCount, MPI_DOUBLE, proc, MPI_COMM_WORLD);
      }
    }

    lead++;
  }
}

/*
Name: rref();
Parameters: Matrix *
Return: void
Description: Puts the matrix into reduced row echelon form.
              This function also uses ref()
*/
void rref(Matrix *m)
{
  if (!isValid(m))
  {
    error("Invalid matrix");
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  ref(m);

  int rowCount = m->rows;
  int colCount = m->cols;

  for (int r = rowCount - 1; r >= 0; r--)
  {
    // skip if row is all zeroes
    if (m->matrix[r][r] == 0)
      continue;

    // Broadcast pivot row
    MPI_Bcast(m->matrix[r], colCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int rows_per_proc = r / size;
    int remaining_rows = r % size;
    int start = rank * rows_per_proc + (rank < remaining_rows ? rank : remaining_rows);
    int end = start + rows_per_proc + (rank < remaining_rows ? 1 : 0);

    for (int i = start; i < end; i++)
    {
      if (m->matrix[i][r] != 0)
      {
        double factor = m->matrix[i][r];
        for (int j = 0; j < colCount; j++)
        {
          m->matrix[i][j] -= factor * m->matrix[r][j];
        }
      }
    }

    // get other rows back
    // Broadcast updated rows back
    for (int proc = 0; proc < size; proc++)
    {
      int proc_start = proc * rows_per_proc + (proc < remaining_rows ? proc : remaining_rows);
      int proc_end = proc_start + rows_per_proc + (proc < remaining_rows ? 1 : 0);

      for (int row = proc_start; row < proc_end; row++)
      {
        MPI_Bcast(m->matrix[row], colCount, MPI_DOUBLE, proc, MPI_COMM_WORLD);
      }
    }
  }
}
