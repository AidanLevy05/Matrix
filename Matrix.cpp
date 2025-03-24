#include "Matrix.h"

/*
Name: isValid()
Parameters: N/A
Return: bool
Description: Check if the matrix exists
*/
bool Matrix::isValid() const
{
  return (rows > 0 && cols > 0);
}

/*
Name: getGCD()
Parameters: long long a, long long b
Return: long long
Description: Returns the greatest common divisor
Note: Not associated with the Matrix Class
*/
long long getGCD(long long a, long long b)
{
  while (b != 0)
  {
    long long temp = b;
    b = a % b;
    a = temp;
  }
  return a;
}

/*
Name greatestCD()
Parameters: vector<double>
Return: double
Description: Returns greatest common divisor
            of doubles. Returns 0 if none found
*/
double greatestCD(vector<double> nums)
{

  if (nums.empty())
  {
    throw runtime_error("Error: Empty vector");
  }

  // 10 decimal point precision since
  // we are working with doubles
  const int precision = 10;

  double scale = pow(10, precision);
  vector<long long> scaled_nums(nums.size());
  transform(nums.begin(), nums.end(), scaled_nums.begin(), [&](double num)
            { return static_cast<long long>(round(num * scale)); });

  long long result = scaled_nums[0];
  for (long unsigned int i = 1; i < scaled_nums.size(); i++)
  {
    result = getGCD(result, scaled_nums[i]);
    if (result == 0)
      return 0;
  }

  return static_cast<double>(result) / scale;
}

/*
Name: infiniteSolutions()
Parameters: int row
Return: bool
Description: Returns true if a row contains all zeros.
*/
bool Matrix::infiniteSolutions(int row)
{

  // Valid matrix
  if (!isValid() || row >= rows)
  {
    throw runtime_error("Error: Invalid matrix operation");
  }

  if (matrix[row][cols - 1] != 0)
    return false;

  for (int i = 0; i < cols; i++)
  {
    if (matrix[row][i] != 0)
    {
      return false;
    }
  }
  return true;
}

/*
Name: noSolutions()
Parameters: int row
Return: bool
Description: Returns true if row contains leading zeroes
             followed by a non-zero number
*/
bool Matrix::noSolutions(int row)
{

  // Valid matrix
  if (!isValid() || row >= rows)
  {
    throw runtime_error("Error: Invalid matrix operation");
  }

  const double EPSILON = 1e-5;

  for (int i = 0; i < cols - 1; i++)
  {
    if (fabs(matrix[row][i]) > EPSILON)
    {
      return false;
    }
  }
  return fabs(matrix[row][cols - 1]) > EPSILON;
}

/*
Name: noSolutions()
Parameters: N/A
Return: bool
Description: Calls other noSolutions(int) and checks all rows.
*/
bool Matrix::noSolutions()
{

  // Valid matrix
  if (!isValid())
  {
    throw runtime_error("Error: Invalid matrix operation");
  }

  for (int i = 0; i < rows; i++)
    if (noSolutions(i))
      return true;

  return false;
}

/*
Name: canSimplifyRow()
Parameters: int row
Return: bool
Description: Returns true if the last element is non zero
            and there is only one non-zero number elsewhere
*/
bool Matrix::canSimplifyRow(int row)
{

  // Valid matrix
  if (!isValid() || row >= rows || cols <= 2 || rows <= 2)
  {
    throw runtime_error("Error: Invalid matrix operation");
  }

  if (matrix[row][cols - 1] == 0)
    return false;

  int nonZeroCount = 0;

  for (int i = 0; i < cols; i++)
  {
    if (matrix[row][i] != 0)
    {
      nonZeroCount++;
      if (nonZeroCount > 1)
        return false;
    }
  }

  return (nonZeroCount == 1);
}

/*
Name: Matrix()
Parameters: N/A
Return: N/A
Description: Default Constructor
*/
Matrix::Matrix()
{
}

/*
Name: Matrix()
Parameters: int r, int c
Return:  N/A
Description: Default Constructor with Matrix Size
*/
Matrix::Matrix(int r, int c)
{
  if (r < 1 || c < 1)
  {
    throw runtime_error("Error: Invalid matrix dimension size.");
  }

  rows = r;
  cols = c;

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      matrix[i][j] = 0;
}

/*
Name: Matrix()
Parameters: int r, int c, int defval
Return: N/A
Description: Default Constructor with Matrix Size and default value
*/
Matrix::Matrix(int r, int c, double defval)
{
  if (r < 1 || c < 1)
  {
    throw runtime_error("Error: Invalid matrix dimension size.");
  }

  rows = r;
  cols = c;

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      matrix[i][j] = defval;
}

/*
Name: Matrix()
Parameters: const Matrix&
Return: N/A
Description: Copy constructor
*/
Matrix::Matrix(const Matrix &other)
{
  rows = other.rows;
  cols = other.cols;

  if (other.isValid())
  {
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < cols; j++)
      {
        matrix[i][j] = other.matrix[i][j];
      }
    }
  }
  else
  {
    throw runtime_error("Error: Could not copy matrix.");
  }
}

/*
Name: Matrix()
Parameters: string file
Return: N/A
Description: Default constructor that reads in from file
*/
Matrix::Matrix(string file)
{

  // Before loading in the information,
  // open the file and count the number
  // of rows and columns.
  // Asume rows and cols = 0.
  rows = cols = 0;
  int numWhitespace = 0;
  int rowIndex = 0;
  int colIndex = 0;
  string line;

  ifstream infile(file);
  if (!infile)
    throw runtime_error("Error opening file");

  while (getline(infile, line))
  {
    // increment rows by one
    rows++;

    // count number of columns
    numWhitespace = 0;
    for (long unsigned int i = 0; i < line.size(); i++)
      if (line[i] == ' ')
        numWhitespace++;

    // columns = # whitespace + 1
    cols = numWhitespace + 1;
  }
  infile.close();

  // Confirm valid dimensions
  // Checking cols <= 1 because if there
  // are 0 whitespaces, cols still equals 1.
  if (rows == 0 || cols <= 1)
    throw runtime_error("Error: Invalid matrix dimensions");

  // Initialize array
  setRandom(10);

  // Reopen file and read in numbers
  infile.open(file);
  if (!infile)
    throw runtime_error("Error opening file");

  while (getline(infile, line))
  {
    colIndex = 0;
    string currentNum = "";

    for (long unsigned int i = 0; i < line.size(); i++)
    {
      if (isdigit(line[i]) || (line[i] == '-' && currentNum.empty()) || line[i] == '.')
        currentNum += line[i];
      else if (line[i] == ' ' || i == line.size() - 1)
      {
        if (!currentNum.empty())
        {
          matrix[rowIndex][colIndex] = stod(currentNum);
          currentNum = "";
          colIndex++;
        }
      }
      else
      {
        throw runtime_error("Error: Invalid character read from constructor");
      }
    }

    // Check for last num if line doesn't end with space
    if (!currentNum.empty())
      matrix[rowIndex][colIndex] = stod(currentNum);

    rowIndex++;
  }
  infile.close();
}

/*
Name: ~Matrix()
Parameters: N/A
Return: N/A
Description: Destructor
*/
Matrix::~Matrix()
{
}

/*
Name: setValue()
Parameters: int r, int c, value
Return: void
Description: Change one element of the matrix
*/
void Matrix::setValue(int r, int c, double value)
{
  if (r < 0 || c < 0 || r >= rows || c >= cols)
  {
    throw runtime_error("Error: Invalid matrix dimension size.");
  }
  matrix[r][c] = value;
}

/*
Name: setRandom()
Parameters: int max
Return: void
Description: Initializes array with random values given a max
*/
void Matrix::setRandom(int max)
{
  if (isValid())
  {
    for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++)
        matrix[i][j] = rand() % max;
  }
  else
  {
    throw runtime_error("Error: Matrix is not correctly initialized.");
  }
}

/*
Name: getValue
Parameters: int r, int c
Return: double
Description: Returns a given value from the matrix
*/
double Matrix::getValue(int r, int c)
{
  if (r < 0 || c < 0 || r >= rows || c >= cols)
  {
    throw runtime_error("Error: Invalid matrix dimension size.");
  }
  if (!isValid())
  {
    throw runtime_error("Error: Cannot return value from invalid matrix.");
  }
  return matrix[r][c];
}

/*
Name: setDimensions()
Parameters: int r, int c
Return: void
Description: Sets the parameters of the matrix
*/
void Matrix::setDimensions(int r, int c)
{
  if (r < 1 || c < 1)
  {
    throw runtime_error("Error: Invalid matrix dimension size.");
  }

  rows = r;
  cols = c;
}

/*
Name: getRows()
Parameters: N/A
Return: int
Description: Get the number of rows in the matrix
*/
int Matrix::getRows() { return rows; }

/*
Name: getCols()
Parameters: N/A
Return: int
Description: Get the number of cols in the matrix
*/
int Matrix::getCols() { return cols; }

/*
Name: dimensions()
Parameters: N/A
Return: void
Description: Prints the dimensions of the matrix
*/
void Matrix::dimensions()
{
  cout << "Dimension: [" << rows << "x" << cols << "]" << endl;
}

/*
Name: display()
Parameters: N/A
Return: void
Description: Prints all elements of the array neatly.
*/
void Matrix::display()
{
  if (isValid())
  {
    for (int i = 0; i < rows; i++)
    {
      cout << "[";
      for (int j = 0; j < cols; j++)
      {
        cout << " " << matrix[i][j] << " ";
      }
      cout << "]" << endl;
    }
  }
  else
  {
    throw runtime_error("You are attempting to print out an invalid matrix.");
  }
}

/*
Name: solve()
Parameters: vector<double>, int size
Return: bool
Description: Returns true if the solution set works; Otherwise, false.
*/
bool Matrix::solve(vector<double> &solutions, int size)
{

  // Check valid array
  if (!isValid() || solutions.size() != static_cast<long unsigned int>(size))
    throw runtime_error("Error: Cannot check solution set on an empty matrix");

  const double EPSILON = 1e-5;
  double value = 0;

  for (int i = 0; i < rows; i++)
  {
    value = 0;
    for (int j = 0; j < cols - 1; j++)
    {
      value += matrix[i][j] * solutions[j];
    }
    if (fabs(value - matrix[i][cols - 1]) > EPSILON)
      return false;
  }
  return true;
}

/*
Name: isRREF()
Parameters: N/A
Returns: bool
Description: Returns true if matrix is in RREF form
*/
bool Matrix::isRREF()
{

  if (!isValid())
    return false;

  // Count the number of zeroes in the
  // coefficient matrix.
  // If this number is less than
  // (rows x (cols-1)) - (rows)
  // then it is not in RREF.
  int totalZeroes = 0;
  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols - 1; j++)
    {
      if (matrix[i][j] == 0)
        totalZeroes++;
    }
  }

  return (totalZeroes == (rows * (cols - 1) - rows));
}

/*
Name: parametricForm()
Parameters: N/A
Return: Void
Description: Solves the matrix using rref(), then
            prints out the variables in parametric form
*/
void Matrix::parametricForm()
{

  // Validity
  if (!isValid())
    throw runtime_error("Error: Cannot do parametric form");

  vector<double> solutions = rref();

  if (rows == cols - 1)
  {
    cout << endl
         << "Parametric Form: " << endl;
    for (long unsigned int i = 0; i < solutions.size(); i++)
    {
      cout << "x" << (i + 1) << ": " << solutions[i] << endl;
    }
    cout << endl;
  }
  else if (rows != cols - 1 || noSolutions())
  {
    throw runtime_error("WARNING: This code is not completed yet.");
  }
}

/*
Name: parametricForm();
Parameters: vector<double>
Return: Void
Description: Same as above except is it designed to be called from the
            rref() function. Since the one above calls the rref()
            functiion, calling it twice would be a waste
            of resources.
*/
void Matrix::parametricForm(vector<double> sols)
{

  // Validity
  if (!isValid())
    throw runtime_error("Error: Cannot find parametric form");

  if (rows == cols - 1)
  {
    cout << endl
         << "Parametric Form: " << endl;
    for (long unsigned int i = 0; i < sols.size(); i++)
    {
      cout << "x" << (i + 1) << ": " << sols[i] << endl;
    }
    cout << endl;
  }
  else if (rows != cols - 1 || noSolutions())
  {
    throw runtime_error("WARNING: This code is not completed yet.");
  }
}

/*
Name: parametricVectorForm()
Parameters: N/A
Return: Void
Description: Prints out the solutions in parametric vector form nicely.
*/
void Matrix::parametricVectorForm()
{

  // Validity
  if (!isValid())
    throw runtime_error("Error: Cannot find parametric form");

  // Save solutions as an n x 1 vector
  if (rows == cols - 1)
  {
    Matrix vector(rows, 1);
    for (int i = 0; i < rows; i++)
      vector.setValue(i, 0, matrix[i][cols - 1]);

    cout << endl
         << "Parametric Vector Form:" << endl;
    vector.display();
    cout << endl;
  }
  else if (rows != cols - 1 || noSolutions())
  {
    throw runtime_error("WARNING: This code is not completed yet.");
  }
}

/*
Name: rref()
Parameters: N/A
Return: vector<double>
Description: Performs RREF operation on matrix, returns solutions
*/
vector<double> Matrix::rref()
{
  if (!isValid())
  {
    throw runtime_error("Error: Invalid matrix.");
  }

  vector<double> solutions;
  int row = 0, col = 0;
  const double EPSILON = 1e-5;
  bool noSols = false;

  while (row < rows && col < cols - 1)
  {

    if (fabs(matrix[row][col]) < EPSILON)
    {
      bool swapped = false;
      for (int i = row + 1; i < rows; i++)
      {
        if (fabs(matrix[i][col]) > EPSILON)
        {
          swapRows(row, i);
          swapped = true;
          break;
        }
      }
      if (!swapped)
      {
        col++;
        continue;
      }
    }

    double pivot = matrix[row][col];
    if (fabs(pivot) > EPSILON)
    {
      multiplyRow(row, 1.0 / pivot);
    }

    for (int i = row + 1; i < rows; i++)
    {
      if (fabs(matrix[i][col]) > EPSILON)
      {
        double scaleFactor = matrix[i][col];
        replaceRow(i, row, -scaleFactor);
      }
    }

    row++;
    col++;
  }

  for (int i = rows - 1; i >= 0; i--)
  {
    int pivotCol = -1;

    for (int j = 0; j < cols - 1; j++)
    {
      if (fabs(matrix[i][j]) > EPSILON)
      {
        pivotCol = j;
        break;
      }
    }

    if (pivotCol == -1)
      continue;

    double pivotValue = matrix[i][pivotCol];
    if (fabs(pivotValue - 1.0) > EPSILON)
    {
      multiplyRow(i, 1.0 / pivotValue);
    }

    for (int k = i - 1; k >= 0; k--)
    {
      if (fabs(matrix[k][pivotCol]) > EPSILON)
      {
        double scaleFactor = matrix[k][pivotCol];
        replaceRow(k, i, -1 * scaleFactor);
      }
    }
  }

  bool uniqueSolution = true;

  for (int i = 0; i < rows; i++)
  {
    bool isZeroRow = true;
    for (int j = 0; j < cols - 1; j++)
    {
      if (fabs(matrix[i][j]) > EPSILON)
      {
        isZeroRow = false;
        break;
      }
    }
    if (isZeroRow && fabs(matrix[i][cols - 1]) > EPSILON)
    {
      uniqueSolution = false;
      break;
    }
  }

  // check if inconsistent
  for (int i = 0; i < rows; i++)
  {
    if (noSolutions(i))
    {
      noSols = true;
    }
  }

  if (uniqueSolution && !noSols)
  {
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < cols - 1; j++)
      {
        if (fabs(matrix[i][j] - 1) < EPSILON)
        {
          solutions.push_back(matrix[i][cols - 1]);
          break;
        }
      }
    }
  }
  else if (!uniqueSolution && !noSols)
  {
    // seperate the basic and free variables
    // by 0.0 in the middle.
    vector<int> pivotCols;
    vector<int> freeCols;

    // find basic cols
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < cols - 1; j++)
      {
        if (fabs(matrix[i][j] - 1) < EPSILON)
        {
          pivotCols.push_back(j);
          break;
        }
      }
    }

    // find free cols
    for (int i = 0; i < cols - 1; i++)
    {
      if (find(pivotCols.begin(), pivotCols.end(), i) == pivotCols.end())
      {
        freeCols.push_back(i);
      }
    }

    // store solution for basic vars
    for (int i = 0; i < rows; i++)
    {
      for (int j = 0; j < cols - 1; j++)
      {
        if (fabs(matrix[i][j] - 1) < EPSILON)
        {
          solutions.push_back(matrix[i][cols - 1]);
          break;
        }
      }
    }

    // insert separator
    solutions.push_back(0.0);

    // store parametric components
    for (int freeCol : freeCols)
    {
      for (int i = 0; i < rows; i++)
      {
        if (fabs(matrix[i][freeCol]) > EPSILON)
        {
          solutions.push_back(-1 * matrix[i][freeCol]);
        }
      }
    }

    parametricForm(solutions);
  }

  if (noSols)
  {
    throw runtime_error("Warning: The matrix has no solutions.");
  }

  return solutions;
}

/*
Name: swapRows()
Parameters: int row1, int row2
Return: void
Description: Swaps two rows in the array
*/
void Matrix::swapRows(int row1, int row2)
{

  // Check for validity
  if (!isValid() || row1 >= rows || row2 >= rows)
  {
    throw runtime_error("Error: Invalid matrix index.");
  }

  // Don't waste resources swapping
  // the row with itself
  if (row1 == row2)
    return;

  // swap rows
  for (int i = 0; i < cols; i++)
  {
    double temp = matrix[row1][i];
    matrix[row1][i] = matrix[row2][i];
    matrix[row2][i] = temp;
  }
}

/*
Name: multiplyRow()
Parameters: int row, double value
Return: void
Description: Multiplies an entire row by a value
*/
void Matrix::multiplyRow(int row, double value)
{

  // Check validity
  if (!isValid() || row >= rows || value == 0)
    throw runtime_error("Error: Invalid matrix index.");

  for (int i = 0; i < cols; i++)
    matrix[row][i] *= value;
}

/*
Name: replaceRow()
Parameters: int replaceRow, int getRow, int value
Return: Void
Description: replaceRow = replaceRow + (value * getRow)
*/
void Matrix::replaceRow(int replaceRow, int getRow, double value)
{

  // Check validity
  if (!isValid() || replaceRow >= rows || getRow >= rows || value == 0)
    throw runtime_error("Error: Invalid matrix index.");

  if (replaceRow == getRow)
    return;

  for (int i = 0; i < cols; i++)
    matrix[replaceRow][i] += value * matrix[getRow][i];
}

/*
Name: operator[]
Parameters: int index
Return: double*
Description: Returns double* to access a given element
*/
double *Matrix::operator[](int index)
{
  if (index < 0 || index >= rows)
  {
    throw runtime_error("Error: Invalid matrix index.");
  }
  return matrix[index];
}

/*
Name: operator=
Parameters: const Matrix& other
Return: Matrix&
Description: Overloaded = operator
*/
Matrix &Matrix::operator=(const Matrix &other)
{
  if (this == &other)
    return *this;

  rows = other.rows;
  cols = other.cols;

  for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols; j++)
    {
      matrix[i][j] = other.matrix[i][j];
    }
  }

  return *this;
}

/*
Name: operator*
Parameters: Matrix& other
Return: Matrix
Description: Returns multiplied matrix
*/
Matrix Matrix::operator*(Matrix &other)
{
  if (cols != other.rows || !isValid() || !other.isValid())
  {
    throw runtime_error("Error: Cannot multiply matricies with invalid dimensions.");
  }

  // If you are multiplying matricies A x B = C
  // and A is size m x n
  // and B is size n x p,
  // Then C will be size m x p.
  int newRows = rows;
  int newCols = other.cols;
  Matrix newMatrix(newRows, newCols);

  for (int i = 0; i < newRows; i++)
    for (int j = 0; j < newCols; j++)
      for (int k = 0; k < cols; k++)
        newMatrix[i][j] += matrix[i][k] * other.matrix[k][j];

  return newMatrix;
}

/*
Name: operator+
Parameters: Matrix& other
Return: Matrix
Description: Returns sum of two matricies
*/
Matrix Matrix::operator+(Matrix &other)
{
  if (!isValid() || !other.isValid() || cols != other.cols || rows != other.rows)
  {
    throw runtime_error("Error: Cannot add matricies with invalid dimensions.");
  }

  int newRows = rows;
  int newCols = cols;
  Matrix newMatrix(newRows, newCols);

  for (int i = 0; i < newRows; i++)
    for (int j = 0; j < newCols; j++)
      newMatrix[i][j] = matrix[i][j] + other.matrix[i][j];

  return newMatrix;
}

/*
Name: operator-
Parameters: Matrix& other
Return: Matrix
Description: Returns difference of two matricies
*/
Matrix Matrix::operator-(Matrix &other)
{
  if (!isValid() || !other.isValid() || cols != other.cols || rows != other.rows)
  {
    throw runtime_error("Error: Cannot add matricies with invalid dimensions.");
  }

  int newRows = rows;
  int newCols = cols;
  Matrix newMatrix(newRows, newCols);

  for (int i = 0; i < newRows; i++)
    for (int j = 0; j < newCols; j++)
      newMatrix[i][j] = matrix[i][j] - other.matrix[i][j];

  return newMatrix;
}
