#include "Matrix.h"

bool Matrix::isValid() const {
    return (rows > 0 && cols > 0 && matrix != nullptr);
}

Matrix::Matrix() {
    matrix = nullptr;
    rows = cols = 0;
}

Matrix::Matrix(int r, int c) {
    if (r < 1 || c < 1) {
        throw string("Error: Invalid matrix dimension size.");
    }

    rows = r;
    cols = c;

    matrix = new double*[rows];

    for (int i = 0; i < rows; i++)
        matrix[i] = new double[cols];

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            matrix[i][j] = 0;
}

Matrix::Matrix(int r, int c, double defval) {
    if (r < 1 || c < 1) {
        throw string("Error: Invalid matrix dimension size.");
    }

    rows = r;
    cols = c;

    matrix = new double*[rows];

    for (int i = 0; i < rows; i++)
        matrix[i] = new double[cols];

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            matrix[i][j] = defval;
}

Matrix::Matrix(const Matrix& other) {
    rows = other.rows;
    cols = other.cols;

    if (other.isValid()) {
        matrix = new double*[rows];
        for (int i = 0; i < rows; i++) {
            matrix[i] = new double[cols];
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = other.matrix[i][j];
            }
        }
    } else {
        matrix = nullptr;
    }
}

Matrix::Matrix(string file) {

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
        throw string("Error opening file");

    while(getline(infile, line)) {
        // increment rows by one
        rows++;

        // count number of columns
        numWhitespace = 0;
        for (int i = 0; i < line.size(); i++)
            if (line[i] == ' ')
                numWhitespace++;

        // columns = # whitespace + 1
        cols = numWhitespace+1;
    }
    infile.close();

    // Confirm valid dimensions
    // Checking cols <= 1 because if there
    // are 0 whitespaces, cols still equals 1.
    if (rows == 0 || cols <= 1)
        throw string("Error: Invalid matrix dimensions");

    // Allocate memory
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++)
        matrix[i] = new double[cols];

    // Initialize array
    setRandom(10);

    // Reopen file and read in numbers
    infile.open(file);
    if (!infile)
        throw string("Error opening file");

    while(getline(infile, line)) {
        colIndex = 0;
        string currentNum = "";

        for (int i = 0; i < line.size(); i++) {
            if (isdigit(line[i]) || (line[i] == '-' && currentNum.empty()) || line[i] == '.')
                currentNum += line[i];
            else if (line[i] == ' ' || i == line.size()-1) {
                if(!currentNum.empty()) {
                    matrix[rowIndex][colIndex] = stod(currentNum);
                    currentNum = "";
                    colIndex++;
                }
            }
        }

        // Check for last num if line doesn't end with space
        if (!currentNum.empty())
            matrix[rowIndex][colIndex] = stod(currentNum);

        rowIndex++;
    }

    infile.close();

}

Matrix::~Matrix() {
    for (int i = 0; i < rows; i++)
        delete[] matrix[i];
    delete[] matrix;
}

void Matrix::setValue(int r, int c, double value) {
    if (r < 0 || c < 0 || r >= rows || c >= cols) {
        throw string("Error: Invalid matrix dimension size.");
    } 
    matrix[r][c] = value;
}

void Matrix::setRandom(int max) {
    if (isValid()) {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                matrix[i][j] = rand() % max;
    } else {
        throw string("Error: Matrix is not correctly initialized.");
    }
}

double Matrix::getValue(int r, int c) {
    if (r < 0 || c < 0 || r >= rows || c >= cols) {
        throw string("Error: Invalid matrix dimension size.");
    } 
    return matrix[r][c];
}

void Matrix::setDimensions(int r, int c) {
    if (r < 1 || c < 1) {
        throw string("Error: Invalid matrix dimension size.");
    }

    if (matrix != nullptr) {
        for (int i = 0; i < rows; i++)
            delete[] matrix[i];
        delete[] matrix;
    }

    rows = r;
    cols = c;

    matrix = new double*[rows];
    for (int i = 0; i < rows; i++)
        matrix[i] = new double[cols];
}

int Matrix::getRows() { return rows; }

int Matrix::getCols() { return cols; }

void Matrix::dimensions() {
    cout << "Dimension: [" << rows << "x" << cols << "]" << endl;
}

void Matrix::display() {
    if (isValid()) {
        for (int i = 0; i < rows; i++) {
            cout << "[";
            for (int j = 0; j < cols; j++) {
                cout << " " << matrix[i][j] << " ";
            }
            cout << "]" << endl;
        }
    } else {
        cout << "You are attempting to print out an invalid array." << endl;
    }
}

bool Matrix::solve(int solutions[], int size) {

    // Check valid size
    if (size != cols-1)
        throw string("Error: Invalid solution set size");

    // Check valid array
    if (!isValid()) 
        throw string("Error: Cannot check solution set on an empty array");

    double value = 0;

    for (int i = 0; i < rows; i++) {
        value = 0;
        for (int j = 0; j < cols-1; j++) {
            value += matrix[i][j] * solutions[j];
        }
        if (value != matrix[i][cols-1])
            return false;
    }
    return true;
}

double* Matrix::operator[](int index) {
    if (index < 0 || index >= rows) {
        throw string("Error: Invalid matrix index.");
    }
    return matrix[index];
}

Matrix& Matrix::operator=(const Matrix& other) {
    if (this == &other) return *this;

    for (int i = 0; i < rows; i++)
        delete[] matrix[i];
    delete[] matrix;

    rows = other.rows;
    cols = other.cols;
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = other.matrix[i][j];
        }
    }

    return *this;

}

Matrix Matrix::operator*(Matrix& other) {
    if (cols != other.rows || !isValid() || !other.isValid()) {
        throw string("Error: Cannot multiply matricies with invalid dimensions.");
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

Matrix Matrix::operator+(Matrix& other) {
    if (!isValid() || !other.isValid() || cols != other.cols || rows != other.rows) {
        throw string("Error: Cannot add matricies with invalid dimensions.");
    }

    int newRows = rows;
    int newCols = cols;
    Matrix newMatrix(newRows, newCols);

    for (int i = 0; i < newRows; i++)
        for(int j = 0; j < newCols; j++)
            newMatrix[i][j] = matrix[i][j] + other.matrix[i][j];

    return newMatrix;
}

Matrix Matrix::operator-(Matrix& other) {
    if (!isValid() || !other.isValid() || cols != other.cols || rows != other.rows) {
        throw string("Error: Cannot add matricies with invalid dimensions.");
    }

    int newRows = rows;
    int newCols = cols;
    Matrix newMatrix(newRows, newCols);

    for (int i = 0; i < newRows; i++)
        for(int j = 0; j < newCols; j++)
            newMatrix[i][j] = matrix[i][j] - other.matrix[i][j];

    return newMatrix;
}