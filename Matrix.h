#ifndef MATRIX_H
#define MATRIX_H 

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;

class Matrix {

private:
    double** matrix;
    int rows;
    int cols;

    bool isValid() const;
    double greatestCD(vector<double>);

public:
    // Constructors and destructor
    Matrix();
    Matrix(int, int);
    Matrix(int, int, double);
    Matrix(string);
    Matrix(const Matrix&);
    ~Matrix();

    // Basic functions
    void setValue(int, int, double);
    void setRandom(int);
    double getValue(int, int);
    void setDimensions(int, int);
    int getRows();
    int getCols();
    void dimensions();
    void display();

    // Solving
    // Functions with a * next to them will become private functions
    bool solve(vector<double>, int);
    vector<double> rref();
    void swapRows(int, int); // *
    void multiplyRow(int, double); // *
    void replaceRow(int, int, double); // *
    bool infiniteSolutions(int); // *
    bool noSolutions(int); // *
    bool canSimplifyRow(int); // *
    bool isRREF();

    // Overloaded operators
    double* operator[](int);
    Matrix& operator=(const Matrix&);
    Matrix operator*(Matrix&);
    Matrix operator+(Matrix&);
    Matrix operator-(Matrix&);

};

#endif /* MATRIX_H */
