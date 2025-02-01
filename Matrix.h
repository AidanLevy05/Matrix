#ifndef MATRIX_H
#define MATRIX_H 

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
using namespace std;

class Matrix {

private:
    double** matrix;
    int rows;
    int cols;

    bool isValid() const;

public:
    // Constructors and destructors
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
    bool solve(int[], int);

    // Overloaded operators
    double* operator[](int);
    Matrix& operator=(const Matrix&);
    Matrix operator*(Matrix&);
    Matrix operator+(Matrix&);
    Matrix operator-(Matrix&);

};

#endif
