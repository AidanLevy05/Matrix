#ifndef MATRIX_H
#define MATRIX_H

#define MAX_COLS 100
#define MAX_ROWS 100

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
using namespace std;

class Matrix
{

private:
  int rows;
  int cols;
  double matrix[MAX_ROWS][MAX_COLS];

  bool isValid() const;
  double greatestCD(vector<double>);
  void swapRows(int, int);
  void multiplyRow(int, double);
  void replaceRow(int, int, double);
  bool infiniteSolutions(int);
  bool noSolutions(int);
  bool canSimplifyRow(int);

public:
  // Constructors and destructor
  Matrix();
  Matrix(int, int);
  Matrix(int, int, double);
  Matrix(string);
  Matrix(const Matrix &);
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
  bool solve(vector<double> &, int);
  void parametricForm();
  void parametricForm(vector<double>);
  void parametricVectorForm();
  vector<double> rref();
  bool isRREF();
  bool noSolutions();

  // Overloaded operators
  double *operator[](int);
  Matrix &operator=(const Matrix &);
  Matrix operator*(Matrix &);
  Matrix operator+(Matrix &);
  Matrix operator-(Matrix &);
};

#endif /* MATRIX_H */
