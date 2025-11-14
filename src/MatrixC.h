#ifndef MATRIXC_H
#define MATRIXC_H

#define MAX_ROWS 2100
#define MAX_COLS 2100

#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  int rows;
  int cols;
  double matrix[MAX_ROWS][MAX_COLS];
} Matrix;

/* Error handling */
void error(const char *);
bool isValid(const Matrix *);

/* Accessors */
int getRows(const Matrix *);
int getCols(const Matrix *);

/* Constructors */
void init(Matrix *);
void initSize(Matrix *, const int, const int);
void initValue(Matrix *, const int, const int, const double);
void copyMatrix(Matrix *, const Matrix *);

/* Elementary row operations */
void rowScale(Matrix *, const int, const double);
void rowSwap(Matrix *, const int, const int);
void rowReplace(Matrix *, const int, const int, const double);

/* Core functions */
void display(const Matrix *);
void setValue(Matrix *, const int, const int, const double);
void setRandom(Matrix *, const int);
double getValue(Matrix *, const int, const int);
void displayDimensions(const Matrix *);
bool isSquare(const Matrix *);
void transpose(const Matrix *, Matrix *);
void writeToFile(const Matrix *, const char *);

/* Solving */
bool solve(const Matrix *, const double[], const int);
void multiplyMatrix(const Matrix *, const Matrix *, Matrix *);
void addMatrix(const Matrix *, const Matrix *, Matrix *);
void subtractMatrix(const Matrix *, const Matrix *, Matrix *);
void ref(Matrix *);
void rref(Matrix *);
void LU(Matrix *, Matrix *, Matrix *);

/* Non-OpenMPI Solving */
void Seq_multiplyMatrix(const Matrix *, const Matrix *, Matrix *);
void Seq_ref(Matrix *);
void Seq_rref(Matrix *);
void Seq_LU(const Matrix *, Matrix *, Matrix *);

#endif /* MATRIXC_H */