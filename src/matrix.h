#ifndef SRC_MATRIX_H_
#define SRC_MATRIX_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SUCCESS 1
#define FAILURE 0
#define OK 0
#define INCORRECT_MATRIX 1
#define CALC_ERROR 2
#define EPS 0.0000001

typedef struct matrix_struct {
  double** matrix;
  int rows;
  int columns;
} matrix_t;

int create_matrix(int rows, int columns, matrix_t* result);
void remove_matrix(matrix_t* A);
int eq_matrix(matrix_t* A, matrix_t* B);
int sum_matrix(matrix_t* A, matrix_t* B, matrix_t* result);
int sub_matrix(matrix_t* A, matrix_t* B, matrix_t* result);
int mult_number(matrix_t* A, double number, matrix_t* result);
int mult_matrix(matrix_t* A, matrix_t* B, matrix_t* result);
int transpose(matrix_t* A, matrix_t* result);
int calc_complements(matrix_t* A, matrix_t* result);
int determinant(matrix_t* A, double* result);
int inverse_matrix(matrix_t* A, matrix_t* result);
matrix_t minor_matrix(matrix_t* A, int row, int column);

#endif  // SRC_MATRIX_H_
