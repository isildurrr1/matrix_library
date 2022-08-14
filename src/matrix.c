#include "matrix.h"

// Получаем матрицу с вычеркнутой строчкой/столбцом (минор)
matrix_t minor_matrix(matrix_t *A, int row, int column) {
  matrix_t minor;
  create_matrix(A->rows - 1, A->columns - 1, &minor);
  for (int r = 0, m_r = 0; r < A->rows; ++r) {
    if (r != row) {
      for (int c = 0, m_c = 0; c < A->columns; ++c) {
        if (c != column) {
          minor.matrix[m_r][m_c++] = A->matrix[r][c];
        }
      }
      ++m_r;
    }
  }
  return minor;
}

/*Создание матриц*/
int create_matrix(int rows, int columns, matrix_t *result) {
  int error = OK;
  memset(result, 0, sizeof(matrix_t));
  if (rows <= 0 || columns <= 0) {
    error = INCORRECT_MATRIX;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows, sizeof(double *));
    for (int r = 0; r < rows; r++) {
      result->matrix[r] = (double *)calloc(columns, sizeof(double *));
    }
  }
  return error;
}

// Очистка матрицы
void remove_matrix(matrix_t *A) {
  if (A) {
    if (A->matrix && A->rows && A->columns) {
      for (int r = 0; r < A->rows; r++)
        if (A->matrix[r]) free(A->matrix[r]);
      free(A->matrix);
      A->rows = 0;
      A->columns = 0;
    }
  }
}

// Сравнение матриц
int eq_matrix(matrix_t *A, matrix_t *B) {
  int result = FAILURE;
  if (A->rows == B->rows && A->columns == B->columns) {
    result = SUCCESS;
    for (int r = 0; r < A->rows; r++) {
      for (int c = 0; c < A->columns; c++) {
        if (fabs(A->matrix[r][c] - B->matrix[r][c]) >= EPS) result = FAILURE;
      }
    }
  }
  return result;
}

// Сложение матриц
int sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  int Ar = A->rows, Ac = A->columns;
  int Br = B->rows, Bc = B->columns;
  if (Ar <= 0 || Ac <= 0 || Br <= 0 || Bc <= 0 || !result) {
    error = INCORRECT_MATRIX;
  } else {
    if (Ar == Br && Ac == Bc) {
      create_matrix(Ar, Ac, result);
      for (int r = 0; r < Ar; r++) {
        for (int c = 0; c < Ac; c++) {
          result->matrix[r][c] = A->matrix[r][c] + B->matrix[r][c];
        }
      }
    } else {
      error = CALC_ERROR;
    }
  }
  return error;
}

// Вычитание матриц
int sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (A->rows <= 0 || A->columns <= 0 || B->rows <= 0 || B->columns <= 0 ||
      !result) {
    error = INCORRECT_MATRIX;
  } else {
    if (A->rows == B->rows && A->columns == B->columns) {
      create_matrix(A->rows, A->columns, result);
      for (int r = 0; r < A->rows; r++) {
        for (int c = 0; c < A->columns; c++) {
          result->matrix[r][c] = A->matrix[r][c] - B->matrix[r][c];
        }
      }
    } else {
      error = CALC_ERROR;
    }
  }
  return error;
}

// Умножение матрицы на число
int mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = OK;
  if (A->rows <= 0 || A->columns <= 0) {
    error = INCORRECT_MATRIX;
  } else {
    if (result) {
      create_matrix(A->rows, A->columns, result);
      for (int r = 0; r < A->rows; r++) {
        for (int c = 0; c < A->columns; c++) {
          result->matrix[r][c] = A->matrix[r][c] * number;
        }
      }
    }
  }
  return error;
}

// Умножение матриц
int mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (A->rows <= 0 || A->columns <= 0 || B->rows <= 0 || B->columns <= 0) {
    error = INCORRECT_MATRIX;
  } else {
    if (A->columns == B->rows) {
      create_matrix(A->rows, B->columns, result);
      for (int r = 0; r < result->rows; r++) {
        for (int c = 0; c < result->columns; c++) {
          for (int k = 0; k < A->columns; k++) {
            result->matrix[r][c] += A->matrix[r][k] * B->matrix[k][c];
          }
        }
      }
    } else {
      error = CALC_ERROR;
    }
  }
  return error;
}

// Транспонирование матрицы
int transpose(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (A->rows <= 0 || A->columns <= 0 || !result) {
    error = INCORRECT_MATRIX;
  } else {
    create_matrix(A->columns, A->rows, result);
    for (int r = 0; r < result->rows; r++) {
      for (int c = 0; c < result->columns; c++) {
        result->matrix[r][c] = A->matrix[c][r];
      }
    }
  }
  return error;
}

// Матрица алгебраических дополнений
int calc_complements(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (A->rows <= 0 || A->columns <= 0) {
    error = INCORRECT_MATRIX;
  } else {
    if (A->rows == A->columns) {
      create_matrix(A->rows, A->columns, result);
      if (A->rows == 1) {
        result->matrix[0][0] = 1.0;
      } else {
        for (int r = 0; r < result->rows; r++) {
          for (int c = 0; c < result->columns; c++) {
            double det = 0.0;
            matrix_t minor = minor_matrix(A, r, c);
            determinant(&minor, &det);
            result->matrix[r][c] = det * pow(-1, r + c);
            remove_matrix(&minor);
          }
        }
      }
    } else {
      error = CALC_ERROR;
    }
  }
  return error;
}

// Вычесляем определитель матрицы
int determinant(matrix_t *A, double *result) {
  double det;
  int error = OK;
  if (A->rows <= 0 || A->columns <= 0) {
    error = INCORRECT_MATRIX;
  } else {
    *result = 0.0;
    if (A->rows == A->columns) {
      if (A->rows == 1) {
        *result = A->matrix[0][0];
      } else if (A->rows == 2) {
        *result = A->matrix[0][0] * A->matrix[1][1] -
                  A->matrix[0][1] * A->matrix[1][0];
      } else {
        det = 0.0;
        for (int c = 0; c < A->columns; c++) {
          matrix_t minor = minor_matrix(A, 0, c);
          determinant(&minor, &det);
          *result += A->matrix[0][c] * pow(-1.0, c) * det;
          remove_matrix(&minor);
        }
      }
    } else {
      error = CALC_ERROR;
    }
  }
  return error;
}

// Вычисляем обратную матрицу
int inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = OK;
  double det = 0.0;
  error = determinant(A, &det);
  if (fabs(det) < EPS && error == OK) {
    error = CALC_ERROR;
  } else if (error == OK) {
    matrix_t m_alg_dop, m_trsp;
    calc_complements(A, &m_alg_dop);
    transpose(&m_alg_dop, &m_trsp);
    mult_number(&m_trsp, 1 / det, result);
    remove_matrix(&m_alg_dop);
    remove_matrix(&m_trsp);
  }
  return error;
}
