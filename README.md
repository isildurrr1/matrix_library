# matrix_library
> Library for operations on matrices, written in C language
## Library Functions:

### 1. Creating matrices (create_matrix)
```c
int create_matrix(int rows, int columns, matrix_t *result);
```
### 2. Cleaning of matrices (remove_matrix)
```c
void remove_matrix(matrix_t *A);
```
### 3. Matrix comparison (eq_matrix)
```c
int eq_matrix(matrix_t *A, matrix_t *B);
```
### 4. Adding (sum_matrix)
```c
int sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
```
### 5. Subtracting matrices (sub_matrix)
```c
int sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
```
### 6. Matrix multiplication by scalar (mult_number)
```c
int mult_number(matrix_t *A, double number, matrix_t *result);
```
### 7. Multiplication of two matrices (mult_matrix)
```c
int mult_number(matrix_t *A, double number, matrix_t *result);
```
### 8. Matrix transpose (transpose)
```c
int transpose(matrix_t *A, matrix_t *result);
```
### 9. Matrix of algebraic complements (calc_complements)
```c
int calc_complements(matrix_t *A, matrix_t *result);
```
### 10. Matrix determinant (determinant)
```c
int determinant(matrix_t *A, double *result);
```
### 11. Inverse of the matrix (inverse_matrix)
```c
int inverse_matrix(matrix_t *A, matrix_t *result);
```
### 12. Matrix minor (minor_matrix)
```c
matrix_t minor_matrix(matrix_t *A, int row, int column);
```