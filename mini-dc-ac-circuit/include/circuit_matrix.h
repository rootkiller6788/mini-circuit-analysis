#ifndef CIRCUIT_MATRIX_H
#define CIRCUIT_MATRIX_H
#include <stddef.h>
#include <complex.h>
#include "circuit_elements.h"

/* Linear Algebra for Circuit Analysis (L3: Mathematical Structures)
 * Provides matrix operations for nodal/mesh analysis, MNA formulation.
 * Reference: Golub & Van Loan "Matrix Computations" 4th ed.
 *            Horn & Johnson "Matrix Analysis"
 *
 * All matrices are stored in row-major order.
 * Complex matrices support AC steady-state analysis. */

/* Real matrix operations */

/* Allocate n x m matrix (row-major). Caller must free(). */
double* matrix_alloc(int rows, int cols);

/* Free matrix allocated by matrix_alloc(). */
void matrix_free(double *m);

/* Zero all elements: A[i][j] = 0 */
void matrix_zero(double *A, int rows, int cols);

/* Set to identity: A[i][i] = 1, A[i][j] = 0 for i != j.
 * Matrix must be square (rows == cols). */
void matrix_identity(double *A, int n);

/* Copy: B = A (rows x cols) */
void matrix_copy(const double *A, double *B, int rows, int cols);

/* C = A + B (rows x cols) */
void matrix_add(const double *A, const double *B, double *C, int rows, int cols);

/* C = A - B (rows x cols) */
void matrix_sub(const double *A, const double *B, double *C, int rows, int cols);

/* C = A * B  where A is m x k, B is k x n, C is m x n.
 * Standard O(m*n*k) triple-loop multiplication. */
void matrix_mul(const double *A, const double *B, double *C,
                int m, int k, int n);

/* y = A * x  matrix-vector multiply, A is m x n. */
void matrix_vec_mul(const double *A, const double *x, double *y, int m, int n);

/* B = A^T (rows x cols -> cols x rows) */
void matrix_transpose(const double *A, double *B, int rows, int cols);

/* Determinant of n x n matrix using LU decomposition.
 * O(n^3) time. Returns 0 on success, sets *det to determinant.
 * Returns -1 if singular. */
int matrix_determinant(const double *A, int n, double *det);

/* LU Decomposition (Doolittle algorithm): A = L * U.
 * L: unit lower triangular (diagonal = 1). U: upper triangular.
 * L and U stored in-place in LU matrix.
 * pivot[n]: row permutation indices for partial pivoting.
 * Returns number of row swaps (for determinant sign), or -1 if singular. */
int matrix_lu_decompose(double *A, int n, int *pivot);

/* Solve Ax = b using pre-computed LU decomposition.
 * Forward substitution L*y = b, then back substitution U*x = y.
 * pivot[n] from LU decomposition. */
void matrix_lu_solve(const double *LU, const int *pivot,
                     const double *b, double *x, int n);

/* Solve Ax = b: full pipeline (LU decompose + solve).
 * A is modified in-place. Returns 0 on success, -1 if singular. */
int matrix_solve(double *A, const double *b, double *x, int n);

/* Matrix inverse: A_inv = A^{-1}.
 * Uses LU decomposition to solve A * A_inv = I column by column.
 * Returns 0 on success, -1 if singular. */
int matrix_inverse(const double *A, double *A_inv, int n);

/* Matrix condition number estimate (1-norm).
 * cond_1(A) = ||A||_1 * ||A^{-1}||_1.
 * Large (>1e12) => ill-conditioned for double precision. */
double matrix_cond_est(const double *A, int n);

/* Frobenius norm: ||A||_F = sqrt(sum(A[i][j]^2)) */
double matrix_norm_frobenius(const double *A, int rows, int cols);

/* 1-norm (max column sum): ||A||_1 = max_j sum_i |A[i][j]| */
double matrix_norm_1(const double *A, int rows, int cols);

/* Infinity-norm (max row sum): ||A||_inf = max_i sum_j |A[i][j]| */
double matrix_norm_inf(const double *A, int rows, int cols);

/* Eigenvalues of 2x2 real matrix: [a b; c d].
 * Characteristic equation: lambda^2 - trace*lambda + det = 0.
 * Returns 1 if real eigenvalues, 0 if complex conjugate pair. */
int matrix_eigen_2x2(double a, double b, double c, double d,
                     double *lambda1_r, double *lambda1_i,
                     double *lambda2_r, double *lambda2_i);

/* Eigenvalues of symmetric 3x3 real matrix using analytical formula.
 * Returns number of real eigenvalues (always 3 for symmetric). */
int matrix_eigen_sym_3x3(const double *A, double *eigenvalues);

/* Complex matrix operations (for AC analysis) */

/* Allocate complex n x m matrix */
double complex* cmatrix_alloc(int rows, int cols);
void cmatrix_free(double complex *m);
void cmatrix_zero(double complex *A, int rows, int cols);

/* C = A + B */
void cmatrix_add(const double complex *A, const double complex *B,
                 double complex *C, int rows, int cols);

/* C = A * B */
void cmatrix_mul(const double complex *A, const double complex *B,
                 double complex *C, int m, int k, int n);

/* y = A * x */
void cmatrix_vec_mul(const double complex *A, const double complex *x,
                     double complex *y, int m, int n);

/* Complex LU decomposition with partial pivoting */
int cmatrix_lu_decompose(double complex *A, int n, int *pivot);

/* Solve complex Ax = b using LU */
void cmatrix_lu_solve(const double complex *LU, const int *pivot,
                      const double complex *b, double complex *x, int n);

/* Full complex solve */
int cmatrix_solve(double complex *A, const double complex *b,
                  double complex *x, int n);

/* Sparse Matrix Operations (L8: for large circuits)
 * Compressed Row Storage (CRS) format for efficient large-circuit solving. */
typedef struct {
    int    n_rows;      /* number of rows */
    int    n_cols;      /* number of columns */
    int    nnz;         /* number of non-zero entries */
    int    max_nnz;     /* allocated capacity */
    double *values;     /* non-zero values [nnz] */
    int    *col_idx;    /* column indices [nnz] */
    int    *row_ptr;    /* row pointers [n_rows+1] */
} SparseMatrix_t;

/* Initialize empty sparse matrix */
void sparse_init(SparseMatrix_t *sp, int rows, int cols, int capacity);

/* Free sparse matrix */
void sparse_free(SparseMatrix_t *sp);

/* Insert value at position (i,j). Updates existing value if present. */
void sparse_insert(SparseMatrix_t *sp, int i, int j, double value);

/* Get value at (i,j). Returns 0 if not stored. */
double sparse_get(const SparseMatrix_t *sp, int i, int j);

/* y = A * x  for sparse A */
void sparse_vec_mul(const SparseMatrix_t *A, const double *x, double *y);

/* Sparse Gaussian elimination (Markowitz pivot for fill-in reduction).
 * Returns 0 on success. */
int sparse_solve(SparseMatrix_t *A, const double *b, double *x);

/* Circuit-specific: Build sparse matrix from netlist for DC analysis.
 * Use CRS format for efficiency with large circuits. */
int sparse_build_from_circuit(const Circuit_t *ckt, SparseMatrix_t *G,
                              double *I_vec);

#endif /* CIRCUIT_MATRIX_H */
