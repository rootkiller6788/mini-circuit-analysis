/**
 * @file spice_matrix.h
 * @brief Sparse matrix data structures for Modified Nodal Analysis (MNA)
 *
 * Knowledge coverage:
 *   L3 (Math Structures): Sparse matrix CSR format, complex matrix, dense vector
 *   L4 (Fundamental Laws): MNA formulation embeds KCL + KVL + branch equations
 *   L5 (Algorithms): LU decomposition with partial pivoting
 *
 * The Modified Nodal Analysis (MNA) formulation (Ho, Ruehli, Brennan 1975)
 * produces matrices of the form:
 *
 *   [ G   B ] [ v ]   [ i ]
 *   [ C   D ] [ j ] = [ e ]
 *
 * where G is the conductance matrix, B and C are source connection matrices,
 * D is zeros (or ω-related for inductors in AC), v is node voltages,
 * and j is branch currents.
 *
 * Course alignment:
 *   Berkeley EE16B (Circuits): MNA formulation, nodal analysis
 *   Georgia Tech ECE 6350 (EM): sparse linear system solvers
 *   ETH 227-0455 (EM): matrix methods for field/circuit problems
 */

#ifndef SPICE_MATRIX_H
#define SPICE_MATRIX_H

#include <stddef.h>
#include <stdint.h>
#include <complex.h>

/** Maximum size of the MNA matrix (pre-allocated) */
#define SPICE_MAX_MNA_SIZE 2048

/** Maximum non-zero entries in the sparse matrix */
#define SPICE_MAX_NNZ      65536

/* ── L3: Mathematical Structures ──────────────────────────────────── */

/**
 * @brief Complex number type using C99 _Complex double
 *
 * Used in AC analysis for frequency-domain phasor calculations.
 * real(z) = creal(z), imag(z) = cimag(z).
 */
typedef double complex spice_complex_t;

/**
 * @brief Single entry in a coordinate-format sparse matrix
 *
 * Coordinates (row, col) with a double value.
 * Used during matrix assembly before conversion to CSR.
 */
typedef struct {
    int32_t row;       /**< Row index (0-based) */
    int32_t col;       /**< Column index (0-based) */
    double  value;     /**< Matrix entry value */
} spice_sparse_entry_t;

/**
 * @brief Compressed Sparse Row (CSR) matrix
 *
 * The dominant format for sparse MNA matrices. Provides O(1)
 * row access and efficient matrix-vector multiplication.
 *
 * For row i:
 *   columns[k] for k in [row_ptr[i], row_ptr[i+1]) gives column indices
 *   values[k]  for k in [row_ptr[i], row_ptr[i+1]) gives entry values
 *
 * Reference: Saad, "Iterative Methods for Sparse Linear Systems" (2003) §3.4
 */
typedef struct {
    int32_t   nrows;        /**< Number of rows                        */
    int32_t   ncols;        /**< Number of columns                     */
    int32_t   nnz;          /**< Number of non-zero entries            */
    int32_t   nnz_alloc;    /**< Allocated capacity for entries        */
    double*   values;       /**< Array of non-zero values (length nnz) */
    int32_t*  columns;      /**< Column index for each value           */
    int32_t*  row_ptr;      /**< Row pointer array (length nrows+1)    */
    int32_t   is_symmetric; /**< 1 if matrix is structurally symmetric  */
} spice_csr_matrix_t;

/**
 * @brief Dense matrix for LU decomposition and small dense sub-blocks
 *
 * MNA matrices for small to medium circuits (<500 nodes) are often
 * solved with dense LU for numerical robustness.
 */
typedef struct {
    int32_t  nrows;         /**< Number of rows                        */
    int32_t  ncols;         /**< Number of columns                     */
    int32_t  ld;            /**< Leading dimension (>= nrows)          */
    double*  data;          /**< Data in column-major order (BLAS convention) */
} spice_dense_matrix_t;

/**
 * @brief Dense vector of real numbers
 *
 * Represents node voltages, RHS vectors, and branch currents.
 */
typedef struct {
    int32_t  size;          /**< Number of elements                    */
    double*  data;          /**< Data array                            */
} spice_vector_t;

/**
 * @brief Dense vector of complex numbers (for AC analysis)
 */
typedef struct {
    int32_t          size;
    spice_complex_t* data;
} spice_complex_vector_t;

/**
 * @brief LU decomposition result
 *
 * L and U factors stored in-place in the original matrix,
 * with pivot information for row interchanges.
 *
 * For a matrix A of size n, the LU decomposition with partial
 * pivoting produces:
 *   P * A = L * U
 * where P is stored implicitly in pivot[]:
 *   swap(row i, pivot[i]) was applied during factorization.
 */
typedef struct {
    spice_dense_matrix_t  matrix; /**< Combined L+U factors (in-place)  */
    int32_t*              pivot;  /**< Pivot indices (row permutations) */
    int32_t               size;   /**< Matrix order                     */
    int32_t               info;   /**< 0 = success, >0 = singular       */
} spice_lu_factor_t;

/* ── L5: Matrix Operations ────────────────────────────────────────── */

/**
 * @brief Allocate a dense matrix
 *
 * @param nrows Number of rows
 * @param ncols Number of columns
 * @return Pointer to allocated matrix, or NULL on failure
 *
 * Leading dimension is set to nrows. Data is zero-initialized.
 * Complexity: O(nrows * ncols)
 */
spice_dense_matrix_t* spice_dense_matrix_alloc(int32_t nrows, int32_t ncols);

/**
 * @brief Free a dense matrix
 */
void spice_dense_matrix_free(spice_dense_matrix_t *mat);

/**
 * @brief Get element A[row][col] from dense (col-major) matrix
 *
 * @param mat Matrix
 * @param row Row index (0-based)
 * @param col Column index (0-based)
 * @return A[row][col]
 *
 * Complexity: O(1)
 */
double spice_dense_get(const spice_dense_matrix_t *mat, int32_t row, int32_t col);

/**
 * @brief Set element A[row][col] in dense (col-major) matrix
 *
 * Complexity: O(1)
 */
void spice_dense_set(spice_dense_matrix_t *mat, int32_t row, int32_t col, double value);

/**
 * @brief Add value to element A[row][col] (used in stamping)
 *
 * Complexity: O(1)
 */
void spice_dense_add(spice_dense_matrix_t *mat, int32_t row, int32_t col, double value);

/**
 * @brief Allocate a real vector
 *
 * @param size Number of elements
 * @return Allocated vector or NULL
 */
spice_vector_t* spice_vector_alloc(int32_t size);

/**
 * @brief Free a real vector
 */
void spice_vector_free(spice_vector_t *vec);

/**
 * @brief Allocate a complex vector for AC analysis
 */
spice_complex_vector_t* spice_complex_vector_alloc(int32_t size);

/**
 * @brief Free a complex vector
 */
void spice_complex_vector_free(spice_complex_vector_t *vec);

/**
 * @brief Allocate a CSR sparse matrix
 *
 * @param nrows Number of rows
 * @param ncols Number of columns
 * @param max_nnz Maximum non-zero entries (for pre-allocation)
 * @return Allocated sparse matrix or NULL
 *
 * Complexity: O(max_nnz)
 */
spice_csr_matrix_t* spice_csr_matrix_alloc(int32_t nrows, int32_t ncols, int32_t max_nnz);

/**
 * @brief Free a CSR sparse matrix
 */
void spice_csr_matrix_free(spice_csr_matrix_t *mat);

/**
 * @brief Add an entry to CSR matrix during assembly phase
 *
 * @param mat   CSR matrix
 * @param row   Row index
 * @param col   Column index
 * @param value Value to add
 * @return 0 on success, -1 if out of space
 *
 * Entries must be added in order of increasing (row, col)
 * for correct CSR conversion. This function accumulates entries
 * in a temporary coordinate list.
 */
int spice_csr_add_entry(spice_csr_matrix_t *mat, int32_t row, int32_t col, double value);

/**
 * @brief Finalize CSR matrix (sort entries, build row_ptr)
 *
 * After all entries are added, this function sorts by (row, col),
 * merges duplicate entries, and builds the row_ptr array.
 *
 * @param mat CSR matrix to finalize
 * @return 0 on success
 */
int spice_csr_finalize(spice_csr_matrix_t *mat);

/**
 * @brief Sparse matrix-vector multiply: y = A * x
 *
 * @param A CSR matrix
 * @param x Input vector
 * @param y Output vector (pre-allocated, size = A->nrows)
 *
 * Complexity: O(nnz)
 *
 * Reference: Saad (2003) Algorithm 3.1
 */
void spice_csr_matvec(const spice_csr_matrix_t *A, const spice_vector_t *x, spice_vector_t *y);

/**
 * @brief Compute LU decomposition with partial pivoting
 *
 * Implements Gaussian elimination with row pivoting.
 * L and U factors overwrite the input matrix.
 *
 * @param A    Dense matrix (n x n) — overwritten with L and U
 * @param pivot Output pivot array (size n), caller-allocated
 * @param n    Matrix order
 * @return 0 on success, k>0 if matrix is singular (zero pivot at step k)
 *
 * Reference: Golub & Van Loan, "Matrix Computations" (2013) §3.4
 * Complexity: O(n^3)
 */
int spice_dense_lu_factor(double *A, int32_t *pivot, int32_t n, int32_t ld);

/**
 * @brief Solve A*x = b using LU factors from spice_dense_lu_factor
 *
 * @param A     LU factors (L+U in-place from factorization)
 * @param pivot Pivot array from factorization
 * @param b     Right-hand side vector (size n) — overwritten with solution x
 * @param n     System order
 * @param ld    Leading dimension of A
 *
 * Step 1: Forward substitution: L * y = P * b
 * Step 2: Backward substitution: U * x = y
 *
 * Reference: Golub & Van Loan (2013) §3.4.1
 * Complexity: O(n^2)
 */
void spice_dense_lu_solve(const double *A, const int32_t *pivot, double *b,
                          int32_t n, int32_t ld);

/**
 * @brief Solve complex linear system A*x = b using LU
 *
 * For AC analysis, the MNA matrix is complex-valued.
 * This solves the complex system directly.
 *
 * @param A Complex matrix (n x n, col-major) — overwritten with LU
 * @param pivot Pivot array (size n)
 * @param b Complex RHS — overwritten with solution
 * @param n System order
 * @return 0 on success
 *
 * Complexity: O(n^3)
 */
int spice_complex_lu_solve(spice_complex_t *A, int32_t *pivot, spice_complex_t *b, int32_t n);

/**
 * @brief Estimate the condition number of a matrix
 *
 * Uses the 1-norm. cond(A) = ||A||_1 * ||A^{-1}||_1
 *
 * @param A Matrix (n x n)
 * @param n Order
 * @return Estimated condition number (>= 1.0), or -1 on error
 *
 * Reference: Higham, "Accuracy and Stability of Numerical Algorithms" (2002) §15
 */
double spice_condition_number(const double *A, int32_t n);

/**
 * @brief Compute dot product of two vectors
 *
 * @param a First vector
 * @param b Second vector
 * @param n Length
 * @return a·b = Σ a[i]*b[i]
 *
 * Complexity: O(n)
 */
double spice_vector_dot(const double *a, const double *b, int32_t n);

/**
 * @brief Compute the L2 (Euclidean) norm of a vector
 *
 * @param v Vector
 * @param n Length
 * @return ||v||_2 = sqrt(Σ v[i]^2)
 *
 * Complexity: O(n)
 */
double spice_vector_norm2(const double *v, int32_t n);

/**
 * @brief Compute the infinity norm of a vector
 *
 * @param v Vector
 * @param n Length
 * @return ||v||_∞ = max |v[i]|
 *
 * Complexity: O(n)
 */
double spice_vector_norm_inf(const double *v, int32_t n);

/**
 * @brief Scale a vector: y = alpha * x
 *
 * @param y Output vector
 * @param x Input vector
 * @param alpha Scalar
 * @param n Length
 *
 * Complexity: O(n)
 */
void spice_vector_scale(double *y, const double *x, double alpha, int32_t n);

/**
 * @brief Vector addition: y = alpha * x + y (daxpy)
 *
 * @param y Output/accumulator vector
 * @param x Input vector
 * @param alpha Scalar multiplier
 * @param n Length
 *
 * Complexity: O(n)
 */
void spice_vector_axpy(double *y, const double *x, double alpha, int32_t n);

/**
 * @brief Copy vector: y = x
 *
 * Complexity: O(n)
 */
void spice_vector_copy(double *y, const double *x, int32_t n);

/**
 * @brief Set all elements of a vector to zero
 *
 * Complexity: O(n)
 */
void spice_vector_zero(double *v, int32_t n);

#endif /* SPICE_MATRIX_H */
