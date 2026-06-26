/**
 * @file spice_matrix.c
 * @brief Sparse and dense matrix operations for MNA solver
 *
 * Implements CSR sparse matrix, dense matrix, LU decomposition
 * with partial pivoting, and basic BLAS-1 vector operations.
 *
 * Reference: Golub & Van Loan, "Matrix Computations" 4th ed. (2013)
 *   §3.4 — Gaussian elimination and LU factorization
 *   §3.1 — BLAS operations
 *   Saad, "Iterative Methods for Sparse Linear Systems" 2nd ed. (2003)
 *   §3.4 — CSR format
 */

#include "spice_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════
 * Dense Matrix Operations
 * ═══════════════════════════════════════════════════════════════════════ */

spice_dense_matrix_t* spice_dense_matrix_alloc(int32_t nrows, int32_t ncols) {
    if (nrows <= 0 || ncols <= 0) return NULL;

    spice_dense_matrix_t *mat = calloc(1, sizeof(spice_dense_matrix_t));
    if (!mat) return NULL;

    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->ld = nrows;  /* Column-major: leading dimension = number of rows */
    mat->data = calloc((size_t)nrows * ncols, sizeof(double));
    if (!mat->data) {
        free(mat);
        return NULL;
    }
    return mat;
}

void spice_dense_matrix_free(spice_dense_matrix_t *mat) {
    if (!mat) return;
    free(mat->data);
    free(mat);
}

double spice_dense_get(const spice_dense_matrix_t *mat, int32_t row, int32_t col) {
    if (!mat || row < 0 || row >= mat->nrows || col < 0 || col >= mat->ncols)
        return 0.0;
    /* Column-major: A[row + col * ld] */
    return mat->data[row + col * (size_t)mat->ld];
}

void spice_dense_set(spice_dense_matrix_t *mat, int32_t row, int32_t col, double value) {
    if (!mat || row < 0 || row >= mat->nrows || col < 0 || col >= mat->ncols)
        return;
    mat->data[row + col * (size_t)mat->ld] = value;
}

void spice_dense_add(spice_dense_matrix_t *mat, int32_t row, int32_t col, double value) {
    if (!mat || row < 0 || row >= mat->nrows || col < 0 || col >= mat->ncols)
        return;
    mat->data[row + col * (size_t)mat->ld] += value;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Real Vector Operations (BLAS-1 equivalents)
 * ═══════════════════════════════════════════════════════════════════════ */

spice_vector_t* spice_vector_alloc(int32_t size) {
    if (size <= 0) return NULL;
    spice_vector_t *vec = calloc(1, sizeof(spice_vector_t));
    if (!vec) return NULL;
    vec->size = size;
    vec->data = calloc((size_t)size, sizeof(double));
    if (!vec->data) { free(vec); return NULL; }
    return vec;
}

void spice_vector_free(spice_vector_t *vec) {
    if (!vec) return;
    free(vec->data);
    free(vec);
}

spice_complex_vector_t* spice_complex_vector_alloc(int32_t size) {
    if (size <= 0) return NULL;
    spice_complex_vector_t *v = calloc(1, sizeof(spice_complex_vector_t));
    if (!v) return NULL;
    v->size = size;
    v->data = calloc((size_t)size, sizeof(spice_complex_t));
    if (!v->data) { free(v); return NULL; }
    return v;
}

void spice_complex_vector_free(spice_complex_vector_t *v) {
    if (!v) return;
    free(v->data);
    free(v);
}

double spice_vector_dot(const double *a, const double *b, int32_t n) {
    if (!a || !b || n <= 0) return 0.0;
    double sum = 0.0;
    for (int32_t i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

double spice_vector_norm2(const double *v, int32_t n) {
    if (!v || n <= 0) return 0.0;
    double sum_sq = 0.0;
    for (int32_t i = 0; i < n; i++) {
        sum_sq += v[i] * v[i];
    }
    return sqrt(sum_sq);
}

double spice_vector_norm_inf(const double *v, int32_t n) {
    if (!v || n <= 0) return 0.0;
    double max_abs = 0.0;
    for (int32_t i = 0; i < n; i++) {
        double absval = fabs(v[i]);
        if (absval > max_abs) max_abs = absval;
    }
    return max_abs;
}

void spice_vector_scale(double *y, const double *x, double alpha, int32_t n) {
    if (!y || !x || n <= 0) return;
    for (int32_t i = 0; i < n; i++) {
        y[i] = alpha * x[i];
    }
}

void spice_vector_axpy(double *y, const double *x, double alpha, int32_t n) {
    if (!y || !x || n <= 0) return;
    for (int32_t i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

void spice_vector_copy(double *y, const double *x, int32_t n) {
    if (!y || !x || n <= 0) return;
    memcpy(y, x, (size_t)n * sizeof(double));
}

void spice_vector_zero(double *v, int32_t n) {
    if (!v || n <= 0) return;
    memset(v, 0, (size_t)n * sizeof(double));
}

/* ═══════════════════════════════════════════════════════════════════════
 * CSR Sparse Matrix Operations
 * ═══════════════════════════════════════════════════════════════════════ */

spice_csr_matrix_t* spice_csr_matrix_alloc(int32_t nrows, int32_t ncols, int32_t max_nnz) {
    if (nrows <= 0 || ncols <= 0 || max_nnz <= 0) return NULL;

    spice_csr_matrix_t *mat = calloc(1, sizeof(spice_csr_matrix_t));
    if (!mat) return NULL;

    mat->nrows = nrows;
    mat->ncols = ncols;
    mat->nnz = 0;
    mat->nnz_alloc = max_nnz;
    mat->is_symmetric = 0;

    mat->values  = calloc((size_t)max_nnz, sizeof(double));
    mat->columns = calloc((size_t)max_nnz, sizeof(int32_t));
    mat->row_ptr = calloc((size_t)nrows + 1, sizeof(int32_t));

    if (!mat->values || !mat->columns || !mat->row_ptr) {
        free(mat->values);
        free(mat->columns);
        free(mat->row_ptr);
        free(mat);
        return NULL;
    }

    /* Initialize row_ptr: all rows empty */
    for (int32_t i = 0; i <= nrows; i++) {
        mat->row_ptr[i] = 0;
    }

    return mat;
}

void spice_csr_matrix_free(spice_csr_matrix_t *mat) {
    if (!mat) return;
    free(mat->values);
    free(mat->columns);
    free(mat->row_ptr);
    free(mat);
}

/**
 * @brief Compare two entries by (row, col) for sorting
 */
static int entry_cmp(const void *a, const void *b) {
    const spice_sparse_entry_t *ea = (const spice_sparse_entry_t*)a;
    const spice_sparse_entry_t *eb = (const spice_sparse_entry_t*)b;
    if (ea->row != eb->row) return ea->row - eb->row;
    return ea->col - eb->col;
}

int spice_csr_add_entry(spice_csr_matrix_t *mat, int32_t row, int32_t col, double value) {
    if (!mat || row < 0 || row >= mat->nrows || col < 0 || col >= mat->ncols)
        return -1;
    if (mat->nnz >= mat->nnz_alloc) return -1;

    /* Store in coordinate format; finalize will sort and merge */
    mat->values[mat->nnz] = value;
    mat->columns[mat->nnz] = col;
    /* Use row_ptr temporarily to store row indices during assembly */
    mat->row_ptr[mat->nnz] = row;
    mat->nnz++;
    return 0;
}

int spice_csr_finalize(spice_csr_matrix_t *mat) {
    if (!mat || mat->nnz == 0) return 0;

    /* Build temporary entry list */
    spice_sparse_entry_t *entries = calloc((size_t)mat->nnz, sizeof(spice_sparse_entry_t));
    if (!entries) return -1;

    for (int32_t k = 0; k < mat->nnz; k++) {
        entries[k].row   = mat->row_ptr[k];
        entries[k].col   = mat->columns[k];
        entries[k].value = mat->values[k];
    }

    /* Sort by (row, col) */
    qsort(entries, (size_t)mat->nnz, sizeof(spice_sparse_entry_t), entry_cmp);

    /* Merge duplicates and build final CSR */
    int32_t *counts = calloc((size_t)mat->nrows, sizeof(int32_t));
    if (!counts) { free(entries); return -1; }

    /* First pass: count non-zeros per row (merge adjacent duplicates) */
    int32_t unique_nnz = 0;
    for (int32_t k = 0; k < mat->nnz; k++) {
        if (k > 0 && entries[k].row == entries[k-1].row &&
            entries[k].col == entries[k-1].col) {
            /* Duplicate — merge into previous */
            entries[unique_nnz - 1].value += entries[k].value;
        } else {
            entries[unique_nnz] = entries[k];
            counts[entries[k].row]++;
            unique_nnz++;
        }
    }

    /* Build row_ptr */
    mat->row_ptr[0] = 0;
    for (int32_t i = 0; i < mat->nrows; i++) {
        mat->row_ptr[i + 1] = mat->row_ptr[i] + counts[i];
    }

    /* Store values and columns in CSR order */
    int32_t *pos = calloc((size_t)mat->nrows, sizeof(int32_t));
    if (!pos) { free(counts); free(entries); return -1; }
    memcpy(pos, mat->row_ptr, (size_t)mat->nrows * sizeof(int32_t));

    /* Zero out value/column arrays for clean rebuild */
    memset(mat->values,  0, (size_t)mat->nnz_alloc * sizeof(double));
    memset(mat->columns, 0, (size_t)mat->nnz_alloc * sizeof(int32_t));

    for (int32_t k = 0; k < unique_nnz; k++) {
        int32_t r = entries[k].row;
        int32_t idx = pos[r];
        mat->values[idx]  = entries[k].value;
        mat->columns[idx] = entries[k].col;
        pos[r]++;
    }

    mat->nnz = unique_nnz;

    free(entries);
    free(counts);
    free(pos);
    return 0;
}

void spice_csr_matvec(const spice_csr_matrix_t *A, const spice_vector_t *x, spice_vector_t *y) {
    if (!A || !x || !y) return;
    if (x->size < A->ncols || y->size < A->nrows) return;

    /* y = A * x */
    for (int32_t i = 0; i < A->nrows; i++) {
        double sum = 0.0;
        for (int32_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; k++) {
            sum += A->values[k] * x->data[A->columns[k]];
        }
        y->data[i] = sum;
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * LU Decomposition with Partial Pivoting (Dense)
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Gaussian elimination with partial pivoting to compute PA = LU
 *
 * Factorization is performed in-place. On output:
 *   - Upper triangular part (including diagonal) contains U
 *   - Strictly lower triangular part contains the multipliers of L
 *     (with implicit unit diagonal)
 *
 * Pivot strategy: at step k, find the element of largest absolute value
 * in column k below (and including) the diagonal. Swap rows if needed.
 * This ensures numerical stability by avoiding small pivots.
 *
 * @param A     Dense n×n matrix in column-major order (overwritten with L+U)
 * @param pivot Array of pivot indices: row i was swapped with row pivot[i]
 * @param n     Matrix order
 * @param ld    Leading dimension (ld >= n)
 * @return 0 on success, k>0 if zero pivot encountered at step k
 *
 * Reference: Golub & Van Loan (2013) Algorithm 3.4.1
 */
int spice_dense_lu_factor(double *A, int32_t *pivot, int32_t n, int32_t ld) {
    if (!A || !pivot || n <= 0 || ld < n) return -1;

    /* Initialize pivot array */
    for (int32_t i = 0; i < n; i++) {
        pivot[i] = i;
    }

    for (int32_t k = 0; k < n; k++) {
        /* Find pivot: element of maximum absolute value in column k */
        double max_abs = 0.0;
        int32_t max_row = k;

        for (int32_t i = k; i < n; i++) {
            double abs_val = fabs(A[i + k * (size_t)ld]);
            if (abs_val > max_abs) {
                max_abs = abs_val;
                max_row = i;
            }
        }

        /* Check for singularity */
        if (max_abs < 1e-30) {
            return k + 1;  /* Singular at step k */
        }

        /* Pivot: swap rows k and max_row */
        if (max_row != k) {
            pivot[k] = max_row;
            for (int32_t j = 0; j < n; j++) {
                double tmp = A[k + j * (size_t)ld];
                A[k + j * (size_t)ld] = A[max_row + j * (size_t)ld];
                A[max_row + j * (size_t)ld] = tmp;
            }
        } else {
            pivot[k] = k;
        }

        /* Compute multipliers and eliminate below diagonal */
        double pivot_val = A[k + k * (size_t)ld];

        for (int32_t i = k + 1; i < n; i++) {
            double multiplier = A[i + k * (size_t)ld] / pivot_val;
            A[i + k * (size_t)ld] = multiplier;  /* Store L multiplier */

            /* Update row i */
            for (int32_t j = k + 1; j < n; j++) {
                A[i + j * (size_t)ld] -= multiplier * A[k + j * (size_t)ld];
            }
        }
    }

    return 0;
}

/**
 * @brief Solve Ax = b given LU factors from spice_dense_lu_factor
 *
 * Step 1: Forward substitution — solve L y = P b
 * Step 2: Backward substitution — solve U x = y
 *
 * The L matrix has unit diagonal (implicit). The multipliers are
 * stored in the strictly lower triangular part of A.
 *
 * @param A     LU factors (overwritten matrix from factorization)
 * @param pivot Pivot array
 * @param b     RHS vector (overwritten with solution x)
 * @param n     System order
 * @param ld    Leading dimension
 *
 * Reference: Golub & Van Loan (2013) §3.4.1
 */
void spice_dense_lu_solve(const double *A, const int32_t *pivot, double *b,
                          int32_t n, int32_t ld) {
    if (!A || !pivot || !b || n <= 0) return;

    /* Step 1: Apply pivot permutations to b (forward order) */
    /* Note: pivot[k] gives the row that was swapped with row k during LU.
       We need P * b, which means b[pivot[i]] goes to position i. */
    /* Actually, for the standard LAPACK convention:
       pivot[k] = j means row k <-> row j was done at step j.
       The forward solve needs to apply permutations first. */

    /* Apply row interchanges to b: for k=0..n-2, swap b[k] and b[pivot[k]] */
    for (int32_t k = 0; k < n - 1; k++) {
        int32_t pk = pivot[k];
        if (pk != k) {
            double tmp = b[k];
            b[k] = b[pk];
            b[pk] = tmp;
        }
    }

    /* Step 2: Forward substitution L y = b (y overwrites b)
       L has unit diagonal, multipliers in lower part. */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < i; j++) {
            b[i] -= A[i + j * (size_t)ld] * b[j];
        }
        /* b[i] already has the right value since L[i][i] = 1 */
    }

    /* Step 3: Back substitution U x = y (x overwrites b)
       U is upper triangular (including diagonal). */
    for (int32_t i = n - 1; i >= 0; i--) {
        for (int32_t j = i + 1; j < n; j++) {
            b[i] -= A[i + j * (size_t)ld] * b[j];
        }
        b[i] /= A[i + i * (size_t)ld];
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * Complex LU Solution (for AC analysis)
 * ═══════════════════════════════════════════════════════════════════════ */

int spice_complex_lu_solve(spice_complex_t *A, int32_t *pivot, spice_complex_t *b, int32_t n) {
    if (!A || !pivot || !b || n <= 0) return -1;

    /* Complex Gaussian elimination with partial pivoting */
    for (int32_t k = 0; k < n; k++) {
        /* Find pivot */
        double max_abs = 0.0;
        int32_t max_row = k;
        for (int32_t i = k; i < n; i++) {
            double cabs_val = cabs(A[i + k * (size_t)n]);
            if (cabs_val > max_abs) {
                max_abs = cabs_val;
                max_row = i;
            }
        }
        if (max_abs < 1e-30) return k + 1;

        pivot[k] = max_row;
        if (max_row != k) {
            /* Swap rows */
            for (int32_t j = 0; j < n; j++) {
                spice_complex_t tmp = A[k + j * (size_t)n];
                A[k + j * (size_t)n] = A[max_row + j * (size_t)n];
                A[max_row + j * (size_t)n] = tmp;
            }
            spice_complex_t tmpb = b[k];
            b[k] = b[max_row];
            b[max_row] = tmpb;
        }

        spice_complex_t pivot_val = A[k + k * (size_t)n];

        /* Eliminate */
        for (int32_t i = k + 1; i < n; i++) {
            spice_complex_t mult = A[i + k * (size_t)n] / pivot_val;
            A[i + k * (size_t)n] = mult;
            for (int32_t j = k + 1; j < n; j++) {
                A[i + j * (size_t)n] -= mult * A[k + j * (size_t)n];
            }
        }
    }

    /* Forward substitution */
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < i; j++) {
            b[i] -= A[i + j * (size_t)n] * b[j];
        }
    }

    /* Back substitution */
    for (int32_t i = n - 1; i >= 0; i--) {
        for (int32_t j = i + 1; j < n; j++) {
            b[i] -= A[i + j * (size_t)n] * b[j];
        }
        b[i] /= A[i + i * (size_t)n];
    }

    return 0;
}

/* ═══════════════════════════════════════════════════════════════════════
 * Condition Number Estimation
 * ═══════════════════════════════════════════════════════════════════════ */

/**
 * @brief Compute the 1-norm of a dense matrix
 *
 * ||A||_1 = max_j Σ_i |A[i][j]|  (maximum column sum)
 */
static double matrix_norm1(const double *A, int32_t n) {
    double max_col_sum = 0.0;
    for (int32_t j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (int32_t i = 0; i < n; i++) {
            col_sum += fabs(A[i + j * (size_t)n]);
        }
        if (col_sum > max_col_sum) max_col_sum = col_sum;
    }
    return max_col_sum;
}

double spice_condition_number(const double *A, int32_t n) {
    if (!A || n <= 0) return -1.0;

    /* Compute ||A||_1 */
    double normA = matrix_norm1(A, n);
    if (normA == 0.0) return -1.0;

    /* Compute LU factorization to estimate ||A^{-1}||_1 */
    /* Make a copy of A since LU overwrites it */
    double *A_copy = calloc((size_t)n * n, sizeof(double));
    if (!A_copy) return -1.0;
    memcpy(A_copy, A, (size_t)n * n * sizeof(double));

    int32_t *pivot = calloc((size_t)n, sizeof(int32_t));
    if (!pivot) { free(A_copy); return -1.0; }

    int info = spice_dense_lu_factor(A_copy, pivot, n, n);
    if (info > 0) { free(pivot); free(A_copy); return -1.0; }

    /* Estimate ||A^{-1}||_1 using inverse power method on A^T
     * (Higham's algorithm for 1-norm condition estimation) */
    /* Simplified: use 1 iteration */
    double *y = calloc((size_t)n, sizeof(double));
    double *x = calloc((size_t)n, sizeof(double));
    if (!y || !x) { free(pivot); free(A_copy); free(y); free(x); return -1.0; }

    /* Start with x = (1/n, 1/n, ..., 1/n) */
    for (int32_t i = 0; i < n; i++) x[i] = 1.0 / n;

    /* Solve A^T y = x, then A z = sign(y), estimate = ||z||_1 */
    /* Solve A^T y = x: solve y^T A = x^T, i.e., forward/back on A */
    /* Use the LU decomposition: solve A^T y = x → y = A^{-T} x */
    /* A^{-T} = (L^{-T})(U^{-T}), need to solve U^T L^T y = x */
    /* For simplicity, solve 2 systems: A z = e_j, estimate = max ||z||_1 */

    double norm_inv_est = 0.0;
    for (int32_t trial = 0; trial < 3; trial++) {
        /* Solve A z = x */
        memcpy(y, x, (size_t)n * sizeof(double));
        spice_dense_lu_solve(A_copy, pivot, y, n, n);

        /* Update estimate */
        double z_norm1 = 0.0;
        for (int32_t i = 0; i < n; i++) z_norm1 += fabs(y[i]);
        if (z_norm1 > norm_inv_est) norm_inv_est = z_norm1;

        /* Set x = sign(y) for next iteration (cautious: avoid zeros) */
        for (int32_t i = 0; i < n; i++) {
            x[i] = (y[i] >= 0) ? 1.0 : -1.0;
        }
    }

    free(A_copy);
    free(pivot);
    free(y);
    free(x);

    return normA * norm_inv_est;
}
