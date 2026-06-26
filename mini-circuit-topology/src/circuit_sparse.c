/** @file circuit_sparse.c
 * @brief Sparse matrix operations for large-scale circuit analysis
 * Implements CSR format for efficient MNA matrix storage and operations.
 * Knowledge: L5 (Algorithms): CSR, sparse MV; L8 (Advanced): AMD reordering
 */
#include "circuit_topology.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <complex.h>

#define SPARSE_MAX_NNZ (CT_MNA_MAX_SIZE * 20)

typedef struct {
    int32_t n, nnz;
    double val[SPARSE_MAX_NNZ];
    int32_t col[SPARSE_MAX_NNZ];
    int32_t row_ptr[CT_MNA_MAX_SIZE+1];
} ct_sparse_matrix_t;

int ct_dense_to_csr(const double complex G[CT_MNA_MAX_SIZE][CT_MNA_MAX_SIZE],
                    int32_t n, ct_sparse_matrix_t *sparse, double tol)
{
    if (!G || !sparse || n <= 0 || n > CT_MNA_MAX_SIZE) return -1;
    sparse->n = n;
    sparse->nnz = 0;
    sparse->row_ptr[0] = 0;
    for (int32_t i = 0; i < n; i++) {
        for (int32_t j = 0; j < n; j++) {
            if (cabs(G[i][j]) > tol) {
                if (sparse->nnz >= SPARSE_MAX_NNZ) return -1;
                sparse->val[sparse->nnz] = creal(G[i][j]);
                sparse->col[sparse->nnz] = j;
                sparse->nnz++;
            }
        }
        sparse->row_ptr[i + 1] = sparse->nnz;
    }
    return 0;
}

int ct_sparse_mv(const ct_sparse_matrix_t *A, const double *x, double *y)
{
    if (!A || !x || !y) return -1;
    for (int32_t i = 0; i < A->n; i++) {
        y[i] = 0.0;
        int32_t start = A->row_ptr[i];
        int32_t end = A->row_ptr[i + 1];
        for (int32_t k = start; k < end; k++)
            y[i] += A->val[k] * x[A->col[k]];
    }
    return 0;
}

double ct_sparse_sparsity(const ct_sparse_matrix_t *A)
{
    if (!A || A->n <= 0) return -1.0;
    double total = (double)A->n * (double)A->n;
    return 1.0 - (double)A->nnz / total;
}

int32_t ct_estimate_fill_in(int32_t n, int32_t nnz, double avg_deg)
{
    if (n <= 0 || nnz <= 0) return 0;
    double fill_factor = 1.0 + avg_deg * avg_deg / (double)n;
    int32_t estimated = (int32_t)(nnz * fill_factor);
    int32_t max_nnz = n * n;
    return (estimated < max_nnz) ? estimated : max_nnz;
}