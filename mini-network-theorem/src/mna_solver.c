/**
 * @file mna_solver.c
 * @brief Modified Nodal Analysis (MNA) — SPICE-Compatible Circuit Solver
 *
 * Implements the Modified Nodal Analysis method, the industry-standard
 * approach used by SPICE and all commercial circuit simulators.
 *
 * Unlike pure nodal analysis, MNA handles voltage sources natively
 * (without source transformation) by adding their currents as unknowns
 * and their defining equations as constraints.
 *
 * System structure:
 *   ┌         ┐ ┌   ┐   ┌   ┐
 *   │ G    B  │ │ v │ = │ i │
 *   │ C    D  │ │ j │   │ e │
 *   └         ┘ └   ┘   └   ┘
 *
 * G (n×n): Conductance matrix (same as pure nodal)
 * B (n×m): Voltage source incidence (±1 entries)
 * C (m×n): Essentially B^T for independent sources
 * D (m×m): Zero matrix for independent sources
 * v (n×1): Unknown node voltages
 * j (m×1): Unknown voltage source currents
 * i (n×1): Known current source injections
 * e (m×1): Known voltage source values
 *
 * Knowledge Coverage:
 *   L5 - Algorithms: Modified Nodal Analysis
 *   L3 - Mathematical Structures: Augmented matrix system
 *   L8 - Advanced Topics: SPICE-compatible simulation core
 *
 * Reference:
 *   Ho, Ruehli, Brennan (1975) "The modified nodal approach to network analysis"
 *   IEEE Trans. Circuits and Systems, vol. CAS-22, pp. 504-509
 *   Nagel & Pederson (1973) "SPICE — Simulation Program with Integrated
 *   Circuit Emphasis"
 *
 * Course Mapping:
 *   Berkeley EE105: SPICE simulation, MNA formulation
 *   MIT 6.002: Circuit simulation methods
 *   Stanford EE247: EDA tools and circuit simulation
 */

#include "network_theorem.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ==========================================================================
 * L5: Build MNA System Matrix
 * ==========================================================================
 *
 * Constructs the MNA system for a given circuit topology.
 *
 * Algorithm:
 *   1. Count non-ground nodes → n
 *   2. Count independent voltage sources → m
 *   3. Allocate system of size dim = n + m
 *   4. Fill G block: conductance matrix (nodal admittance)
 *   5. Fill B block: voltage source connection stamps
 *      - B[node_from][vsrc_idx] = +1
 *      - B[node_to][vsrc_idx]   = -1
 *   6. Fill C block: B^T (for independent sources)
 *   7. Fill D block: zeros (for independent sources)
 *   8. Fill RHS b vector: current source contributions + voltage source values
 *
 * Returns: 0 on success, -1 on error.
 */
int mna_build_system(const CircuitTopology *ckt, MNASystem *sys) {
    if (!ckt || !sys) return -1;

    uint32_t n_nodes = ckt->num_nodes;
    if (n_nodes == 0) return -1;

    /* Count voltage sources */
    uint32_t n_vsrc = 0;
    for (uint32_t i = 0; i < ckt->num_elements; i++) {
        if (ckt->elements[i].type == ELEM_VOLTAGE_SOURCE) {
            n_vsrc++;
        }
    }

    uint32_t dim = n_nodes + n_vsrc;
    sys->n_nodes = n_nodes;
    sys->n_vsources = n_vsrc;
    sys->dim = dim;

    /* Allocate system matrix and vectors */
    sys->A = (double *)calloc((size_t)dim * dim, sizeof(double));
    sys->b = (double *)calloc(dim, sizeof(double));
    sys->x = (double *)calloc(dim, sizeof(double));

    if (!sys->A || !sys->b || !sys->x) {
        mna_free(sys);
        return -1;
    }

    /* --- Fill G block (n×n conductance matrix) --- */
    for (uint32_t e = 0; e < ckt->num_elements; e++) {
        const CircuitElement *elem = &ckt->elements[e];
        uint32_t i = elem->node_from;
        uint32_t j = elem->node_to;

        double G = 0.0;
        switch (elem->type) {
            case ELEM_RESISTOR:
                if (fabs(elem->value) > 1e-15) {
                    G = 1.0 / elem->value;
                }
                break;
            case ELEM_CAPACITOR:
                /* For DC: open circuit, G=0. For AC: would use jωC → complex.
                 * This implementation focuses on DC resistive circuits. */
                continue;
            case ELEM_INDUCTOR:
                /* For DC: short circuit → needs special handling.
                 * For now, inductor shorts are treated as very large G. */
                G = 1e12;
                break;
            case ELEM_SHORT:
                G = 1e12;
                break;
            default:
                continue;
        }

        if (fabs(G) < 1e-15) continue;

        /* Stamp into conductance matrix */
        if (i < n_nodes) sys->A[i * dim + i] += G;
        if (j < n_nodes) sys->A[j * dim + j] += G;
        if (i < n_nodes && j < n_nodes) {
            sys->A[i * dim + j] -= G;
            sys->A[j * dim + i] -= G;
        }
    }

    /* --- Stamp voltage sources into B and C blocks --- */
    uint32_t vsrc_idx = 0;
    for (uint32_t e = 0; e < ckt->num_elements; e++) {
        const CircuitElement *elem = &ckt->elements[e];
        if (elem->type != ELEM_VOLTAGE_SOURCE) continue;

        uint32_t i = elem->node_from;
        uint32_t j = elem->node_to;
        double V_src = elem->value;

        uint32_t col = n_nodes + vsrc_idx;

        /* B block: stamp ±1 for voltage source connections */
        if (i < n_nodes) {
            sys->A[i * dim + col] = 1.0;   /* B[i][vsrc] = +1 */
            sys->A[col * dim + i] = 1.0;   /* C[vsrc][i] = +1 (B^T) */
        }
        if (j < n_nodes) {
            sys->A[j * dim + col] = -1.0;   /* B[j][vsrc] = -1 */
            sys->A[col * dim + j] = -1.0;   /* C[vsrc][j] = -1 */
        }

        /* RHS: voltage source value */
        sys->b[col] = V_src;

        vsrc_idx++;
    }

    /* --- Stamp current sources into RHS --- */
    for (uint32_t e = 0; e < ckt->num_elements; e++) {
        const CircuitElement *elem = &ckt->elements[e];
        if (elem->type != ELEM_CURRENT_SOURCE) continue;

        uint32_t i = elem->node_from;
        uint32_t j = elem->node_to;
        double Is = elem->value;

        /* Current source: enters node j, leaves node i */
        if (i < n_nodes) sys->b[i] -= Is;
        if (j < n_nodes) sys->b[j] += Is;
    }

    /* Ground node 0: set row 0 = [1,0,...,0], b[0] = 0 */
    for (uint32_t col = 0; col < dim; col++) {
        sys->A[0 * dim + col] = 0.0;
    }
    for (uint32_t row = 0; row < dim; row++) {
        sys->A[row * dim + 0] = 0.0;
    }
    sys->A[0 * dim + 0] = 1.0;
    sys->b[0] = 0.0;

    return 0;
}

/* ==========================================================================
 * L5: Solve MNA System
 * ==========================================================================
 *
 * Solves the augmented linear system Ax = b using Gaussian elimination.
 * The solution x contains:
 *   x[0..n-1] = node voltages
 *   x[n..n+m-1] = voltage source currents
 *
 * Returns: 0 on success, -1 if singular.
 */
int mna_solve(MNASystem *sys) {
    if (!sys || !sys->A || !sys->b || !sys->x) return -1;

    /* Copy A and b since Gaussian elimination modifies them */
    uint32_t dim = sys->dim;
    double *A_copy = (double *)malloc((size_t)dim * dim * sizeof(double));
    double *b_copy = (double *)malloc(dim * sizeof(double));

    if (!A_copy || !b_copy) {
        free(A_copy); free(b_copy);
        return -1;
    }

    memcpy(A_copy, sys->A, (size_t)dim * dim * sizeof(double));
    memcpy(b_copy, sys->b, dim * sizeof(double));

    int ret = gaussian_elimination(A_copy, b_copy, dim, sys->x);

    free(A_copy);
    free(b_copy);
    return ret;
}

/* ==========================================================================
 * L5: Free MNA System Memory
 * ========================================================================== */
void mna_free(MNASystem *sys) {
    if (!sys) return;
    free(sys->A);
    free(sys->b);
    free(sys->x);
    sys->A = NULL;
    sys->b = NULL;
    sys->x = NULL;
    sys->dim = 0;
    sys->n_nodes = 0;
    sys->n_vsources = 0;
}

/* ==========================================================================
 * L8: Sparse Matrix Representation for Large Circuits
 * ==========================================================================
 *
 * For circuits with thousands of nodes, the MNA matrix is highly sparse
 * (each node typically connects to only 2-5 other nodes). Using dense
 * Gaussian elimination wastes O(n³) time and O(n²) memory.
 *
 * This section defines a Compressed Sparse Row (CSR) format for efficient
 * storage and partial solution for large-scale circuit analysis.
 *
 * CSR format:
 *   - values[]: nonzero entries in row-major order
 *   - col_idx[]: column index for each nonzero
 *   - row_ptr[i]: start index in values/col_idx for row i
 *
 * Reference: Saad "Iterative Methods for Sparse Linear Systems" (2003)
 */

/**
 * Convert a dense MNA matrix to CSR format.
 * Only stores entries with |value| > 1e-15.
 */
SparseMatrixCSR *dense_to_csr(const double *A, uint32_t n) {
    if (!A || n == 0) return NULL;

    /* First pass: count nonzeros */
    uint32_t nnz = 0;
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < n; j++) {
            if (fabs(A[i * n + j]) > 1e-15) nnz++;
        }
    }

    SparseMatrixCSR *csr = (SparseMatrixCSR *)malloc(sizeof(SparseMatrixCSR));
    if (!csr) return NULL;

    csr->nnz = nnz;
    csr->n_rows = n;
    csr->n_cols = n;
    csr->values  = (double *)malloc(nnz * sizeof(double));
    csr->col_idx = (uint32_t *)malloc(nnz * sizeof(uint32_t));
    csr->row_ptr = (uint32_t *)calloc(n + 1, sizeof(uint32_t));

    if (!csr->values || !csr->col_idx || !csr->row_ptr) {
        free(csr->values); free(csr->col_idx); free(csr->row_ptr);
        free(csr);
        return NULL;
    }

    /* Second pass: fill CSR arrays */
    uint32_t pos = 0;
    for (uint32_t i = 0; i < n; i++) {
        csr->row_ptr[i] = pos;
        for (uint32_t j = 0; j < n; j++) {
            double val = A[i * n + j];
            if (fabs(val) > 1e-15) {
                csr->values[pos] = val;
                csr->col_idx[pos] = j;
                pos++;
            }
        }
    }
    csr->row_ptr[n] = pos;

    return csr;
}

/**
 * Sparse matrix-vector multiplication: y = A_csr * x
 */
void csr_mat_vec_mul(const SparseMatrixCSR *csr, const double *x, double *y) {
    if (!csr || !x || !y) return;

    for (uint32_t i = 0; i < csr->n_rows; i++) {
        double sum = 0.0;
        uint32_t start = csr->row_ptr[i];
        uint32_t end   = csr->row_ptr[i + 1];
        for (uint32_t p = start; p < end; p++) {
            sum += csr->values[p] * x[csr->col_idx[p]];
        }
        y[i] = sum;
    }
}

/**
 * Free CSR matrix memory.
 */
void csr_free(SparseMatrixCSR *csr) {
    if (!csr) return;
    free(csr->values);
    free(csr->col_idx);
    free(csr->row_ptr);
    free(csr);
}

/* ==========================================================================
 * L8: Nonlinear Circuit Analysis — Newton-Raphson Method
 * ==========================================================================
 *
 * For circuits with nonlinear elements (diodes, transistors), the MNA
 * system becomes:
 *   F(x) = G(x) * x - b(x) = 0
 *
 * Newton-Raphson iteration:
 *   x_{k+1} = x_k - J^{-1}(x_k) * F(x_k)
 *
 * where J = ∂F/∂x is the Jacobian matrix.
 *
 * For a diode: I_d = I_s * (exp(V_d / (n*V_T)) - 1)
 *   Conductance: g_d = ∂I_d/∂V_d = I_s/(n*V_T) * exp(V_d/(n*V_T))
 *
 * This implementation provides the framework; specific device models
 * would stamp their contributions into the Jacobian.
 *
 * Reference: McCalla "Fundamentals of Computer-Aided Circuit Simulation" (1987)
 */
int newton_raphson_step(double *x, uint32_t n,
                        int (*eval_F)(const double *x, double *F, void *ctx),
                        int (*eval_J)(const double *x, double *J, void *ctx),
                        void *ctx, double *dx_norm) {
    if (!x || !eval_F || !eval_J || n == 0) return -1;

    /* Allocate workspace */
    double *F = (double *)malloc(n * sizeof(double));
    double *J = (double *)malloc((size_t)n * n * sizeof(double));
    double *dx = (double *)malloc(n * sizeof(double));

    if (!F || !J || !dx) {
        free(F); free(J); free(dx);
        return -1;
    }

    /* Evaluate F(x_k) */
    if (eval_F(x, F, ctx) != 0) {
        free(F); free(J); free(dx);
        return -1;
    }

    /* Evaluate Jacobian J(x_k) */
    if (eval_J(x, J, ctx) != 0) {
        free(F); free(J); free(dx);
        return -1;
    }

    /* Solve J * dx = -F */
    for (uint32_t i = 0; i < n; i++) F[i] = -F[i];

    int ret = gaussian_elimination(J, F, n, dx);
    if (ret != 0) {
        free(F); free(J); free(dx);
        return -1;
    }

    /* Compute ||dx|| */
    double norm = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        x[i] += dx[i];
        norm += dx[i] * dx[i];
    }
    *dx_norm = sqrt(norm);

    free(F);
    free(J);
    free(dx);
    return 0;
}
