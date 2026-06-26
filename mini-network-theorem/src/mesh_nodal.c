/**
 * @file mesh_nodal.c
 * @brief Nodal Analysis, Mesh Analysis, and Matrix Solvers
 *
 * Implements the two fundamental systematic circuit analysis methods:
 *   - Nodal Analysis (node-voltage method): Solve Y*V = I for node voltages
 *   - Mesh Analysis (mesh-current method): Solve Z*I = V for mesh currents
 *
 * Also provides linear algebra primitives: Gaussian elimination with partial
 * pivoting, LU decomposition (Doolittle), and forward/backward substitution.
 *
 * Knowledge Coverage:
 *   L3 - Mathematical Structures: Matrix-based circuit equations
 *   L5 - Algorithms/Methods: Gaussian elimination, LU decomposition
 *
 * Reference:
 *   Hayt, Kemmerly, Durbin "Engineering Circuit Analysis" (2019) Ch 4, 10
 *   Golub & Van Loan "Matrix Computations" (2013)
 *
 * Course Mapping:
 *   MIT 6.002: Node and mesh analysis
 *   Berkeley EE16A: Nodal analysis, matrix methods
 *   Stanford EE102A: Linear system solvers
 */

#include "network_theorem.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ==========================================================================
 * L2: Nodal Admittance Matrix Construction
 * ==========================================================================
 *
 * For a circuit with N non-ground nodes:
 *
 * Y[i][i] = sum of all admittances connected to node i
 * Y[i][j] = -sum of all admittances directly connecting node i to node j
 *
 * I[i] = algebraic sum of current sources entering node i
 *        (+ for sources directed into node, - for sources directed away)
 *
 * Voltage sources are handled via modified nodal analysis (separate module).
 * This function treats them as open for the admittance matrix (deactivated).
 */

double *build_nodal_admittance_matrix(const CircuitTopology *ckt, uint32_t *n_out) {
    if (!ckt || !n_out) return NULL;
    uint32_t n = ckt->num_nodes;
    if (n == 0) return NULL;
    *n_out = n;

    /* Allocate n×n matrix (row-major), initialize to zero */
    double *Y = (double *)calloc((size_t)n * n, sizeof(double));
    if (!Y) return NULL;

    for (uint32_t e = 0; e < ckt->num_elements; e++) {
        const CircuitElement *elem = &ckt->elements[e];
        uint32_t i = elem->node_from;
        uint32_t j = elem->node_to;

        /* Skip voltage sources in pure nodal (handled by MNA) */
        if (elem->type == ELEM_VOLTAGE_SOURCE) continue;
        if (elem->type == ELEM_OPEN) continue;

        /* Determine admittance contribution */
        double G = 0.0;
        switch (elem->type) {
            case ELEM_RESISTOR:
                if (fabs(elem->value) < 1e-15) {
                    /* Zero resistance: treat as short (infinite admittance,
                     * skip for basic nodal — needs MNA) */
                    continue;
                }
                G = 1.0 / elem->value;
                break;
            case ELEM_CAPACITOR:
                /* Admittance of capacitor: Y_C = jωC. For DC analysis,
                 * capacitor is open circuit (G=0). Skip. */
                continue;
            case ELEM_INDUCTOR:
                /* Admittance of inductor: Y_L = 1/(jωL). For DC,
                 * inductor is short (infinite G). Needs MNA. */
                continue;
            case ELEM_CURRENT_SOURCE:
                /* Current sources contribute to RHS, not Y matrix */
                continue;
            case ELEM_SHORT:
                /* Short circuit: nodes i and j are the same node.
                 * In a real implementation, we'd collapse them.
                 * Here, set a large conductance to approximate. */
                G = 1e12;
                break;
            default:
                continue;
        }

        /* Add to diagonal entries */
        if (i < n) Y[i * n + i] += G;
        if (j < n) Y[j * n + j] += G;

        /* Subtract from off-diagonal entries (mutual coupling) */
        if (i < n && j < n) {
            Y[i * n + j] -= G;
            Y[j * n + i] -= G;
        }
    }

    return Y;
}

/* ==========================================================================
 * L3: Gaussian Elimination with Partial Pivoting
 * ==========================================================================
 *
 * Solves Ax = b for x, modifying A and b in-place.
 *
 * Algorithm:
 *   For k = 0..n-2:
 *     1. Find pivot: row p with max |A[p][k]| for p >= k
 *     2. Swap rows k and p in both A and b
 *     3. For i = k+1..n-1:
 *          m = A[i][k] / A[k][k]
 *          For j = k..n-1: A[i][j] -= m * A[k][j]
 *          b[i] -= m * b[k]
 *   Back-substitution from row n-1 down to 0
 *
 * Complexity: O(n³) time, O(1) extra space (in-place)
 *
 * Returns 0 on success, -1 if singular matrix detected.
 */
int gaussian_elimination(double *A, double *b, uint32_t n, double *x) {
    if (!A || !b || !x || n == 0) return -1;

    /* Forward elimination */
    for (uint32_t k = 0; k < n; k++) {
        /* Partial pivoting: find row with maximum |A[i][k]| */
        uint32_t pivot = k;
        double max_val = fabs(A[k * n + k]);
        for (uint32_t i = k + 1; i < n; i++) {
            double val = fabs(A[i * n + k]);
            if (val > max_val) {
                max_val = val;
                pivot = i;
            }
        }

        /* Singularity check */
        if (max_val < 1e-15) return -1;

        /* Swap rows k and pivot */
        if (pivot != k) {
            for (uint32_t j = k; j < n; j++) {
                double tmp = A[k * n + j];
                A[k * n + j] = A[pivot * n + j];
                A[pivot * n + j] = tmp;
            }
            double tmp = b[k];
            b[k] = b[pivot];
            b[pivot] = tmp;
        }

        /* Eliminate below */
        double pivot_val = A[k * n + k];
        for (uint32_t i = k + 1; i < n; i++) {
            double factor = A[i * n + k] / pivot_val;
            A[i * n + k] = 0.0; /* explicitly zero for clarity */
            for (uint32_t j = k + 1; j < n; j++) {
                A[i * n + j] -= factor * A[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }

    /* Back substitution */
    for (int k = (int)n - 1; k >= 0; k--) {
        double sum = b[k];
        for (uint32_t j = k + 1; j < n; j++) {
            sum -= A[k * n + j] * x[j];
        }
        x[k] = sum / A[k * n + k];
    }

    return 0;
}

/* ==========================================================================
 * L3: Nodal Voltage Solver
 * ==========================================================================
 *
 * Constructs the RHS current vector from independent current sources
 * and solves for all node voltages.
 *
 * Current source contributions:
 *   +I at node_from (current entering)
 *   -I at node_to   (current leaving)
 */
int solve_nodal_voltages(double *Y, double *I, uint32_t n, double *V) {
    if (!Y || !I || !V || n == 0) return -1;

    /* For pure nodal analysis, we need to ground one node.
     * Assume node 0 is ground — remove its row and column.
     * This function assumes caller has already done this, or
     * the matrix is nonsingular because ground is properly set.
     *
     * Strategy: use Gaussian elimination directly.
     * If the matrix is singular (no ground), we try to fix it
     * by grounding node 0.
     */
    uint32_t sz = n;
    double *A = (double *)malloc((size_t)sz * sz * sizeof(double));
    double *b = (double *)malloc(sz * sizeof(double));
    if (!A || !b) { free(A); free(b); return -1; }

    memcpy(A, Y, (size_t)sz * sz * sizeof(double));
    memcpy(b, I, sz * sizeof(double));

    /* Ground node 0: set row 0 and col 0 to identity */
    for (uint32_t j = 0; j < sz; j++) A[0 * sz + j] = 0.0;
    for (uint32_t i = 0; i < sz; i++) A[i * sz + 0] = 0.0;
    A[0 * sz + 0] = 1.0;
    b[0] = 0.0; /* V[0] = 0 (ground) */

    int ret = gaussian_elimination(A, b, sz, V);
    free(A);
    free(b);
    return ret;
}

/* ==========================================================================
 * L3: Mesh Impedance Matrix Construction
 * ==========================================================================
 *
 * For a planar circuit with M meshes:
 *
 * Z[i][i] = sum of all impedances in mesh i (self-impedance)
 * Z[i][j] = -(impedance shared between mesh i and j) for i≠j
 *           Sign depends on relative current directions through shared branch
 *
 * V[i] = algebraic sum of voltage sources in mesh i
 *        (+ if voltage rise in mesh current direction)
 *
 * Only applicable to planar circuits with clearly defined meshes.
 */
double *build_mesh_impedance_matrix(const CircuitTopology *ckt, uint32_t *m_out) {
    if (!ckt || !m_out) return NULL;
    uint32_t m = ckt->num_meshes;
    if (m == 0) return NULL;
    *m_out = m;

    double *Z = (double *)calloc((size_t)m * m, sizeof(double));
    if (!Z) return NULL;

    /* For each mesh, sum the resistances of all branches it contains.
     * This simplified version assumes all elements are resistors.
     * Mutual impedance between meshes is the negative of shared resistance. */
    for (uint32_t mi = 0; mi < m; mi++) {
        const CircuitMesh *mesh_i = &ckt->meshes[mi];

        for (uint32_t bi = 0; bi < mesh_i->num_branches; bi++) {
            uint32_t bidx = mesh_i->branch_ids[bi];
            if (bidx >= ckt->num_branches) continue;
            const CircuitBranch *branch = &ckt->branches[bidx];
            if (branch->element_idx >= ckt->num_elements) continue;
            const CircuitElement *elem = &ckt->elements[branch->element_idx];

            double R = 0.0;
            if (elem->type == ELEM_RESISTOR) {
                R = elem->value;
            } else if (elem->type == ELEM_SHORT) {
                R = 0.0;
            } else {
                continue; /* skip non-resistive for DC mesh analysis */
            }

            /* Self-impedance: always positive contribution */
            Z[mi * m + mi] += R;

            /* Check which other meshes share this branch */
            for (uint32_t mj = mi + 1; mj < m; mj++) {
                const CircuitMesh *mesh_j = &ckt->meshes[mj];
                for (uint32_t bj = 0; bj < mesh_j->num_branches; bj++) {
                    if (mesh_j->branch_ids[bj] == bidx) {
                        /* Shared branch — subtract mutual impedance */
                        Z[mi * m + mj] -= R;
                        Z[mj * m + mi] -= R;
                    }
                }
            }
        }
    }

    return Z;
}

/* ==========================================================================
 * L3: Mesh Current Solver
 * ==========================================================================
 *
 * Constructs the RHS voltage vector and solves for mesh currents.
 * V_mesh[i] = sum of voltage source rises in mesh i
 */
int solve_mesh_currents(double *Z, double *V, uint32_t n, double *I) {
    if (!Z || !V || !I || n == 0) return -1;
    return gaussian_elimination(Z, V, n, I);
}

/* ==========================================================================
 * L5: LU Decomposition — Doolittle Algorithm
 * ==========================================================================
 *
 * Decomposes A into L*U where L is unit lower triangular and U is upper
 * triangular. Stored in-place: the upper triangle + diagonal holds U,
 * the strictly lower triangle holds L (L has 1s on diagonal, not stored).
 *
 * A[i][j] → { U[i][j] if i ≤ j
 *            { L[i][j] if i > j
 *
 * Algorithm:
 *   For k = 0..n-1:
 *     U[k][j] = A[k][j] - Σ_{m=0}^{k-1} L[k][m]*U[m][j]  for j ≥ k
 *     L[i][k] = (A[i][k] - Σ_{m=0}^{k-1} L[i][m]*U[m][k]) / U[k][k]  for i > k
 *
 * Complexity: O(n³)
 * Returns: 0 on success, -1 if zero pivot (singular matrix)
 */
int lu_decompose(double *A, uint32_t n) {
    if (!A || n == 0) return -1;

    for (uint32_t k = 0; k < n; k++) {
        /* Compute U[k][j] for j >= k */
        for (uint32_t j = k; j < n; j++) {
            double sum = 0.0;
            for (uint32_t m = 0; m < k; m++) {
                sum += A[k * n + m] * A[m * n + j];
            }
            A[k * n + j] -= sum;
        }

        /* Pivot check */
        if (fabs(A[k * n + k]) < 1e-15) return -1;

        /* Compute L[i][k] for i > k */
        for (uint32_t i = k + 1; i < n; i++) {
            double sum = 0.0;
            for (uint32_t m = 0; m < k; m++) {
                sum += A[i * n + m] * A[m * n + k];
            }
            A[i * n + k] = (A[i * n + k] - sum) / A[k * n + k];
        }
    }
    return 0;
}

/* ==========================================================================
 * L5: LU Solve — Forward and Backward Substitution
 * ==========================================================================
 *
 * Given LU decomposition of A (stored in-place in LU), solves Ax = b.
 *
 * Step 1 — Forward substitution (Ly = b):
 *   y[0] = b[0]
 *   y[i] = b[i] - Σ_{j=0}^{i-1} L[i][j] * y[j]
 *
 * Step 2 — Backward substitution (Ux = y):
 *   x[n-1] = y[n-1] / U[n-1][n-1]
 *   x[i] = (y[i] - Σ_{j=i+1}^{n-1} U[i][j] * x[j]) / U[i][i]
 */
void lu_solve(const double *LU, const double *b, uint32_t n, double *x) {
    if (!LU || !b || !x || n == 0) return;

    /* Allocate temporary y vector */
    double *y = (double *)malloc(n * sizeof(double));
    if (!y) return;

    /* Forward substitution: Ly = b */
    for (uint32_t i = 0; i < n; i++) {
        double sum = b[i];
        for (uint32_t j = 0; j < i; j++) {
            sum -= LU[i * n + j] * y[j];
        }
        y[i] = sum; /* L[i][i] = 1, no division needed */
    }

    /* Backward substitution: Ux = y */
    for (int i = (int)n - 1; i >= 0; i--) {
        double sum = y[i];
        for (uint32_t j = i + 1; j < n; j++) {
            sum -= LU[(uint32_t)i * n + j] * x[j];
        }
        x[i] = sum / LU[(uint32_t)i * n + (uint32_t)i];
    }

    free(y);
}

/* ==========================================================================
 * L5: Matrix-Vector Multiply — y = A*x
 * ==========================================================================
 *
 * Computes y[i] = Σ_{j=0}^{n-1} A[i*n + j] * x[j]
 * Used for residual computation and iterative refinement.
 */
static void mat_vec_mul(const double *A, const double *x, uint32_t n, double *y) {
    for (uint32_t i = 0; i < n; i++) {
        double sum = 0.0;
        for (uint32_t j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        y[i] = sum;
    }
}

/* ==========================================================================
 * L5: Compute Residual Norm — ||b - A*x||₂
 * ==========================================================================
 *
 * Residual r = b - A*x. Returns the Euclidean norm.
 * Used to assess solution quality from Gaussian elimination.
 */
double compute_residual_norm(const double *A, const double *x,
                             const double *b, uint32_t n) {
    if (!A || !x || !b || n == 0) return -1.0;

    double *Ax = (double *)malloc(n * sizeof(double));
    if (!Ax) return -1.0;

    mat_vec_mul(A, x, n, Ax);

    double norm = 0.0;
    for (uint32_t i = 0; i < n; i++) {
        double diff = b[i] - Ax[i];
        norm += diff * diff;
    }
    free(Ax);
    return sqrt(norm);
}

/* ==========================================================================
 * L5: Matrix Condition Number Estimate (1-norm)
 * ==========================================================================
 *
 * Estimates the condition number κ₁(A) = ||A||₁ * ||A⁻¹||₁.
 * A well-conditioned matrix has κ ≈ 1-100; poorly conditioned has κ >> 1000.
 *
 * Uses the Hager-Higham estimator for ||A⁻¹||₁ without explicit inversion.
 * This is a simplified version using power iteration on A⁻¹.
 */
double matrix_condition_estimate(const double *A, uint32_t n) {
    if (!A || n == 0) return -1.0;

    /* Compute ||A||₁ (max column sum) */
    double norm_A = 0.0;
    for (uint32_t j = 0; j < n; j++) {
        double col_sum = 0.0;
        for (uint32_t i = 0; i < n; i++) {
            col_sum += fabs(A[i * n + j]);
        }
        if (col_sum > norm_A) norm_A = col_sum;
    }

    /* Estimate ||A⁻¹||₁ via 5 iterations of power method on A⁻¹.
     * This is a crude estimate but sufficient for diagnostic purposes. */
    double *y = (double *)malloc(n * sizeof(double));
    double *z = (double *)malloc(n * sizeof(double));
    if (!y || !z) { free(y); free(z); return -1.0; }

    /* Initial vector: all ones */
    for (uint32_t i = 0; i < n; i++) y[i] = 1.0;

    double norm_Ainv = 0.0;
    int max_iter = 5;
    for (int iter = 0; iter < max_iter; iter++) {
        /* Solve A*z = y using LU */
        double *LU = (double *)malloc((size_t)n * n * sizeof(double));
        if (!LU) { free(y); free(z); return -1.0; }
        memcpy(LU, A, (size_t)n * n * sizeof(double));

        if (lu_decompose(LU, n) != 0) {
            free(LU);
            break;
        }
        lu_solve(LU, y, n, z);
        free(LU);

        /* ||z||₁ */
        double z_norm = 0.0;
        for (uint32_t i = 0; i < n; i++) z_norm += fabs(z[i]);

        if (z_norm > norm_Ainv) norm_Ainv = z_norm;

        /* Update y: sign pattern for next iteration */
        for (uint32_t i = 0; i < n; i++) {
            y[i] = (z[i] >= 0) ? 1.0 : -1.0;
        }
    }

    free(y);
    free(z);

    return (norm_A > 1e-15) ? norm_A * norm_Ainv : -1.0;
}

/* ==========================================================================
 * L3: Matrix Determinant via LU Decomposition
 * ==========================================================================
 *
 * det(A) = product of U's diagonal elements (since det(L)=1).
 * Sign must account for row swaps if pivoting is used.
 *
 * Here we use plain LU (no pivoting), so det(A) = Π U[i][i].
 */
double matrix_determinant(const double *A, uint32_t n) {
    if (!A || n == 0) return 0.0;

    double *LU = (double *)malloc((size_t)n * n * sizeof(double));
    if (!LU) return 0.0;
    memcpy(LU, A, (size_t)n * n * sizeof(double));

    if (lu_decompose(LU, n) != 0) {
        free(LU);
        return 0.0; /* singular → det = 0 */
    }

    double det = 1.0;
    for (uint32_t i = 0; i < n; i++) {
        det *= LU[i * n + i];
    }
    free(LU);
    return det;
}

/* ==========================================================================
 * L3: Matrix Inverse via Gaussian Elimination
 * ==========================================================================
 *
 * Computes A⁻¹ by solving A*X = I column by column.
 * inv_A is stored in row-major order with n×n elements.
 *
 * Complexity: O(n³) — solves n linear systems
 */
int matrix_inverse(const double *A, uint32_t n, double *inv_A) {
    if (!A || !inv_A || n == 0) return -1;

    /* Create identity matrix for RHS */
    double *I = (double *)calloc((size_t)n * n, sizeof(double));
    if (!I) return -1;
    for (uint32_t i = 0; i < n; i++) I[i * n + i] = 1.0;

    /* Solve for each column */
    for (uint32_t col = 0; col < n; col++) {
        double *A_copy = (double *)malloc((size_t)n * n * sizeof(double));
        double *b_copy = (double *)malloc(n * sizeof(double));
        if (!A_copy || !b_copy) {
            free(A_copy); free(b_copy); free(I); return -1;
        }
        memcpy(A_copy, A, (size_t)n * n * sizeof(double));
        for (uint32_t i = 0; i < n; i++) b_copy[i] = I[i * n + col];

        double *x = (double *)malloc(n * sizeof(double));
        if (!x) {
            free(A_copy); free(b_copy); free(I); return -1;
        }

        int ret = gaussian_elimination(A_copy, b_copy, n, x);
        if (ret != 0) {
            free(A_copy); free(b_copy); free(x); free(I);
            return -1;
        }

        for (uint32_t i = 0; i < n; i++) {
            inv_A[i * n + col] = x[i];
        }

        free(A_copy); free(b_copy); free(x);
    }

    free(I);
    return 0;
}

/* ==========================================================================
 * L2: Build RHS Current Vector for Nodal Analysis
 * ==========================================================================
 *
 * I[i] = sum of independent current sources entering node i
 *        minus those leaving node i.
 *
 * Also handles current from voltage sources that have been source-transformed.
 */
int build_current_vector(const CircuitTopology *ckt, double *I, uint32_t n) {
    if (!ckt || !I || n == 0) return -1;

    memset(I, 0, n * sizeof(double));

    for (uint32_t e = 0; e < ckt->num_elements; e++) {
        const CircuitElement *elem = &ckt->elements[e];
        if (elem->type != ELEM_CURRENT_SOURCE) continue;

        uint32_t from = elem->node_from;
        uint32_t to   = elem->node_to;
        double Is = elem->value;

        /* Current source directed from node 'from' to node 'to':
         * Current enters node 'to', leaves node 'from' */
        if (from < n) I[from] -= Is;
        if (to < n)   I[to]   += Is;
    }
    return 0;
}
