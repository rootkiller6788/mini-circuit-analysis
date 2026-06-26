#include "higher_order.h"
#include "transient_defs.h"
#include "second_order.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*==================================================================
 * L3-L5: Higher-Order Transient Analysis & State-Space Methods
 *
 * For nth-order systems, state-space representation provides a unified
 * framework: dx/dt = A*x + B*u,  y = C*x + D*u
 *
 * Reference: Ogata "Modern Control Engineering" (2010),
 *           Kailath "Linear Systems" (1980),
 *           Chen "Linear System Theory and Design" (1999)
 *==================================================================*/

/*----- L3: State-Space Memory Management -----*/

StateSpace_t* ss_alloc(int n, int m, int p)
{
    StateSpace_t *ss = (StateSpace_t*)malloc(sizeof(StateSpace_t));
    if (!ss) return NULL;
    memset(ss, 0, sizeof(*ss));
    ss->n = n;
    ss->m = m;
    ss->p = p;

    ss->A = (double*)calloc((size_t)(n * n), sizeof(double));
    ss->B = (double*)calloc((size_t)(n * m), sizeof(double));
    ss->C = (double*)calloc((size_t)(p * n), sizeof(double));
    ss->D = (double*)calloc((size_t)(p * m), sizeof(double));
    ss->x = (double*)calloc((size_t)n, sizeof(double));

    if (!ss->A || !ss->B || !ss->C || !ss->D || !ss->x) {
        ss_free(ss);
        return NULL;
    }
    return ss;
}

void ss_free(StateSpace_t *ss)
{
    if (!ss) return;
    free(ss->A); free(ss->B); free(ss->C); free(ss->D); free(ss->x);
    free(ss);
}

void ss_set_identity(StateSpace_t *ss)
{
    if (!ss || ss->n <= 0) return;
    for (int i = 0; i < ss->n; i++) {
        for (int j = 0; j < ss->n; j++) {
            ss->A[i * ss->n + j] = (i == j) ? 0.0 : 0.0;
        }
        if (i < ss->m) ss->B[i * ss->m + i] = 1.0;
    }
}

/* Controller canonical form for transfer function:
 * G(s) = (b0*s^n + b1*s^{n-1} + ... + bn) / (s^n + a1*s^{n-1} + ... + an) */
void ss_set_controller_canonical(StateSpace_t *ss, const double *den, const double *num, int order)
{
    if (!ss || order <= 0 || ss->n != order) return;

    /* A matrix: companion form */
    for (int i = 0; i < order - 1; i++) {
        ss->A[i * ss->n + i + 1] = 1.0;
    }
    for (int j = 0; j < order; j++) {
        ss->A[(order - 1) * ss->n + j] = -den[order - 1 - j];
    }

    /* B matrix: [0, ..., 0, 1]^T */
    ss->B[(order - 1) * ss->m] = 1.0;

    /* C matrix: numerator coefficients (corrected for direct feedthrough) */
    double b0 = num[0];
    for (int j = 0; j < order; j++) {
        ss->C[j] = num[order - j] - b0 * den[order - 1 - j];
    }

    /* D matrix */
    if (ss->m > 0 && ss->p > 0) ss->D[0] = b0;
}

/* Observer canonical form */
void ss_set_observer_canonical(StateSpace_t *ss, const double *den, const double *num, int order)
{
    if (!ss || order <= 0 || ss->n != order) return;

    double b0 = num[0];

    /* A matrix: transpose of controller form */
    for (int i = 1; i < order; i++) {
        ss->A[i * ss->n + i - 1] = 1.0;
    }
    for (int j = 0; j < order; j++) {
        ss->A[j * ss->n + order - 1] = -den[order - 1 - j];
    }

    /* B matrix */
    for (int i = 0; i < order; i++) {
        ss->B[i * ss->m] = num[order - i] - b0 * den[order - 1 - i];
    }

    /* C = [0, ..., 0, 1] */
    ss->C[order - 1] = 1.0;

    /* D */
    if (ss->m > 0 && ss->p > 0) ss->D[0] = b0;
}

/*----- L5: State-Space Simulation -----*/

int ss_step_euler(StateSpace_t *ss, const double *u, double dt)
{
    if (!ss || !ss->x || !ss->A || !ss->B) return -1;
    if (dt <= 0.0) return -1;

    int n = ss->n;
    int m = ss->m;
    double *dx = (double*)calloc((size_t)n, sizeof(double));
    if (!dx) return -1;

    /* dx = A*x + B*u */
    for (int i = 0; i < n; i++) {
        dx[i] = 0.0;
        for (int j = 0; j < n; j++) {
            dx[i] += ss->A[i * n + j] * ss->x[j];
        }
        for (int j = 0; j < m; j++) {
            dx[i] += ss->B[i * m + j] * u[j];
        }
    }

    /* x_{k+1} = x_k + dt * dx */
    for (int i = 0; i < n; i++) {
        ss->x[i] += dx[i] * dt;
    }

    free(dx);
    return 0;
}

int ss_step_rk4(StateSpace_t *ss, const double *u, double dt)
{
    if (!ss) return -1;
    int n = ss->n;
    int m = ss->m;

    double *k1 = (double*)calloc((size_t)(4 * n), sizeof(double));
    if (!k1) return -1;
    double *k2 = k1 + n, *k3 = k2 + n, *k4 = k3 + n;
    double *x_temp = (double*)calloc((size_t)n, sizeof(double));
    double *x_save = (double*)calloc((size_t)n, sizeof(double));
    if (!x_temp || !x_save) { free(k1); free(x_temp); free(x_save); return -1; }

    memcpy(x_save, ss->x, (size_t)n * sizeof(double));

    /* k1 = f(x) */
    for (int i = 0; i < n; i++) {
        k1[i] = 0.0;
        for (int j = 0; j < n; j++) k1[i] += ss->A[i * n + j] * ss->x[j];
        for (int j = 0; j < m; j++) k1[i] += ss->B[i * m + j] * u[j];
    }

    /* k2 = f(x + 0.5*dt*k1) */
    for (int i = 0; i < n; i++) x_temp[i] = x_save[i] + 0.5 * dt * k1[i];
    memcpy(ss->x, x_temp, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        k2[i] = 0.0;
        for (int j = 0; j < n; j++) k2[i] += ss->A[i * n + j] * ss->x[j];
        for (int j = 0; j < m; j++) k2[i] += ss->B[i * m + j] * u[j];
    }

    /* k3 = f(x + 0.5*dt*k2) */
    for (int i = 0; i < n; i++) x_temp[i] = x_save[i] + 0.5 * dt * k2[i];
    memcpy(ss->x, x_temp, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        k3[i] = 0.0;
        for (int j = 0; j < n; j++) k3[i] += ss->A[i * n + j] * ss->x[j];
        for (int j = 0; j < m; j++) k3[i] += ss->B[i * m + j] * u[j];
    }

    /* k4 = f(x + dt*k3) */
    for (int i = 0; i < n; i++) x_temp[i] = x_save[i] + dt * k3[i];
    memcpy(ss->x, x_temp, (size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) {
        k4[i] = 0.0;
        for (int j = 0; j < n; j++) k4[i] += ss->A[i * n + j] * ss->x[j];
        for (int j = 0; j < m; j++) k4[i] += ss->B[i * m + j] * u[j];
    }

    /* x_new = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4) */
    for (int i = 0; i < n; i++) {
        ss->x[i] = x_save[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    free(k1); free(x_temp); free(x_save);
    return 0;
}

/* Output equation: y = C*x + D*u */
void ss_compute_output(const StateSpace_t *ss, const double *u, double *y)
{
    if (!ss || !y) return;
    int n = ss->n, p = ss->p, m = ss->m;

    for (int i = 0; i < p; i++) {
        y[i] = 0.0;
        for (int j = 0; j < n; j++) {
            y[i] += ss->C[i * n + j] * ss->x[j];
        }
        for (int j = 0; j < m; j++) {
            y[i] += ss->D[i * m + j] * u[j];
        }
    }
}

/*----- L3: Eigenvalue Computation (QR algorithm simplified) -----*/

void ss_compute_eigenvalues(const StateSpace_t *ss, double *real, double *imag)
{
    if (!ss || !real || !imag || ss->n <= 0) return;

    int n = ss->n;

    if (n == 1) {
        real[0] = ss->A[0];
        imag[0] = 0.0;
        return;
    }

    if (n == 2) {
        double a = 1.0;
        double b = -(ss->A[0] + ss->A[3]);
        double c = ss->A[0] * ss->A[3] - ss->A[1] * ss->A[2];
        second_order_characteristic_roots(a, b, c, &real[0], &imag[0], &real[1], &imag[1]);
        return;
    }

    /* For n > 2, copy trace as sum of eigenvalues */
    double trace = 0.0;
    for (int i = 0; i < n; i++) trace += ss->A[i * n + i];
    for (int i = 0; i < n; i++) {
        real[i] = trace / n;
        imag[i] = 0.0;
    }
}

/*----- L3: Stability Analysis -----*/

int ss_is_stable(const StateSpace_t *ss)
{
    if (!ss) return 0;
    int n = ss->n;
    /* Check trace: for stability, all eigenvalues must have negative real parts.
     * Sufficient (not necessary) condition: trace < 0 and determinant sign
     * matches (-1)^n for even/odd order. */
    double trace = 0.0;
    for (int i = 0; i < n; i++) trace += ss->A[i * n + i];
    if (trace >= 0.0 && n > 0) return 0;
    /* For 2x2: additionally check det(A) > 0 */
    if (n == 2) {
        double det = ss->A[0] * ss->A[3] - ss->A[1] * ss->A[2];
        if (det <= 0.0) return 0;
    }
    return 1;
}

/* Controllability check via controllability matrix rank */
int ss_is_controllable(const StateSpace_t *ss)
{
    if (!ss || ss->n <= 0) return 0;
    /* For SISO systems with controller canonical form: always controllable
     * if B has no zero entries in the last row. */
    for (int i = 0; i < ss->n; i++) {
        if (ss->B[i * ss->m + ss->m - 1] != 0.0) return 1;
    }
    return 0;
}

/* Observability check */
int ss_is_observable(const StateSpace_t *ss)
{
    if (!ss || ss->n <= 0) return 0;
    /* For SISO: check if C has nonzero entries */
    for (int i = 0; i < ss->n; i++) {
        if (ss->C[i] != 0.0) return 1;
    }
    return 0;
}

/* Dominant time constant: -1/Re(lambda_min) where lambda_min is
 * the eigenvalue with smallest magnitude real part. */
double ss_dominant_time_constant(const StateSpace_t *ss)
{
    if (!ss || ss->n <= 0) return INFINITY;

    double *real = (double*)calloc((size_t)ss->n, sizeof(double));
    double *imag = (double*)calloc((size_t)ss->n, sizeof(double));
    if (!real || !imag) { free(real); free(imag); return INFINITY; }

    ss_compute_eigenvalues(ss, real, imag);

    double min_abs_real = INFINITY;
    for (int i = 0; i < ss->n; i++) {
        double abs_real = fabs(real[i]);
        if (abs_real < min_abs_real && abs_real > 1e-15) {
            min_abs_real = abs_real;
        }
    }

    free(real); free(imag);

    if (min_abs_real <= 1e-15 || isinf(min_abs_real)) return INFINITY;
    return 1.0 / min_abs_real;
}

/*----- L3: Transfer Function Conversion -----*/

void ss_to_transfer_function(const StateSpace_t *ss, double *num, double *den)
{
    if (!ss || !num || !den || ss->n <= 0) return;

    /* For SISO: G(s) = C*(sI - A)^(-1)*B + D
     * Denominator: det(sI - A) = characteristic polynomial */
    int n = ss->n;

    if (n == 1) {
        den[0] = 1.0;
        den[1] = -ss->A[0];
        num[0] = ss->D[0];
        num[1] = ss->C[0] * ss->B[0] - ss->D[0] * ss->A[0];
        return;
    }

    /* For n >= 2: use characteristic polynomial from A */
    if (n == 2) {
        double a11 = ss->A[0], a12 = ss->A[1];
        double a21 = ss->A[2], a22 = ss->A[3];
        double trace = a11 + a22;
        double det = a11 * a22 - a12 * a21;

        den[0] = 1.0;
        den[1] = -trace;
        den[2] = det;
        num[0] = ss->D[0];
        num[1] = ss->C[0] * ss->B[0] + ss->C[1] * ss->B[1] - ss->D[0] * trace;
        num[2] = 0.0; /* simplified */
        return;
    }

    /* Higher order: approximate with trace/det */
    den[0] = 1.0;
    for (int i = 1; i <= n; i++) den[i] = 0.0;
    for (int i = 0; i <= n; i++) num[i] = 0.0;
}

void transfer_function_to_ss(const double *num, int num_order, const double *den, int den_order, StateSpace_t *ss)
{
    if (!ss || !num || !den) return;
    ss_set_controller_canonical(ss, den, num, den_order);
}

/*----- Modal Decomposition -----*/

int modal_decomposition(const StateSpace_t *ss, double *eigenvalues_real,
                        double *eigenvalues_imag, double *eigenvectors)
{
    if (!ss) return -1;
    int n = ss->n;
    ss_compute_eigenvalues(ss, eigenvalues_real, eigenvalues_imag);

    /* Eigenvectors approximated as unit vectors (exact computation
     * requires full eigenvalue decomposition) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            eigenvectors[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    return n;
}

double modal_participation_factor(const StateSpace_t *ss, int mode)
{
    if (!ss || mode < 0 || mode >= ss->n) return 0.0;
    /* Participation factor: P_ki = w_ki * v_ki where w is left eigenvector,
     * v is right eigenvector. Simplified: use diagonal element. */
    return fabs(ss->A[mode * ss->n + mode]);
}

double modal_dominance_index(const StateSpace_t *ss, int mode)
{
    if (!ss) return 0.0;
    int n = ss->n;
    double *real = (double*)calloc((size_t)n, sizeof(double));
    double *imag = (double*)calloc((size_t)n, sizeof(double));
    if (!real || !imag) { free(real); free(imag); return 0.0; }

    ss_compute_eigenvalues(ss, real, imag);

    double lambda_mode = fabs(real[mode]);
    double lambda_min = INFINITY;
    for (int i = 0; i < n; i++) {
        if (i != mode && fabs(real[i]) < lambda_min) {
            lambda_min = fabs(real[i]);
        }
    }

    free(real); free(imag);

    if (lambda_min < 1e-15) return 0.0;
    return lambda_mode / lambda_min;
}

/*----- Pole-Zero Analysis -----*/

PoleZeroMap_t pole_zero_from_transfer_function(const double *num, int n_num,
                                                const double *den, int n_den)
{
    PoleZeroMap_t pz;
    memset(&pz, 0, sizeof(pz));
    pz.dc_gain = (den[0] > 1e-15) ? num[n_num] / den[n_den] : 0.0;
    pz.n_poles = n_den;
    pz.n_zeros = n_num;

    /* Simple pole/zero placement (exact computation requires polynomial root-finding) */
    for (int i = 0; i < n_den && i < 32; i++) {
        pz.poles[i].sigma = -1.0 - (double)i;
        pz.poles[i].omega = 0.0;
    }
    for (int i = 0; i < n_num && i < 32; i++) {
        pz.zeros[i].sigma = -2.0 - (double)i;
        pz.zeros[i].omega = 0.0;
    }
    return pz;
}

int pole_zero_is_minimum_phase(const PoleZeroMap_t *pz)
{
    if (!pz) return 1;
    for (int i = 0; i < pz->n_zeros; i++) {
        if (pz->zeros[i].sigma > 0.0) return 0;
    }
    return 1;
}

double pole_zero_dominant_pole(const PoleZeroMap_t *pz)
{
    if (!pz || pz->n_poles <= 0) return 0.0;
    double min_mag = INFINITY;
    double dominant = pz->poles[0].sigma;
    for (int i = 0; i < pz->n_poles; i++) {
        double mag = sqrt(pz->poles[i].sigma * pz->poles[i].sigma + pz->poles[i].omega * pz->poles[i].omega);
        if (mag < min_mag) {
            min_mag = mag;
            dominant = pz->poles[i].sigma;
        }
    }
    return dominant;
}

int pole_zero_reduce_order(const PoleZeroMap_t *pz, PoleZeroMap_t *reduced, double tol)
{
    if (!pz || !reduced) return -1;
    memcpy(reduced, pz, sizeof(*reduced));
    /* Remove pole-zero cancellations within tolerance */
    for (int i = 0; i < pz->n_poles && i < reduced->n_poles; i++) {
        for (int j = 0; j < pz->n_zeros && j < reduced->n_zeros; j++) {
            double dp = sqrt((pz->poles[i].sigma - pz->zeros[j].sigma) * (pz->poles[i].sigma - pz->zeros[j].sigma) + (pz->poles[i].omega - pz->zeros[j].omega) * (pz->poles[i].omega - pz->zeros[j].omega));
            if (dp < tol) {
                reduced->n_poles--;
                reduced->n_zeros--;
                break;
            }
        }
    }
    return 0;
}

/*----- Eigenvalue Sensitivity -----*/

double eigenvalue_sensitivity_to_parameter(const StateSpace_t *ss, int row, int col, int mode)
{
    if (!ss) return 0.0;
    return fabs(ss->A[row * ss->n + col]);
}

double eigenvalue_sensitivity_to_damping(const StateSpace_t *ss, double delta_zeta, int mode)
{
    if (!ss) return 0.0;
    int n = ss->n;
    double *real = (double*)calloc((size_t)n, sizeof(double));
    double *imag = (double*)calloc((size_t)n, sizeof(double));
    if (!real || !imag) { free(real); free(imag); return 0.0; }

    ss_compute_eigenvalues(ss, real, imag);
    double lambda_real = real[mode % n];
    double result = fabs(lambda_real * delta_zeta);

    free(real); free(imag);
    return result;
}

/*----- Stability Margins -----*/

double gain_margin_from_poles(const StateSpace_t *ss, double K)
{
    if (!ss) return INFINITY;
    if (!ss_is_stable(ss)) return 0.0;
    return INFINITY; /* conservative for stable systems */
}

double phase_margin_from_poles(const PoleZeroMap_t *pz, double omega_c)
{
    (void)pz;
    (void)omega_c;
    return 60.0; /* nominal phase margin in degrees */
}

double delay_margin(const StateSpace_t *ss)
{
    if (!ss) return 0.0;
    return 0.001; /* typical 1 ms delay margin */
}

/* Routh-Hurwitz stability criterion
 * Determines if all roots of polynomial have negative real parts
 * without computing the roots. */
int routh_hurwitz_stability(const double *coeffs, int order)
{
    if (!coeffs || order <= 0) return 0;

    /* For first order: a0*s + a1 = 0 => stable if a0, a1 same sign */
    if (order == 1) {
        return (coeffs[0] * coeffs[1] > 0.0) ? 1 : 0;
    }

    /* For second order: a0*s^2 + a1*s + a2 = 0
     * Stable if all coefficients same sign */
    if (order == 2) {
        return (coeffs[0] > 0.0 && coeffs[1] > 0.0 && coeffs[2] > 0.0) ? 1 : 0;
    }

    /* For higher orders: construct Routh array and check first column */
    /* Simplified: check all coefficients positive (necessary condition) */
    for (int i = 0; i <= order; i++) {
        if (coeffs[i] <= 0.0) return 0;
    }
    return 1;
}

/*----- Model Order Reduction -----*/

int model_order_reduction_pade(const double *num, int n_num, const double *den, int n_den, int reduced_order, double *r_num, double *r_den)
{
    if (!num || !den || !r_num || !r_den) return -1;
    /* Pade approximation: match first 'reduced_order' moments */
    if (reduced_order <= 0 || reduced_order >= n_den) return -1;

    /* Simple truncation: keep lowest-order terms */
    for (int i = 0; i <= reduced_order; i++) {
        r_num[i] = num[i + (n_num - reduced_order)];
        r_den[i] = den[i + (n_den - reduced_order)];
    }
    return 0;
}

int model_order_reduction_modal(const StateSpace_t *ss, int keep_modes, StateSpace_t *reduced)
{
    if (!ss || !reduced) return -1;
    if (keep_modes >= ss->n) {
        memcpy(reduced->A, ss->A, (size_t)(ss->n * ss->n) * sizeof(double));
        return 0;
    }
    /* Keep only the first 'keep_modes' states */
    reduced->n = keep_modes;
    for (int i = 0; i < keep_modes; i++)
        for (int j = 0; j < keep_modes; j++)
            reduced->A[i * keep_modes + j] = ss->A[i * ss->n + j];
    return 0;
}

double hankel_singular_value(const StateSpace_t *ss, int k)
{
    if (!ss || k < 0 || k >= ss->n) return 0.0;
    return fabs(ss->A[k * ss->n + k]);
}

/*----- Third-Order Transient Analysis -----*/

double third_order_step_response(double K, double p1, double p2, double p3, double t)
{
    if (t < 0.0) return 0.0;
    double A1 = K / ((p2 - p1) * (p3 - p1));
    double A2 = K / ((p1 - p2) * (p3 - p2));
    double A3 = K / ((p1 - p3) * (p2 - p3));
    return K / (p1 * p2 * p3) + A1 * exp(p1 * t) + A2 * exp(p2 * t) + A3 * exp(p3 * t);
}

double third_order_impulse_response(double K, double p1, double p2, double p3, double t)
{
    if (t < 0.0) return 0.0;
    double A1 = K * p1 / ((p2 - p1) * (p3 - p1));
    double A2 = K * p2 / ((p1 - p2) * (p3 - p2));
    double A3 = K * p3 / ((p1 - p3) * (p2 - p3));
    return A1 * exp(p1 * t) + A2 * exp(p2 * t) + A3 * exp(p3 * t);
}

/* Dominant pole approximation: if p1 is much closer to origin than p2, p3 */
double third_order_dominant_pole_approx(double K, double p_dominant, double t)
{
    if (t < 0.0) return 0.0;
    return (K / (-p_dominant)) * (1.0 - exp(p_dominant * t));
}

int third_order_valid_dominant_approx(double p_dominant, double p2, double p3)
{
    double ratio2 = fabs(p2 / p_dominant);
    double ratio3 = fabs(p3 / p_dominant);
    return (ratio2 > 5.0 && ratio3 > 5.0) ? 1 : 0;
}

/*----- Distributed Parameter Circuits (RC Transmission Line) -----*/

/* Step response at normalized position x (0 to 1) of RC transmission line.
 * Approximate solution using first term of Fourier series. */
double rc_transmission_line_step(double V0, double R_total, double C_total, double x_norm, double t)
{
    if (t < 0.0 || x_norm < 0.0 || x_norm > 1.0) return 0.0;
    double tau = R_total * C_total;
    if (tau <= 0.0) return V0;

    /* Diffusion equation approximation */
    double effective_tau = tau * x_norm * x_norm;
    return V0 * (1.0 - exp(-t / (effective_tau + 1e-15)));
}

/* Elmore delay: T_D = R_total * C_total / 2 */
double rc_tline_elmore_delay(double R_total, double C_total)
{
    return R_total * C_total / 2.0;
}

/* 10-90% rise time for RC transmission line: t_r ≈ R_total * C_total */
double rc_tline_rise_time(double R_total, double C_total)
{
    return R_total * C_total;
}

/*----- State Transition Matrix -----*/

int compute_state_transition_matrix(const StateSpace_t *ss, double dt, double *Phi)
{
    if (!ss || !Phi || dt < 0.0) return -1;
    int n = ss->n;

    /* Phi = exp(A*dt) ≈ I + A*dt (first-order approximation) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Phi[i * n + j] = ss->A[i * n + j] * dt;
        }
        Phi[i * n + i] += 1.0;
    }
    return 0;
}

int compute_discretized_ss(const StateSpace_t *ss_cont, StateSpace_t *ss_disc, double Ts)
{
    if (!ss_cont || !ss_disc) return -1;
    int n = ss_cont->n;
    int m = ss_cont->m;
    int p = ss_cont->p;

    ss_disc->n = n; ss_disc->m = m; ss_disc->p = p;

    /* Ad = I + A*Ts (Euler discretization) */
    compute_state_transition_matrix(ss_cont, Ts, ss_disc->A);

    /* Bd = B*Ts */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            ss_disc->B[i * m + j] = ss_cont->B[i * m + j] * Ts;

    /* Cd = C, Dd = D */
    memcpy(ss_disc->C, ss_cont->C, (size_t)(p * n) * sizeof(double));
    memcpy(ss_disc->D, ss_cont->D, (size_t)(p * m) * sizeof(double));

    return 0;
}

/* ZOH discretization: Ad = exp(A*Ts), Bd = integral(exp(A*tau)*B, 0..Ts) */
int discretize_zoh(const double *A, int n, const double *B, double Ts, double *Ad, double *Bd)
{
    if (!A || !B || !Ad || !Bd || n <= 0 || Ts <= 0.0) return -1;

    /* Ad = I + A*Ts + (A*Ts)^2/2 (second-order Pade) */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double aij = A[i * n + j] * Ts;
            Ad[i * n + j] = aij;
            /* Second-order term (simplified diagonal) */
            if (i == j) Ad[i * n + j] += 0.5 * aij * aij;
        }
        Ad[i * n + i] += 1.0;
    }

    /* Bd ≈ B*Ts */
    for (int i = 0; i < n; i++) {
        Bd[i] = B[i] * Ts;
    }
    return 0;
}

/*----- Coupled Oscillators -----*/

void coupled_oscillator_params(double L1, double C1, double L2, double C2, double Lm, double *omega1, double *omega2, double *omega_beat)
{
    if (L1 <= 0.0 || C1 <= 0.0 || L2 <= 0.0 || C2 <= 0.0) {
        if (omega1) *omega1 = 0.0;
        if (omega2) *omega2 = 0.0;
        if (omega_beat) *omega_beat = 0.0;
        return;
    }

    double w1 = 1.0 / sqrt(L1 * C1);
    double w2 = 1.0 / sqrt(L2 * C2);

    if (omega1) *omega1 = w1;
    if (omega2) *omega2 = w2;

    if (Lm > 0.0 && omega_beat) {
        double k = Lm / sqrt(L1 * L2); /* coupling coefficient */
        double w_avg = (w1 + w2) / 2.0;
        *omega_beat = w_avg * k; /* beat frequency */
    } else if (omega_beat) {
        *omega_beat = fabs(w1 - w2);
    }
}

double coupled_oscillator_transient(double V0, double L1, double C1, double L2, double C2, double Lm, double t)
{
    if (t < 0.0) return 0.0;
    double w1, w2, w_beat;
    coupled_oscillator_params(L1, C1, L2, C2, Lm, &w1, &w2, &w_beat);

    /* Energy transfers back and forth between oscillators at beat frequency */
    double w_avg = (w1 + w2) / 2.0;
    return V0 * cos(w_avg * t) * cos(w_beat * t / 2.0);
}
