#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "circuit_matrix.h"

int main(void) {
  int n = 100;
  double *A = matrix_alloc(n, n);
  double *b = (double*)calloc((size_t)n, sizeof(double));
  double *x = (double*)calloc((size_t)n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      A[i * n + j] = (i == j) ? 2.0 : ((abs(i - j) == 1) ? -1.0 : 0);
    b[i] = 1.0;
  }
  clock_t start = clock();
  matrix_solve(A, b, x, n);
  clock_t end = clock();
  printf("Solved %dx%d system in %.3f ms\n", n, n,
         (double)(end - start) * 1000 / CLOCKS_PER_SEC);
  matrix_free(A); free(b); free(x);
  return 0;
}
