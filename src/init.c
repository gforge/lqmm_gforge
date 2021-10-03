#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */


extern void C_gradientSh(double theta[], double *x, double *y, double *z, double *weights, double *Tq, double *V, double *W, double *sigma, float *quantile, int *p, int *q, int *m, int *M, int *N, int *Kq, int *start, int *end, double *STEP, double *beta, double *gamma, int *RESET_STEP, double *TOL_LL, double *TOL_THETA, int *CHECK_THETA, int *MAXIT, int *verbose, int CONVERGE[], double loglik[], double deriv[], double optimum[]);
extern void C_gradientSi(double theta[], double *x, double *y, float *quantile, int *N, int *p, double *STEP, double *beta, double *gamma, int *RESET_STEP, double *TOL_LL, double *TOL_THETA, int *CHECK_THETA, int *MAXIT, int *verbose, int CONVERGE[], double deriv[], double optimum[]);
extern void C_ll_h(double theta[], double *x, double *y, double *z, double *weights, double *Tq, double *V, double *W, double *sigma, float *quantile, int *p, int *q, int *m, int *M, int *N, int *Kq, int *start, int *end, double *loglik);

static const R_CMethodDef CEntries[] = {
    {"C_gradientSh", (DL_FUNC) &C_gradientSh, 31},
    {"C_gradientSi", (DL_FUNC) &C_gradientSi, 18},
    {"C_ll_h",       (DL_FUNC) &C_ll_h,       19},
    {NULL, NULL, 0}
};


void R_init_lqmm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
