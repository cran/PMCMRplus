#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Fortran calls */
extern void F77_NAME(dstat)(double *x, double *d, int *n);
extern void F77_NAME(pd)(double *q, int *n, int*m, double *p);
extern void F77_NAME(pava)(double *y, double *w, int *kt, int *n);


static const R_FortranMethodDef FortranEntries[] = {
    {"dstat", (DL_FUNC) &F77_NAME(dstat), 3},
    {"pava",  (DL_FUNC) &F77_NAME(pava),  4},
    {"pd",    (DL_FUNC) &F77_NAME(pd),    4},
    {NULL, NULL, 0}
};

void R_init_PMCMRplus(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
