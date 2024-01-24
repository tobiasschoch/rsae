// RegisteringDynamic Symbols
#include <R_ext/RS.h>
#include <R_ext/Rdynload.h>

// Declaration of the Fortran subroutines as void C functions (Fortran
// uses pass by reference, thus the arguments of the C functions are pointers;
// the macro F77_NAME takes care of naming conventions)
extern void F77_NAME(drlm)(void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *);
extern void F77_NAME(drsaehub)(void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *);
extern void F77_NAME(drsaehubpredict)(void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(drsaehubvariance)(void *, void *, void *, void *, void *,
    void *, void *, void *, void *);
extern void F77_NAME(drsaeresid)(void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *);

// Create array describing each routine
static const R_CMethodDef CMethods[] = {
    {"drlm",             (DL_FUNC) &F77_NAME(drlm),             10},
    {"drsaehub",         (DL_FUNC) &F77_NAME(drsaehub),         20},
    {"drsaehubpredict",  (DL_FUNC) &F77_NAME(drsaehubpredict),  14},
    {"drsaehubvariance", (DL_FUNC) &F77_NAME(drsaehubvariance),  9},
    {"drsaeresid",       (DL_FUNC) &F77_NAME(drsaeresid),       13},
    {NULL, NULL, 0}
};

// Register the routines to R as R_CMethodDef; they are called using .C()
void R_init_rsae(DllInfo *dll)
{
    R_registerRoutines(dll, CMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
