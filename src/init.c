#include <R_ext/RS.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern void F77_NAME(drlm)(void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *);
extern void F77_NAME(drsaehub)(void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *);

extern void F77_NAME(drsaehubpredict)(void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(drsaehubvariance)(void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *);
extern void F77_NAME(drsaeresid)(void *, void *, void *, void *, void *,
    void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"drlm",             (DL_FUNC) &F77_NAME(drlm),             10},
    {"drsaehub",         (DL_FUNC) &F77_NAME(drsaehub),         20},
    {"drsaehubpredict",  (DL_FUNC) &F77_NAME(drsaehubpredict),  14},
    {"drsaehubvariance", (DL_FUNC) &F77_NAME(drsaehubvariance), 10},
    {"drsaeresid",       (DL_FUNC) &F77_NAME(drsaeresid),       13},
    {NULL, NULL, 0}
};

void R_init_rsae(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
