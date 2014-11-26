#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "genoset.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(binary_bound, 3),
    CALLDEF(binary_bound_by_chr, 8),
    CALLDEF(rangeMeans_numeric, 3),
    CALLDEF(rangeMeans_rle, 5),
    CALLDEF(numCallable_rle, 5),
    {NULL, NULL, 0}
};


void R_init_genoset(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
