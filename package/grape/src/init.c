
// See Sections 5.4 and 6.16 of Writing R Extensions

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "grape.h"

static R_NativePrimitiveArgType grape_types[7] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_CMethodDef cMethods[] = {
    {"grape", (DL_FUNC) &grape, 7, grape_types},
    {NULL, NULL, 0, NULL}
};

void attribute_visible R_init_grape(DllInfo *info)
{
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

