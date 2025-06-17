#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP C_ccAFv2(SEXP inp_vec);

static const R_CallMethodDef CallEntries[] = {
    {"C_ccAFv2", (DL_FUNC) &C_ccAFv2, 1},
    {NULL, NULL, 0}
};

void R_init_ccAFv2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
