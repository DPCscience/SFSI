#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP crossprod2distance(SEXP, SEXP);
extern SEXP delete_col(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gaussian_kernel(SEXP, SEXP, SEXP);
extern SEXP getCorrelated(SEXP, SEXP, SEXP);
extern SEXP laplacian_kernel(SEXP, SEXP, SEXP);
extern SEXP polynomial_kernel(SEXP, SEXP, SEXP, SEXP);
extern SEXP readBinFile(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scaleXtX(SEXP, SEXP);
extern SEXP updatebeta_lambda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP writeBinFile(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"crossprod2distance", (DL_FUNC) &crossprod2distance, 2},
    {"delete_col",         (DL_FUNC) &delete_col,         5},
    {"gaussian_kernel",    (DL_FUNC) &gaussian_kernel,    3},
    {"getCorrelated",      (DL_FUNC) &getCorrelated,      3},
    {"laplacian_kernel",   (DL_FUNC) &laplacian_kernel,   3},
    {"polynomial_kernel",  (DL_FUNC) &polynomial_kernel,  4},
    {"readBinFile",        (DL_FUNC) &readBinFile,        5},
    {"scaleXtX",           (DL_FUNC) &scaleXtX,           2},
    {"updatebeta_lambda",  (DL_FUNC) &updatebeta_lambda,  9},
    {"writeBinFile",       (DL_FUNC) &writeBinFile,       5},
    {NULL, NULL, 0}
};

void R_init_SFSI(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}