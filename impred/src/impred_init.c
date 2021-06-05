  
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP impred_randsetsMCMC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_randsetspred(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_sigmaSolvej(SEXP, SEXP, SEXP, SEXP);


static const R_CallMethodDef CallEntries[] = {
    {"BayesBD_BayesBDbinary", (DL_FUNC) &BayesBD_BayesBDbinary, 8},
    {"BayesBD_BayesBDnormal", (DL_FUNC) &BayesBD_BayesBDnormal, 9},
    {"BayesBD_eigenfun",      (DL_FUNC) &BayesBD_eigenfun,      2},
    {"BayesBD_unisliceL",     (DL_FUNC) &BayesBD_unisliceL,     8},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"impred_randsetsMCMC", (DL_FUNC) &impred_randsetsMCMC, 5},
    {"impred_randsetspred", (DL_FUNC) &impred_randsetspred, 8},
    {"impred_sigmaSolvej", (DL_FUNC) &impred_sigmaSolvej, 4}, 
    {NULL, NULL, 0}
};

void R_init_impred(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

