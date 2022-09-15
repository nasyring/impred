  
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP impred_plaus_balanced_aov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_unbalanced_aov_full(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_two_stage_full(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_IMTS_mh_sampler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);



static const R_CallMethodDef CallEntries[] = {
  {"impred_plaus_balanced_aov", (DL_FUNC) &impred_plaus_balanced_aov, 8},
  {"impred_plaus_unbalanced_aov_full", (DL_FUNC) &impred_plaus_unbalanced_aov_full, 8},
  {"impred_plaus_two_stage_full", (DL_FUNC) &impred_plaus_two_stage_full, 7},
   {"impred_IMTS_mh_sampler", (DL_FUNC) &impred_IMTS_mh_sampler, 7},
  {NULL, NULL, 0}
};


void R_init_impred(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
