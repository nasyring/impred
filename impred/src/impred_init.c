  
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP impred_randsetsMCMC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_randsetspred(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_randsetspredlmer(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_genIM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_randsetspreddens(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_sigmaSolve(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_zeroin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_root_function(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_balanced(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_unbalanced(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_balanced_marginal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_unbalanced_marginal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_marginal_lmer(SEXP, SEXP, SEXP);
extern SEXP impred_plaus_balanced_aov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_unbalanced_aov(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_unbalanced_aov_full(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_two_stage(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_plaus_two_stage_full(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);



static const R_CallMethodDef CallEntries[] = {
    {"impred_randsetsMCMC", (DL_FUNC) &impred_randsetsMCMC, 5},
    {"impred_randsetspred", (DL_FUNC) &impred_randsetspred, 8},
    {"impred_randsetspredlmer", (DL_FUNC) &impred_randsetspredlmer, 8},
    {"impred_genIM", (DL_FUNC) &impred_genIM, 6},
    {"impred_randsetspreddens", (DL_FUNC) &impred_randsetspreddens, 11},
    {"impred_sigmaSolve", (DL_FUNC) &impred_sigmaSolve, 5}, 
    {"impred_zeroin", (DL_FUNC) &impred_zeroin, 8},
    {"impred_root_function", (DL_FUNC) &impred_root_function, 8}, 
    {"impred_plaus_balanced", (DL_FUNC) &impred_plaus_balanced, 10}, 
    {"impred_plaus_unbalanced", (DL_FUNC) &impred_plaus_unbalanced, 12},
    {"impred_plaus_balanced_marginal", (DL_FUNC) &impred_plaus_balanced_marginal, 10}, 
    {"impred_plaus_unbalanced_marginal", (DL_FUNC) &impred_plaus_unbalanced_marginal, 10}, 
    {"impred_plaus_marginal_lmer", (DL_FUNC) &impred_plaus_marginal_lmer, 3}, 
  {"impred_plaus_balanced_aov", (DL_FUNC) &impred_plaus_balanced_aov, 8},
  {"impred_plaus_unbalanced_aov", (DL_FUNC) &impred_plaus_unbalanced_aov, 10},
    {"impred_plaus_unbalanced_aov_full", (DL_FUNC) &impred_plaus_unbalanced_aov_full, 8},
  {"impred_plaus_two_stage", (DL_FUNC) &impred_plaus_two_stage, 9},
    {"impred_plaus_two_stage_full", (DL_FUNC) &impred_plaus_two_stage_full, 7},
    {NULL, NULL, 0}
};


void R_init_impred(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
