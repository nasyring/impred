#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;



// randsetsMCMC
Rcpp::List randsetsMCMC(SEXP & H, SEXP & A, SEXP & rL, SEXP & dimH, SEXP & M_samp);
RcppExport SEXP impred_randsetsMCMC(SEXP HSEXP, SEXP SEXPA, SEXP SEXPrL, SEXP SEXPdimH, SEXP SEXPM_samp){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP & >::type H(HSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type A(ASEXP);
    Rcpp::traits::input_parameter< SEXP & >::type rL(rLSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type dimH(dimHSEXP);
    Rcpp::traits::input_parameter< SEXP & >::type M_samp(M_sampSEXP);
    __result = Rcpp::wrap(randsetsMCMC(H,A,rL,dimH,M_samp));
    return __result;
END_RCPP
}



static const R_CallMethodDef CallEntries[] = {
    {"impred_randsetsMCMC", (DL_FUNC) &impred_randsetsMCMC, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_impred(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
    
