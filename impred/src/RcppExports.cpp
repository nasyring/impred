#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;



// randsetsMCMC
Rcpp::List randsetsMCMC(RMatrix<double> H, RMatrix<double> A, RVector<double> rL, RVector<double> dimH, RVector<double> M_samp);
RcppExport SEXP impred_randsetsMCMC(SEXP HSEXP, SEXP ASEXP, SEXP rLSEXP, SEXP dimHSEXP, SEXP M_sampSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< RMatrix<double> >::type H(HSEXP);
    Rcpp::traits::input_parameter< RMatrix<double> >::type A(ASEXP);
    Rcpp::traits::input_parameter< RVector<double> >::type rL(rLSEXP);
    Rcpp::traits::input_parameter< RVector<double> >::type dimH(dimHSEXP);
    Rcpp::traits::input_parameter< RVector<double> >::type M_samp(M_sampSEXP);
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
    
