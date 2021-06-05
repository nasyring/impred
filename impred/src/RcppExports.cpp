#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;



// randsetsMCMC
Rcpp::List randsetsMCMC(NumericMatrix H, NumericMatrix A, NumericVector rL, NumericVector dimH, NumericVector M_samp);
RcppExport SEXP impred_randsetsMCMC(SEXP HSEXP, SEXP ASEXP, SEXP rLSEXP, SEXP dimHSEXP, SEXP M_sampSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type H(HSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rL(rLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dimH(dimHSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_samp(M_sampSEXP);
    __result = Rcpp::wrap(randsetsMCMC(H,A,rL,dimH,M_samp));
    return __result;
END_RCPP
}

// randsetspred
Rcpp::List randsetspred(NumericMatrix S, NumericVector dimS, NumericVector nsize, NumericVector n_i, NumericVector dimn_i, NumericVector k, NumericVector U, NumericVector Ybar);
RcppExport SEXP impred_randsetspred(SEXP SSEXP, SEXP dimSSEXP, SEXP nsizeSEXP, SEXP n_iSEXP, SEXP dimn_iSEXP, SEXP kSEXP, SEXP USEXP, SEXP YbarSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dimS(dimSSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nsize(nsizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dimn_i(dimn_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type U(USEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    __result = Rcpp::wrap(randsetspred(S,dimS,nsize,n_i,dimn_i,k,U,Ybar));
    return __result;
END_RCPP
}

Rcpp::NumericVector zeroin(NumericVector ax, NumericVector bx, NumericVector u, NumericVector v, NumericVector y, NumericVector z, NumericVector(*f)(NumericVector x, NumericVector uu, NumericVector vv, NumericVector yy, NumericVector zz), NumericVector tol);
RcppExport SEXP impred_zeroin(SEXP axSEXP, SEXP bxSEXP, SEXP uSEXP, SEXP vSEXP, SEXP ySEXP, SEXP zSEXP, SEXP fSEXP, SEXP tolSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type ax(axSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bx(bxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector (*)(Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::NumericVector) >::type f(fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tol(tolSEXP);
    __result = Rcpp::wrap(zeroin(ax, bx, u, v, y, z, f, tol));
    return __result;
END_RCPP
}

// sigmaSolvej
Rcpp::List sigmaSolvej(NumericVector Sampsj, NumericVector SL, NumericVector aL, NumericVector lambdaL);
RcppExport SEXP impred_sigmaSolvej(SEXP SampsjSEXP, SEXP SLSEXP, SEXP aLSEXP, SEXP lambdaLSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type Sampsj(SampsjSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SL(SLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type aL(aLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambdaL(lambdaLSEXP);
    __result = Rcpp::wrap(sigmaSolvej(Sampsj, SL, aL, lambdaL));
    return __result;
END_RCPP
}


static const R_CallMethodDef CallEntries[] = {
    {"impred_randsetsMCMC", (DL_FUNC) &impred_randsetsMCMC, 5},
    {"impred_randsetspred", (DL_FUNC) &impred_randsetspred, 8},
    {"impred_sigmaSolvej", (DL_FUNC) &impred_sigmaSolvej, 4}, 
    {NULL, NULL, 0}
};

RcppExport void R_init_impred(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
    
