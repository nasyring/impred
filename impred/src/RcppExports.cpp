#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;



    
Rcpp::List IMTS_mh_sampler(NumericVector lU0, NumericVector V0, NumericVector H0, NumericMatrix Minv, NumericVector rL, NumericVector thetau0s, NumericVector prosd1, NumericVector prosd2);
RcppExport SEXP impred_IMTS_mh_sampler(SEXP lU0SEXP, SEXP V0SEXP, SEXP H0SEXP, SEXP MinvSEXP, SEXP rLSEXP, SEXP thetau0sSEXP, SEXP propsd1SEXP, SEXP propsd2SEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type lU0(lU0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type V0(V0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type H0(H0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Minv(MinvSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rL(rLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thetau0s(thetau0sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsd1(propsd1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propsd2(propsd2SEXP);
    __result = Rcpp::wrap(IMTS_mh_sampler(lU0, V0, H0, Minv, rL, thetau0s, propsd1, propsd2));
    return __result;
END_RCPP
}    
    

Rcpp::List plaus_unbalanced_aov_full(NumericVector theta, NumericVector Ybar, NumericVector S, NumericVector lambda, NumericVector r, NumericVector n, NumericVector n_i, NumericVector ratio);
RcppExport SEXP impred_plaus_unbalanced_aov_full(SEXP thetaSEXP, SEXP YbarSEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP rSEXP, SEXP nSEXP, SEXP n_iSEXP, SEXP ratioSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ratio(ratioSEXP);
    __result = Rcpp::wrap(plaus_unbalanced_aov_full(theta, Ybar, S, lambda, r, n, n_i, ratio));
    return __result;
END_RCPP
}

Rcpp::List plaus_two_stage_full(NumericVector theta, NumericVector xBy, NumericVector S, NumericVector lambda, NumericVector r, NumericVector csigma, NumericVector ratio);
RcppExport SEXP impred_plaus_two_stage_full(SEXP thetaSEXP, SEXP xBySEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP rSEXP, SEXP csigmaSEXP, SEXP ratioSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xBy(xBySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type csigma(csigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ratio(ratioSEXP);
    __result = Rcpp::wrap(plaus_two_stage_full(theta, xBy, S, lambda, r, csigma, ratio));
    return __result;
END_RCPP
}




Rcpp::List plaus_balanced_aov(NumericVector theta, NumericVector Ybar, NumericVector S, NumericVector lambda, NumericVector r, NumericVector n, NumericVector n_i, NumericVector eta);
RcppExport SEXP impred_plaus_balanced_aov(SEXP thetaSEXP, SEXP YbarSEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP rSEXP, SEXP nSEXP, SEXP n_iSEXP, SEXP etaSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);    
    __result = Rcpp::wrap(plaus_balanced_aov(theta, Ybar, S, lambda, r, n, n_i, eta));
    return __result;
END_RCPP
}
