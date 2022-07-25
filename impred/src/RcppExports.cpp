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

Rcpp::List plaus_unbalanced_aov(NumericVector theta, NumericVector Ybar, NumericVector S, NumericVector lambda, NumericVector n, NumericVector n_i, NumericMatrix auxiliary, NumericVector s2a, NumericVector s2e);
RcppExport SEXP impred_plaus_unbalanced_aov(SEXP thetaSEXP, SEXP YbarSEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP nSEXP, SEXP n_iSEXP, SEXP auxiliarySEXP, SEXP s2aSEXP, SEXP s2eSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type auxiliary(auxiliarySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2a(s2aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2e(s2eSEXP);
    __result = Rcpp::wrap(plaus_unbalanced_aov(theta, Ybar, S, lambda, n, n_i, auxiliary, s2a, s2e));
    return __result;
END_RCPP
}



Rcpp::List plaus_two_stage(NumericVector theta, NumericVector xBy, NumericVector S, NumericVector lambda, NumericMatrix auxiliary, NumericVector csigma, NumericVector s2a, NumericVector s2e);
RcppExport SEXP impred_plaus_two_stage(SEXP thetaSEXP, SEXP xBySEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP auxiliarySEXP, SEXP csigmaSEXP, SEXP s2aSEXP, SEXP s2eSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xBy(xBySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type auxiliary(auxiliarySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type csigma(csigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2a(s2aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type s2e(s2eSEXP);
    __result = Rcpp::wrap(plaus_two_stage(theta, xBy, S, lambda, auxiliary, csigma, s2a, s2e));
    return __result;
END_RCPP
}


Rcpp::List plaus_balanced(NumericVector thetaseq, NumericVector saseq, NumericVector seseq, NumericVector n, NumericVector n_i, NumericVector S, NumericVector lambda, NumericVector r, NumericVector Ybar, NumericVector numk);
RcppExport SEXP impred_plaus_balanced(SEXP thetaseqSEXP, SEXP saseqSEXP, SEXP seseqSEXP, SEXP nSEXP, SEXP n_iSEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP rSEXP, SEXP YbarSEXP, SEXP numkSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaseq(thetaseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type saseq(saseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type seseq(seseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type numk(numkSEXP);
    __result = Rcpp::wrap(plaus_balanced(thetaseq,saseq,seseq,n,n_i,S,lambda,r,Ybar,numk));
    return __result;
END_RCPP
}

Rcpp::List plaus_balanced_marginal(NumericVector thetaseq, NumericVector n, NumericVector n_i, NumericVector S, NumericVector lambda, NumericVector r, NumericVector Ybar, NumericVector numk, NumericVector sa2, NumericVector se2);
RcppExport SEXP impred_plaus_balanced_marginal(SEXP thetaseqSEXP, SEXP nSEXP, SEXP n_iSEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP rSEXP, SEXP YbarSEXP, SEXP numkSEXP, SEXP sa2SEXP, SEXP se2SEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaseq(thetaseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type numk(numkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sa2(sa2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se2(se2SEXP);
    __result = Rcpp::wrap(plaus_balanced_marginal(thetaseq,n,n_i,S,lambda,r,Ybar,numk,sa2,se2));
    return __result;
END_RCPP
}

Rcpp::List plaus_unbalanced_marginal(NumericVector thetaseq, NumericVector n, NumericVector n_i, NumericVector S, NumericVector lambda, NumericVector r, NumericVector Ybar, NumericVector numk, NumericVector sa2, NumericVector se2);
RcppExport SEXP impred_plaus_unbalanced_marginal(SEXP thetaseqSEXP, SEXP nSEXP, SEXP n_iSEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP rSEXP, SEXP YbarSEXP, SEXP numkSEXP, SEXP sa2SEXP, SEXP se2SEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaseq(thetaseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type numk(numkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sa2(sa2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type se2(se2SEXP);
    __result = Rcpp::wrap(plaus_unbalanced_marginal(thetaseq,n,n_i,S,lambda,r,Ybar,numk,sa2,se2));
    return __result;
END_RCPP
}

Rcpp::List plaus_marginal_lmer(NumericVector thetaseq, NumericVector total_sigma, NumericVector xBy);
RcppExport SEXP impred_plaus_marginal_lmer(SEXP thetaseqSEXP, SEXP total_sigmaSEXP, SEXP xBySEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type thetaseq(thetaseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type total_sigma(total_sigmaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xBy(xBySEXP);
    __result = Rcpp::wrap(plaus_marginal_lmer(thetaseq,total_sigma,xBy));
    return __result;
END_RCPP
}



Rcpp::List plaus_unbalanced(NumericVector thetaseq, NumericVector saseq, NumericVector seseq, NumericVector n, NumericVector n_i, NumericVector S, NumericVector lambda, NumericVector r, NumericVector Ybar, NumericVector numk, NumericVector samples1, NumericVector samples2);
RcppExport SEXP impred_plaus_unbalanced(SEXP thetaseqSEXP, SEXP saseqSEXP, SEXP seseqSEXP, SEXP nSEXP, SEXP n_iSEXP, SEXP SSEXP, SEXP lambdaSEXP, SEXP rSEXP, SEXP YbarSEXP, SEXP numkSEXP, SEXP samples1SEXP, SEXP samples2SEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type thetaseq(thetaseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type saseq(saseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type seseq(seseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type numk(numkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samples1(samples1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type samples2(samples2SEXP);
    __result = Rcpp::wrap(plaus_unbalanced(thetaseq,saseq,seseq,n,n_i,S,lambda,r,Ybar, numk, samples1, samples2));
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




Rcpp::List randsetspredlmer(NumericMatrix S, NumericVector dimS, NumericVector U, NumericMatrix C1, NumericMatrix C2, NumericVector By, NumericVector x, NumericVector ztz);
RcppExport SEXP impred_randsetspredlmer(SEXP SSEXP, SEXP dimSSEXP, SEXP USEXP, SEXP C1SEXP, SEXP C2SEXP, SEXP BySEXP, SEXP xSEXP, SEXP ztzSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type S(SSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dimS(dimSSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type U(USEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type C1(C1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type C2(C2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type By(BySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ztz(ztzSEXP);
     __result = Rcpp::wrap(randsetspredlmer(S,dimS,U,C1,C2,By,x,ztz));
    return __result;
END_RCPP
}


Rcpp::List genIM(NumericVector Y, NumericMatrix Z, NumericVector museq, NumericVector saseq, NumericVector seseq, NumericVector M);
RcppExport SEXP impred_genIM(SEXP YSEXP, SEXP ZSEXP, SEXP museqSEXP, SEXP saseqSEXP, SEXP seseqSEXP, SEXP MSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type museq(museqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type saseq(saseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type seseq(seseqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M(MSEXP);
    __result = Rcpp::wrap(genIM( Y,  Z,  museq,  saseq, seseq,  M));
    return __result;
END_RCPP
}

Rcpp::List randsetspreddens(NumericVector sigsampdens, NumericVector dimS, NumericVector nsize, NumericVector n_i, NumericVector dimn_i, NumericVector k, NumericVector Ybar, NumericVector predgrid, NumericVector dim_predgrid, NumericVector localpt, NumericVector logdenslocalpt);
RcppExport SEXP impred_randsetspreddens(SEXP sigsampdensSEXP, SEXP dimSSEXP, SEXP nsizeSEXP, SEXP n_iSEXP, SEXP dimn_iSEXP, SEXP kSEXP, SEXP YbarSEXP, SEXP predgridSEXP, SEXP dim_predgridSEXP, SEXP localptSEXP, SEXP logdenslocalptSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type sigsampdens(sigsampdensSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dimS(dimSSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nsize(nsizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_i(n_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dimn_i(dimn_iSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ybar(YbarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type predgrid(predgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dim_predgrid(dim_predgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type localpt(localptSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type logdenslocalpt(logdenslocalptSEXP);
    __result = Rcpp::wrap(randsetspreddens(sigsampdens,dimS,nsize,n_i,dimn_i,k,Ybar,predgrid,dim_predgrid,localpt,logdenslocalpt));
    return __result;
END_RCPP
}

Rcpp::NumericVector zeroin(NumericVector ax, NumericVector bx, NumericVector u, NumericVector v, NumericVector y, NumericVector z, Function f, NumericVector tol);
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
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type tol(tolSEXP);
    __result = Rcpp::wrap(zeroin(ax, bx, u, v, y, z, f, tol));
    return __result;
END_RCPP
}



    
Rcpp::NumericVector root_function(NumericVector x, NumericVector Sampsj, NumericVector SL, NumericVector aL, NumericVector lambdaL);
RcppExport SEXP impred_root_function(SEXP xSEXP, SEXP SampsjSEXP, SEXP SLSEXP, SEXP aLSEXP, SEXP lambdaLSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Sampsj(SampsjSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SL(SLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type aL(aLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambdaL(lambdaLSEXP);
    __result = Rcpp::wrap(root_function(x, Sampsj, SL, aL, lambdaL));
    return __result;
END_RCPP
}

// sigmaSolvej
Rcpp::List sigmaSolve(NumericMatrix Samps, NumericVector SL, NumericVector aL, NumericVector aM, NumericVector lambdaL);
RcppExport SEXP impred_sigmaSolve(SEXP SampsSEXP, SEXP SLSEXP, SEXP aLSEXP,SEXP aMSEXP, SEXP lambdaLSEXP){
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type Samps(SampsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type SL(SLSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type aL(aLSEXP);
        Rcpp::traits::input_parameter< NumericVector >::type aM(aMSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambdaL(lambdaLSEXP);
    __result = Rcpp::wrap(sigmaSolve(Samps, SL, aL, aM, lambdaL));
    return __result;
END_RCPP
}

/*
static const R_CallMethodDef CallEntries[] = {
    {"impred_randsetsMCMC", (DL_FUNC) &impred_randsetsMCMC, 5},
    {"impred_randsetspred", (DL_FUNC) &impred_randsetspred, 8},
    {"impred_zeroin", (DL_FUNC) &impred_zeroin, 8},
    {"impred_root_function", (DL_FUNC) &impred_root_function, 5},
    {"impred_sigmaSolvej", (DL_FUNC) &impred_sigmaSolvej, 4}, 
    {NULL, NULL, 0}
};

RcppExport void R_init_impred(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

*/
