#include "RcppArmadillo.h"
#include <RcppParallel.h>
#include <Rcpp.h>
#include <math.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>




Rcpp::List randsetsMCMC(NumericMatrix H, NumericMatrix A, NumericVector rL, NumericVector dimH, NumericVector M_samp) {
	
	List result;
	int M = int(M_samp[0]);
	int H2 = int(dimH[0]+2);
	int H1 = H2-2;
	NumericVector propsd(1,0.0);
	NumericVector u(2,1.0);
	NumericVector uprop(2,0.0);
	NumericVector logjointold(1,0.0);
	NumericVector logjointnew(1,0.0);
	NumericVector logjointdiff(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector zeroes = NumericVector(H2*1, 0.0); 
        NumericMatrix U = NumericMatrix(H2, 1, zeroes.begin());
	arma::mat Aa = as<arma::mat>(A);
	arma::mat arg;
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	arma::mat Ua = as<arma::mat>(U);
	
	for(int h = 0; h<H1; h++){
		U(h,0) = H(h,0);	
	}

for(int j=0; j<M; j++) {
		for(int i=0; i<2; i++){
			if( i==0 ){
				propsd[0] = 10.0;	
			}else {
				propsd[0] = 0.5;	
			}
			uprop(0) = u(0);uprop(1) = u(1);
			U(H1,0) = uprop[0];
			U(H1+1,0) = uprop[1];
			Ua = as<arma::mat>(U);
			arg = Aa*Ua;
			for(int h = 0; h<H2; h++){
				logjointold[0] = logjointold[0] + 0.5*rL[h]*arg(h,0)-0.5*exp(arg(h,0));	
			}
			uprop[i] = R::rnorm(u[i], propsd[0]);
			U(H1,0) = uprop[0];
			U(H1+1,0) = uprop[1];
			Ua = as<arma::mat>(U);
			arg = Aa*Ua;
			for(int h = 0; h<H2; h++){
				logjointnew[0] = logjointnew[0] + 0.5*rL[h]*arg(h,0)-0.5*exp(arg(h,0));	
			}
			logjointdiff[0] = logjointnew[0] - logjointold[0];
			logjointdiff[0] = fmin(std::exp(logjointdiff[0]), 1.0);
			uu[0] = R::runif(0.0,1.0);
			if(uu[0] <= logjointdiff[0]) {
				if(i==0){
					postsamples0[j] = uprop[0];	
				}else {
					postsamples1[j] = uprop[1];
				}
				u(0)=uprop(0);u(1)=uprop(1);
			}else {
				if(i==0){
					postsamples0[j] = u[0];	
				}else {
					postsamples1[j] = u[1];
				}				
			}
			logjointold[0] = 0.0; logjointnew[0] = 0.0;
		}
	}
result = Rcpp::List::create(Rcpp::Named("samples1") = postsamples0,Rcpp::Named("samples2") = postsamples1);

	return result;
	
	
}



Rcpp::List randsetspred(NumericMatrix S, NumericVector dimS, NumericVector nsize, NumericVector n_i, NumericVector dimn_i, NumericVector k, NumericVector U, NumericVector Ybar) {
	
	List result;
	int M = int(dimS[0]);
	int n = int(nsize[0]);
	int dn_i = int(dimn_i[0]);
	NumericVector sumn_i2(1,0.0);
	NumericVector randnorms = rnorm(M);
	NumericVector Qsamps(M,0.0);
	NumericVector Ql(1,0.0);
	NumericVector Qu(1,0.0);
	NumericVector Ul(1,0.0);
	NumericVector Uu(1,0.0);
	NumericVector zeroes = NumericVector(M*2, 0.0); 
        NumericMatrix randsetpred = NumericMatrix(M, 2, zeroes.begin());
	
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
	
	/*
	for(int j=0; j < M; j++){
		Qsamps[j] = randnorms[j]*std::sqrt(S(j,0)*(1+(1/(n*n))*sumn_i2[0])+S(j,1)*((1/n)+1/(k[0])));	
	}

	std::sort(Qsamps.begin(), Qsamps.end());
	
	for(int j=0; j < M; j++){
		Ul[0] = 0.5-std::fabs(U[j]-0.5);
		Uu[0] = 1.0-Ul[0];
		Ql[0] = Qsamps[std::round(M*Ul[0])];
		Qu[0] = Qsamps[std::round(M*Uu[0])];
		randsetpred(j,0) = Ybar[0]+Ql[0];
		randsetpred(j,1) = Ybar[0]+Qu[0];
	}
	*/

result = Rcpp::List::create(Rcpp::Named("randsetpred") = randsetpred);

	return result;
	
	
}


