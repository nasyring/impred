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
	NumericVector postsamplesdens(M,0.0);
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
				postsamplesdens[j] = logjointnew[0];
			}else {
				if(i==0){
					postsamples0[j] = u[0];	
				}else {
					postsamples1[j] = u[1];
				}
				postsamplesdens[j] = logjointold[0];
			}
			logjointold[0] = 0.0; logjointnew[0] = 0.0;
		}
	}
	result = Rcpp::List::create(Rcpp::Named("samples1") = postsamples0,Rcpp::Named("samples2") = postsamples1,Rcpp::Named("logdens") = postsamplesdens);

	return result;
	
	
}

Rcpp::List plaus_balanced_aov(NumericVector theta, NumericVector Ybar, NumericVector S, NumericVector lambda, NumericVector r, NumericVector n, NumericVector n_i, NumericVector eta){

	List result;
	int m_the = theta.length();
	int dn_i = n_i.length();
	int m_eta = eta.length();

	NumericVector sumn_i2(1, 0.0);
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
	
	NumericVector c1t(1,0.0); NumericVector c2t(1,0.0); NumericVector ct(1,0.0);
	NumericVector c1n(1,0.0); NumericVector c2n(1,0.0); NumericVector cn(1,0.0);
	NumericVector c1e(1,0.0); NumericVector c2e(1,0.0); NumericVector ce(1,0.0);
	c2t[0] = 1/n[0];
	c1t[0] = 1+sumn_i2[0]*(c2t[0]*c2t[0]);
	ct[0] = (c2t[0]*lambda[0]/c1t[0]) - 1.0;
	c2n[0] = (1.0/n[0] + 1.0);
	c1n[0] = 1+sumn_i2[0]/(n[0]*n[0]);
	cn[0] = (c2n[0]*lambda[0]/c1n[0]) - 1.0;
	c2e[0] = (1.0/n[0] + 1.0);
	c1e[0] = 1+(sumn_i2[0]/(n[0]*n[0])) - (2.0*n_i[dn_i-1]/n[0]);
	ce[0] = (c2e[0]*lambda[0]/c1e[0]) - 1.0;
	
	NumericVector plausseq_t(m_the, 0.0);
	NumericVector plausseq_n(m_the, 0.0);
	NumericVector plausseq_e(m_the, 0.0);
	NumericVector den_t(1, 0.0);NumericVector den_n(1, 0.0);NumericVector den_e(1, 0.0);
	NumericVector num(1, 0.0);
	den_t[0] = ((c1t[0]/lambda[0])*((S[0]/r[0]) - (S[1]/r[1])) + (c2t[0]*S[1]/r[1]));
	den_n[0] = ((c1n[0]/lambda[0])*((S[0]/r[0]) - (S[1]/r[1])) + (c2n[0]*S[1]/r[1]));
	den_e[0] = ((c1e[0]/lambda[0])*((S[0]/r[0]) - (S[1]/r[1])) + (c2e[0]*S[1]/r[1]));
	for(int j=0; j<m_the; j++){
		num[0] = (theta[j] - Ybar[0])*(theta[j] - Ybar[0]);
		plausseq_t[j] = num[0]/den_t[0];
		plausseq_n[j] = num[0]/den_n[0];
		plausseq_e[j] = num[0]/den_e[0];
	}
	
	NumericVector Z2(10000,0.0); NumericVector V1(10000,0.0); NumericVector V2(10000,0.0); 
	Z2 = Rcpp::rchisq(10000,1.0);
	V1 = Rcpp::rchisq(10000,r[0]);
	V2 = Rcpp::rchisq(10000,r[1]);
	
	NumericVector zeroes(10000*101,0.0);
	NumericMatrix MCt = NumericMatrix(10000, m_eta, zeroes.begin());
	NumericMatrix MCn = NumericMatrix(10000, m_eta, zeroes.begin());
	NumericMatrix MCe = NumericMatrix(10000, m_eta, zeroes.begin());

	for(int i = 0; i < 10000; i++){
		for(int j = 0; j < m_eta; j++){
			MCt(i,j) = Z2[i]/((V1[i]/r[0])*(1.0/(1.0+ct[0]*eta[j])) + (V2[i]/r[1])*(ct[0]*eta[j]/(1.0+ct[0]*eta[j])));
			MCn(i,j) = Z2[i]/((V1[i]/r[0])*(1.0/(1.0+cn[0]*eta[j])) + (V2[i]/r[1])*(cn[0]*eta[j]/(1.0+cn[0]*eta[j])));
			MCe(i,j) = Z2[i]/((V1[i]/r[0])*(1.0/(1.0+ce[0]*eta[j])) + (V2[i]/r[1])*(ce[0]*eta[j]/(1.0+ce[0]*eta[j])));
		}
	}
	
	
	NumericVector zeroes2(m_the*101,0.0);
	NumericMatrix F_eta_t = NumericMatrix(m_the, m_eta, zeroes.begin());
	NumericMatrix F_eta_n = NumericMatrix(m_the, m_eta, zeroes.begin());
	NumericMatrix F_eta_e = NumericMatrix(m_the, m_eta, zeroes.begin());

	for(int i = 0; i < m_the; i++){
		for(int j = 0; j < m_eta; j++){
			for(int k = 0; k < 10000; k++){
				if(MCt(k,j) < plausseq_t[i]){
					F_eta_t(i,j) = F_eta_t(i,j) + 0.0001;
				}
				if(MCn(k,j) < plausseq_n[i]){
					F_eta_n(i,j) = F_eta_n(i,j) + 0.0001;
				}
				if(MCe(k,j) < plausseq_e[i]){
					F_eta_e(i,j) = F_eta_e(i,j) + 0.0001;
				}
			}
		}
	}
	
	NumericVector F_max_t(m_the, 0.0); 
	NumericVector F_min_t(m_the, 1.0);
	NumericVector F_max_n(m_the, 0.0); 
	NumericVector F_min_n(m_the, 1.0);
	NumericVector F_max_e(m_the, 0.0); 
	NumericVector F_min_e(m_the, 1.0);
	for(int i = 0; i < m_the; i++){
		for(int j = 0; j < m_eta; j++){
			F_max_t[i] = std::max(F_max_t[i], F_eta_t(i, j));
			F_min_t[i] = std::min(F_min_t[i], F_eta_t(i, j));
			F_max_n[i] = std::max(F_max_n[i], F_eta_n(i, j));
			F_min_n[i] = std::min(F_min_n[i], F_eta_n(i, j));
			F_max_e[i] = std::max(F_max_e[i], F_eta_e(i, j));
			F_min_e[i] = std::min(F_min_e[i], F_eta_e(i, j));
		}
	}
	
	NumericVector plaus_t(m_the, 0.0);
	NumericVector maxplaus_t(1, 0.0);
	NumericVector plaus_n(m_the, 0.0);
	NumericVector maxplaus_n(1, 0.0);
	NumericVector plaus_e(m_the, 0.0);
	NumericVector maxplaus_e(1, 0.0);

	if(den_t[0] < 0){
		for(int i = 0; i < m_the; i++){
			plaus_t[i] = F_max_t[i];
			maxplaus_t[0] = std::max(maxplaus_t[0],plaus_t[i]);
		}
	}else {
		for(int i = 0; i < m_the; i++){
			plaus_t[i] = 1.0-F_min_t[i];
			maxplaus_t[0] = std::max(maxplaus_t[0],plaus_t[i]);
		}		
	}

	if(den_n[0] < 0){
		for(int i = 0; i < m_the; i++){
			plaus_n[i] = F_max_n[i];
			maxplaus_n[0] = std::max(maxplaus_n[0],plaus_n[i]);
		}
	}else {
		for(int i = 0; i < m_the; i++){
			plaus_n[i] = 1.0-F_min_n[i];
			maxplaus_n[0] = std::max(maxplaus_n[0],plaus_n[i]);
		}		
	}
	
	if(den_e[0] < 0){
		for(int i = 0; i < m_the; i++){
			plaus_e[i] = F_max_e[i];
			maxplaus_e[0] = std::max(maxplaus_e[0],plaus_e[i]);
		}
	}else {
		for(int i = 0; i < m_the; i++){
			plaus_e[i] = 1.0-F_min_e[i];
			maxplaus_e[0] = std::max(maxplaus_e[0],plaus_e[i]);
		}		
	}	
	
	for(int i = 0; i < m_the; i++){
		plaus_t[i] = plaus_t[i]/maxplaus_t[0];
		plaus_n[i] = plaus_n[i]/maxplaus_n[0];
		plaus_e[i] = plaus_e[i]/maxplaus_e[0];
	}
	
	result = Rcpp::List::create(Rcpp::Named("plauses.theta") = plaus_t, Rcpp::Named("plauses.new") = plaus_n, Rcpp::Named("plauses.exs") = plaus_e);
	return result;
	
}


Rcpp::List plaus_unbalanced_aov(NumericVector theta, NumericVector Ybar, NumericVector S, NumericVector lambda, NumericVector n, NumericVector n_i, NumericMatrix auxiliary, NumericVector s2a, NumericVector s2e, NumericVector assoc){

	List result;
	int L = S.length();
	int m_the = theta.length();
	int m_samps = auxiliary.nrow();
	int dn_i = n_i.length();

	NumericVector sumn_i2(1, 0.0);
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
	
	NumericVector c1t(1,0.0); NumericVector c2t(1,0.0); 
	NumericVector c1n(1,0.0); NumericVector c2n(1,0.0); 
	NumericVector c1e(1,0.0); NumericVector c2e(1,0.0); 
	c2t[0] = 1/n[0];
	c1t[0] = 1+sumn_i2[0]*(c2t[0]*c2t[0]);
	c2n[0] = (1.0/n[0] + 1.0);
	c1n[0] = 1+sumn_i2[0]/(n[0]*n[0]);
	c2e[0] = (1.0/n[0] + 1.0);
	c1e[0] = 1+(sumn_i2[0]/(n[0]*n[0])) - (2.0*n_i[dn_i-1]/n[0]);

	NumericVector MC(1, 0.0);
	NumericVector MCt(m_samps, 0.0);
	NumericVector MCn(m_samps, 0.0);
	NumericVector MCe(m_samps, 0.0);
	NumericVector Z2(1, 0.0);
	NumericVector plausseq(m_the, 0.0);
	
	if(assoc[0]==0.0){

		NumericVector prodS(1, 1.0);
		NumericVector prodvar(1, 1.0);
		for(int j = 0; j < (L-1); j++){
			prodS[0] = prodS[0]*S[j];
			prodvar[0] = prodvar[0]*(lambda[j]*s2a[0] + s2e[0]); 
		}	
		prodS[0] = prodS[0] + S[L-1];	
		
		for(int j = 0; j < m_the; j++){
			plausseq[j] = (theta[j] - Ybar[0])*(theta[j] - Ybar[0])/(prodS[0]);
		}	
		
		for(int j = 0; j < m_samps; j++){
			Z2[0] = R::rchisq(1.0);
			MC[0] = Z2[0]/(prodvar[0]*std::exp(auxiliary(j,0)) + s2e[0]*std::exp(auxiliary(j,1)));
			MCt[j] = MC[0]*(c1t[0]*s2a[0] + c2t[0]*s2e[0]);
			MCn[j] = MC[0]*(c1n[0]*s2a[0] + c2n[0]*s2e[0]);
			MCe[j] = MC[0]*(c1e[0]*s2a[0] + c2e[0]*s2e[0]);
		}
		
	}else {
		NumericVector slogS(1, 0.0);
		NumericVector slogvar(1, 0.0);
		for(int j = 0; j < (L-1); j++){
			slogS[0] = slogS[0] + std::log(S[j]);
			slogvar[0] = slogvar[0]+std::log(lambda[j]*s2a[0] + s2e[0]); 
		}	
		slogS[0] = slogS[0] + std::log(S[L-1]);		
		
		for(int j = 0; j < m_the; j++){
			plausseq[j] = (theta[j] - Ybar[0])/slogS[0];
		}	
		
		for(int j = 0; j < m_samps; j++){
			Z2[0] = R::rnorm(0.0,1.0);
			MC[0] = Z2[0]/(slogvar[0] + auxiliary(j,0) + std::log(s2e[0]) + auxiliary(j,1));	
			MCt[j] = MC[0]*std::sqrt(c1t[0]*s2a[0] + c2t[0]*s2e[0]);
			MCn[j] = MC[0]*std::sqrt(c1n[0]*s2a[0] + c2n[0]*s2e[0]);
			MCe[j] = MC[0]*std::sqrt(c1e[0]*s2a[0] + c2e[0]*s2e[0]);
		}
	}
	
	
	NumericVector Ft(m_the, 0.0); 
	NumericVector Fn(m_the, 0.0); 
	NumericVector Fe(m_the, 0.0);	
	for(int i = 0; i < m_the; i++){
		for(int j = 0; j < m_samps; j++){
			if(MCt[j] < plausseq[i]){
				Ft[i] = Ft[i] + (1.0 / m_samps);	
			}
			if(MCn[j] < plausseq[i]){
				Fn[i] = Fn[i] + (1.0 / m_samps);	
			}
			if(MCe[j] < plausseq[i]){
				Fe[i] = Fe[i] + (1.0 / m_samps);	
			}
		}
	}
	
	NumericVector plaus_t(m_the, 0.0); 
	NumericVector plaus_n(m_the, 0.0); 
	NumericVector plaus_e(m_the, 0.0);	

	if(assoc[0] == 0.0){
		for(int i = 0; i < m_the; i++){
			plaus_t[i] = 1.0 - Ft[i];
			plaus_n[i] = 1.0 - Fn[i];
			plaus_e[i] = 1.0 - Fe[i];
		}
	}else {
		for(int i = 0; i < m_the; i++){
			plaus_t[i] = 1.0 - std::abs(2.0*Ft[i]-1.0);
			plaus_n[i] = 1.0 - std::abs(2.0*Fn[i]-1.0);
			plaus_e[i] = 1.0 - std::abs(2.0*Fe[i]-1.0);
		}
	}
	
	result = Rcpp::List::create(Rcpp::Named("plauses.theta") = plaus_t, Rcpp::Named("plauses.new") = plaus_n, Rcpp::Named("plauses.exs") = plaus_e);
	return result;
	
}

Rcpp::List plaus_unbalanced_aov_full(NumericVector theta, NumericVector Ybar, NumericVector S, NumericVector lambda, NumericVector r, NumericVector n, NumericVector n_i, NumericVector ratio){

	List result;
	int L = S.length();
	int m_the = theta.length();
	int dn_i = n_i.length();
	int m_samps = 10000;
	
	NumericVector sumn_i2(1, 0.0);
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
	
	NumericVector c1t(1,0.0); NumericVector c2t(1,0.0); 
	NumericVector c1n(1,0.0); NumericVector c2n(1,0.0); 
	NumericVector c1e(1,0.0); NumericVector c2e(1,0.0); 
	c2t[0] = 1/n[0];
	c1t[0] = 1+sumn_i2[0]*(c2t[0]*c2t[0]);
	c2n[0] = (1.0/n[0] + 1.0);
	c1n[0] = 1+sumn_i2[0]/(n[0]*n[0]);
	c2e[0] = (1.0/n[0] + 1.0);
	c1e[0] = 1+(sumn_i2[0]/(n[0]*n[0])) - (2.0*n_i[dn_i-1]/n[0]);

	NumericVector MC(1, 0.0);
	NumericVector MCt(m_samps, 0.0);
	NumericVector MCn(m_samps, 0.0);
	NumericVector MCe(m_samps, 0.0);
	NumericVector Z2(1, 0.0);
	NumericVector plausseq(m_the, 0.0);
	
	NumericVector sumS(1, 0.0);
	for(int j = 0; j < L; j++){
		sumS[0] = sumS[0] + S[j];
	}	

	for(int j = 0; j < m_the; j++){
		plausseq[j] = (theta[j] - Ybar[0])*(theta[j] - Ybar[0])/(sumS[0]);
	}	
	NumericVector den(1, 0.0);
	NumericVector chisamp(1, 0.0);
	for(int j = 0; j < m_samps; j++){
		Z2[0] = R::rchisq(1.0);
		den[0]=0.0;
		for(int i = 0; i < L; i++){
			chisamp[0] = R::rchisq(r[i]);
			den[0] = den[0] + (lambda[i]*ratio[0]+1.0)*chisamp[0];
		}
		MC[0] = Z2[0]/den[0];
		MCt[j] = MC[0]*(c1t[0]*ratio[0] + c2t[0]);
		MCn[j] = MC[0]*(c1n[0]*ratio[0] + c2n[0]);
		MCe[j] = MC[0]*(c1e[0]*ratio[0] + c2e[0]);
	}
	
	NumericVector Ft(m_the, 0.0); 
	NumericVector Fn(m_the, 0.0); 
	NumericVector Fe(m_the, 0.0);	
	for(int i = 0; i < m_the; i++){
		for(int j = 0; j < m_samps; j++){
			if(MCt[j] < plausseq[i]){
				Ft[i] = Ft[i] + (1.0 / m_samps);	
			}
			if(MCn[j] < plausseq[i]){
				Fn[i] = Fn[i] + (1.0 / m_samps);	
			}
			if(MCe[j] < plausseq[i]){
				Fe[i] = Fe[i] + (1.0 / m_samps);	
			}
		}
	}
	
	NumericVector plaus_t(m_the, 0.0); 
	NumericVector plaus_n(m_the, 0.0); 
	NumericVector plaus_e(m_the, 0.0);	
	for(int i = 0; i < m_the; i++){
		plaus_t[i] = 1.0 - Ft[i];
		plaus_n[i] = 1.0 - Fn[i];
		plaus_e[i] = 1.0 - Fe[i];
	}

	
	result = Rcpp::List::create(Rcpp::Named("plauses.theta") = plaus_t, Rcpp::Named("plauses.new") = plaus_n, Rcpp::Named("plauses.exs") = plaus_e);
	return result;
	
}

Rcpp::List plaus_two_stage_full(NumericVector theta, NumericVector xBy, NumericVector S, NumericVector lambda,  NumericVector r, NumericVector csigma, NumericVector ratio){

	List result;
	int L = S.length();
	int m_the = theta.length();
	int m_samps = 10000;
	
	NumericVector MC(1, 0.0);
	NumericVector MCt(m_samps, 0.0);
	NumericVector MCn(m_samps, 0.0);
	NumericVector Z2(1, 0.0);
	NumericVector plausseq(m_the, 0.0);
	
	NumericVector sumS(1, 0.0);
	for(int j = 0; j < L; j++){
		sumS[0] = sumS[0] + S[j];
	}	

	for(int j = 0; j < m_the; j++){
		plausseq[j] = (theta[j] - xBy[0])*(theta[j] - xBy[0])/(sumS[0]);
	}	
	NumericVector den(1, 0.0);
	NumericVector chisamp(1, 0.0);
	for(int j = 0; j < m_samps; j++){
		Z2[0] = R::rchisq(1.0);
		den[0]=0.0;
		for(int i = 0; i < L; i++){
			chisamp[0] = R::rchisq(r[i]);
			den[0] = den[0] + (lambda[i]*ratio[0]+1.0)*chisamp[0];
		}
		MC[0] = Z2[0]/den[0];
		MCt[j] = MC[0]*csigma[0];
		MCn[j] = MC[0]*(csigma[0]+1.0);
	}
	
	NumericVector Ft(m_the, 0.0); 
	NumericVector Fn(m_the, 0.0); 	
	for(int i = 0; i < m_the; i++){
		for(int j = 0; j < m_samps; j++){
			if(MCt[j] < plausseq[i]){
				Ft[i] = Ft[i] + (1.0 / m_samps);	
			}
			if(MCn[j] < plausseq[i]){
				Fn[i] = Fn[i] + (1.0 / m_samps);	
			}
		}
	}
	
	NumericVector plaus_t(m_the, 0.0); 
	NumericVector plaus_n(m_the, 0.0); 	
	for(int i = 0; i < m_the; i++){
		plaus_t[i] = 1.0 - Ft[i];
		plaus_n[i] = 1.0 - Fn[i];
	}

	
	result = Rcpp::List::create(Rcpp::Named("plauses.theta") = plaus_t, Rcpp::Named("plauses.new") = plaus_n);
	return result;
	
}



Rcpp::List plaus_two_stage(NumericVector theta, NumericVector xBy, NumericVector S, NumericVector lambda, NumericMatrix auxiliary, NumericVector csigma, NumericVector s2a, NumericVector s2e, NumericVector assoc){

	List result;
	int L = S.length();
	int m_the = theta.length();
	int m_samps = auxiliary.nrow();
	
	NumericVector MC(1, 0.0);
	NumericVector MCt(m_samps, 0.0);
	NumericVector MCn(m_samps, 0.0);
	NumericVector Z2(1, 0.0);
	NumericVector plausseq(m_the, 0.0);	

	if(assoc[0]==0.0){

		NumericVector prodS(1, 1.0);
		NumericVector prodvar(1, 1.0);
		for(int j = 0; j < (L-1); j++){
			prodS[0] = prodS[0]*S[j];
			prodvar[0] = prodvar[0]*(lambda[j]*s2a[0] + s2e[0]); 
		}	
		prodS[0] = prodS[0] + S[L-1];	
		
		for(int j = 0; j < m_the; j++){
			plausseq[j] = (theta[j] - xBy[0])*(theta[j] - xBy[0])/(prodS[0]);
		}	
		
		for(int j = 0; j < m_samps; j++){
			Z2[0] = R::rchisq(1.0);
			MC[0] = Z2[0]/(prodvar[0]*std::exp(auxiliary(j,0)) + s2e[0]*std::exp(auxiliary(j,1)));
			MCt[j] = MC[0]*csigma[0];
			MCn[j] = MC[0]*(csigma[0]+s2e[0]);
		}
		
	}else {
		NumericVector slogS(1, 0.0);
		NumericVector slogvar(1, 0.0);
		for(int j = 0; j < (L-1); j++){
			slogS[0] = slogS[0] + std::log(S[j]);
			slogvar[0] = slogvar[0]+std::log(lambda[j]*s2a[0] + s2e[0]); 
		}	
		slogS[0] = slogS[0] + std::log(S[L-1]);		
		
		for(int j = 0; j < m_the; j++){
			plausseq[j] = (theta[j] - xBy[0])/slogS[0];
		}	
		
		for(int j = 0; j < m_samps; j++){
			Z2[0] = R::rnorm(0.0,1.0);
			MC[0] = Z2[0]/(slogvar[0] + auxiliary(j,0) + std::log(s2e[0]) + auxiliary(j,1));	
			MCt[j] = MC[0]*std::sqrt(csigma[0]);
			MCn[j] = MC[0]*std::sqrt(csigma[0]+s2e[0]);
		}
	}
	
	
	NumericVector Ft(m_the, 0.0); 
	NumericVector Fn(m_the, 0.0); 
	for(int i = 0; i < m_the; i++){
		for(int j = 0; j < m_samps; j++){
			if(MCt[j] < plausseq[i]){
				Ft[i] = Ft[i] + (1.0 / m_samps);	
			}
			if(MCn[j] < plausseq[i]){
				Fn[i] = Fn[i] + (1.0 / m_samps);	
			}
		}
	}
	
	NumericVector plaus_t(m_the, 0.0); 
	NumericVector plaus_n(m_the, 0.0); 

	if(assoc[0] == 0.0){
		for(int i = 0; i < m_the; i++){
			plaus_t[i] = 1.0 - Ft[i];
			plaus_n[i] = 1.0 - Fn[i];
		}
	}else {
		for(int i = 0; i < m_the; i++){
			plaus_t[i] = 1.0 - std::abs(2.0*Ft[i]-1.0);
			plaus_n[i] = 1.0 - std::abs(2.0*Fn[i]-1.0);
		}
	}
	
	result = Rcpp::List::create(Rcpp::Named("plauses.theta") = plaus_t, Rcpp::Named("plauses.new") = plaus_n);
	return result;
	
}





Rcpp::List plaus_balanced(NumericVector thetaseq, NumericVector saseq, NumericVector seseq, NumericVector n, NumericVector n_i, NumericVector S, NumericVector lambda, NumericVector r, NumericVector Ybar, NumericVector numk){
	
	List result;
	int m_the = thetaseq.length();
	int m_se = seseq.length();
	int m_sa = saseq.length();
	int dn_i = n_i.length();

	NumericVector sumn_i2(1, 0.0);
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}

	
	NumericVector U1(10000,0.0); NumericVector U2(10000,0.0); NumericVector U3(10000,0.0);	
	NumericVector H(10000,0.0);
	U1 = Rcpp::runif(10000,0.0,1.0); U2 = Rcpp::runif(10000,0.0,1.0); U3 = Rcpp::runif(10000,0.0,1.0);
	for(int j = 0; j < 10000; j++){
		H[j] = (1-std::abs(2*U1[j] - 1))*(1-std::abs(2*U2[j] - 1))*(1-std::abs(2*U3[j] - 1));	
	}
	
	NumericVector F_sase(m_sa*m_se, 0.0);
	int index = 0;
	for(int j = 0; j < m_sa; j++){
		for(int i = 0; i < m_se; i++){	
			F_sase[index] = (1.0-std::abs(2.0*R::pchisq(S[0]/(lambda[0]*saseq[j] + seseq[i]), r[0], 1, 0) - 1.0))*(1.0-std::abs(2.0*R::pchisq(S[1]/(lambda[1]*saseq[j] + seseq[i]), r[1], 1, 0) - 1.0));
			index = index+1;	
		}
	}
	
	NumericVector plausestheta(m_the, 0.0);
	NumericVector plaus(1, 0.0);
	NumericVector plausesystar(m_the, 0.0);
	NumericVector plausystar(1, 0.0);
	NumericVector plausesystarexs(m_the, 0.0);
	NumericVector plausystarexs(1, 0.0);
	NumericVector F_the(1, 0.0);
	NumericVector F_ystark(1, 0.0);
	NumericVector F_ystark_exs(1, 0.0);
	
	for(int k = 0; k < m_the; k++){
		index = 0;
		for(int j = 0; j < m_sa; j++){
			for(int i = 0; i < m_se; i++){	
				F_the[0] = (1-std::abs(2.0*R::pnorm((Ybar[0] - thetaseq[k])/std::sqrt(saseq[j]*(1+(1/(n[0]*n[0]))*sumn_i2[0]) + seseq[i]/n[0]),0.0,1.0,1,0)-1))*F_sase[index];
				F_ystark[0] = (1-std::abs(2.0*R::pnorm((Ybar[0] - thetaseq[k])/std::sqrt(saseq[j]*(1+(1/(n[0]*n[0]))*sumn_i2[0]) + seseq[i]*(1.0/n[0] + 1.0/numk[0])),0.0,1.0,1,0)-1))*F_sase[index];
				F_ystark_exs[0] = (1-std::abs(2.0*R::pnorm((Ybar[0] - thetaseq[k])/std::sqrt(saseq[j]*(1-2*n_i[dn_i-1]/n[0]+(1/(n[0]*n[0]))*sumn_i2[0]) + seseq[i]*(1.0/n[0] + 1.0/numk[0])),0.0,1.0,1,0)-1))*F_sase[index];
				index = index+1;
				plaus[0] = 0.0;plausystar[0] = 0.0;plausystarexs[0] = 0.0;
				for(int h = 0; h <  10000; h++){
					if(H[h] <= F_the[0]){
						plaus[0] = plaus[0] + 1/10000.0;
					}
					if(H[h] <= F_ystark[0]){
						plausystar[0] = plausystar[0] + 1/10000.0;
					}
					if(H[h] <= F_ystark_exs[0]){
						plausystarexs[0] = plausystarexs[0] + 1/10000.0;
					}
				}
				plausestheta[k] = std::max(plausestheta[k], plaus[0]);
				plausesystar[k] = std::max(plausesystar[k], plausystar[0]);
				plausesystarexs[k] = std::max(plausesystarexs[k], plausystarexs[0]);
			}
		}
	}	
		
		
	NumericVector maxplaus(1, 0.0);
	NumericVector maxplausstar(1, 0.0);
	NumericVector maxplausstarexs(1, 0.0);
	for(int k = 0; k < m_the; k++){
		maxplaus[0] = std::max(maxplaus[0], plausestheta[k]);
		maxplausstar[0] = std::max(maxplausstar[0], plausesystar[k]);
		maxplausstarexs[0] = std::max(maxplausstarexs[0], plausesystarexs[k]);
	}
	
	for(int k = 0; k < m_the; k++){
		plausestheta[k] = plausestheta[k]/maxplaus[0];
		plausesystar[k] = plausesystar[k]/maxplausstar[0];
		plausesystarexs[k] = plausesystarexs[k]/maxplausstarexs[0];
	}
	
	result = Rcpp::List::create(Rcpp::Named("plausestheta") = plausestheta, Rcpp::Named("plausesnew") = plausesystar, Rcpp::Named("plausesexs") = plausesystarexs);
	return result;
	
}

Rcpp::List plaus_unbalanced(NumericVector thetaseq, NumericVector saseq, NumericVector seseq, NumericVector n, NumericVector n_i, NumericVector S, NumericVector lambda, NumericVector r, NumericVector Ybar, NumericVector numk, NumericVector samples1, NumericVector samples2){
	
	List result;
	int m_the = thetaseq.length();
	int m_se = seseq.length();
	int m_sa = saseq.length();
	int dn_i = n_i.length();
	int m_samples = samples1.length();
	int L = S.length();
	
	NumericVector sumn_i2(1, 0.0);
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}

	
	NumericVector U1(10000,0.0); NumericVector U2(10000,0.0); NumericVector U3(10000,0.0);	
	NumericVector H(10000,0.0);
	U1 = Rcpp::runif(10000,0.0,1.0); U2 = Rcpp::runif(10000,0.0,1.0); U3 = Rcpp::runif(10000,0.0,1.0);
	for(int j = 0; j < 10000; j++){
		H[j] = (1-std::abs(2*U1[j] - 1))*(1-std::abs(2*U2[j] - 1))*(1-std::abs(2*U3[j] - 1));	
	}
	
	NumericVector F_sase(m_sa*m_se, 0.0);
	NumericVector F_sa(1, 0.0);
	NumericVector F_se(1, 0.0);
	NumericVector omega1(1, 0.0);
	NumericVector omega2(1, 0.0);
	
	
	int index = 0;
	for(int j = 0; j < m_sa; j++){
		for(int i = 0; i < m_se; i++){
			omega1[0] = 0.0; omega2[0] = 0.0; F_sa[0] = 0.0; F_se[0] = 0.0;
			for(int l = 0; l < (L-1); l++){
				omega1[0] = omega1[0] + std::log(S[l]) - std::log(lambda[l]*saseq[j] + seseq[i]);	
			}
			omega2[0] = std::log(S[L-1]) - std::log(lambda[L-1]*saseq[j] + seseq[i]); 
			for(int h = 0; h < m_samples; h++){
				if(samples1[h] < omega1[0]){
					F_sa[0] = F_sa[0] + 1.0/10000.0;
				}
				if(samples2[h] < omega2[0]){
					F_se[0] = F_se[0] + 1.0/10000.0;
				}
			}
			F_sase[index] = (1.0-std::abs(2.0*F_sa[0] - 1.0))*(1.0-std::abs(2.0*F_se[0] - 1.0));
			index = index+1;	
		}
	}
	
	NumericVector plausestheta(m_the, 0.0);
	NumericVector plaus(1, 0.0);
	NumericVector plausesystar(m_the, 0.0);
	NumericVector plausystar(1, 0.0);
	NumericVector plausesystarexs(m_the, 0.0);
	NumericVector plausystarexs(1, 0.0);
	NumericVector F_the(1, 0.0);
	NumericVector F_ystark(1, 0.0);
	NumericVector F_ystark_exs(1, 0.0);
	
	for(int k = 0; k < m_the; k++){
		index = 0;
		for(int j = 0; j < m_sa; j++){
			for(int i = 0; i < m_se; i++){	
				F_the[0] = (1-std::abs(2.0*R::pnorm((Ybar[0] - thetaseq[k])/std::sqrt(saseq[j]*(1+(1/(n[0]*n[0]))*sumn_i2[0]) + seseq[i]/n[0]),0.0,1.0,1,0)-1))*F_sase[index];
				F_ystark[0] = (1-std::abs(2.0*R::pnorm((Ybar[0] - thetaseq[k])/std::sqrt(saseq[j]*(1+(1/(n[0]*n[0]))*sumn_i2[0]) + seseq[i]*(1.0/n[0] + 1.0/numk[0])),0.0,1.0,1,0)-1))*F_sase[index];
				F_ystark_exs[0] = (1-std::abs(2.0*R::pnorm((Ybar[0] - thetaseq[k])/std::sqrt(saseq[j]*(1-2*n_i[dn_i-1]/n[0]+(1/(n[0]*n[0]))*sumn_i2[0]) + seseq[i]*(1.0/n[0] + 1.0/numk[0])),0.0,1.0,1,0)-1))*F_sase[index];
				index = index+1;
				plaus[0] = 0.0;plausystar[0] = 0.0;plausystarexs[0] = 0.0;
				for(int h = 0; h <  10000; h++){
					if(H[h] <= F_the[0]){
						plaus[0] = plaus[0] + 1/10000.0;
					}
					if(H[h] <= F_ystark[0]){
						plausystar[0] = plausystar[0] + 1/10000.0;
					}
					if(H[h] <= F_ystark_exs[0]){
						plausystarexs[0] = plausystarexs[0] + 1/10000.0;
					}
				}
				plausestheta[k] = std::max(plausestheta[k], plaus[0]);
				plausesystar[k] = std::max(plausesystar[k], plausystar[0]);
				plausesystarexs[k] = std::max(plausesystarexs[k], plausystarexs[0]);
			}
		}
	}	
		

	result = Rcpp::List::create(Rcpp::Named("omega1") = omega1,Rcpp::Named("omega2") = omega2,Rcpp::Named("plausestheta") = plausestheta, Rcpp::Named("plausesnew") = plausesystar, Rcpp::Named("plausesexs") = plausesystarexs);
	return result;
	
}

Rcpp::List plaus_balanced_marginal(NumericVector thetaseq, NumericVector n, NumericVector n_i, NumericVector S, NumericVector lambda, NumericVector r, NumericVector Ybar, NumericVector numk, NumericVector sa2, NumericVector se2){
	
	List result;
	int m_the = thetaseq.length();
	int dn_i = n_i.length();

	NumericVector sumn_i2(1, 0.0);
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
	
	NumericVector H(10000,0.0);
	H = Rcpp::runif(10000,0.0,1.0); 
	
	
	NumericVector plausestheta(m_the, 0.0);
	NumericVector plausesystar(m_the, 0.0);
	NumericVector plausesystarexs(m_the, 0.0);
	NumericVector sd_the(10000, 0.0);
	NumericVector sd_ystark(10000, 0.0);
	NumericVector sd_ystark_exs(10000, 0.0);
	NumericVector F_the(1, 0.0);
	NumericVector F_ystark(1, 0.0);
	NumericVector F_ystark_exs(1, 0.0);
	NumericVector Z(1, 0.0);
	
	for(int j = 0; j < 10000; j++){
		Z[0] = R::rnorm(0.0, 1.0);
		sd_the[j] = Z[0]*std::sqrt(sa2[j]*(1+(1/(n[0]*n[0]))*sumn_i2[0]) + se2[j]/n[0]);
		sd_ystark[j] = Z[0]*std::sqrt(sa2[j]*(1+(1/(n[0]*n[0]))*sumn_i2[0]) + se2[j]*(1.0/numk[0] + 1.0/n[0]));		
		sd_ystark_exs[j] = Z[0]*std::sqrt(sa2[j]*(1-2.0*n_i[dn_i-1]/n[0]+(1/(n[0]*n[0]))*sumn_i2[0]) + se2[j]*(1.0/numk[0] + 1.0/n[0]));	
	}
	
	
	for(int j = 0; j < m_the; j++){
		F_the[0] = 0.0;F_ystark[0] = 0.0;F_ystark_exs[0] = 0.0;
		for(int k = 0; k < 10000; k++){
			if(sd_the[k] < (Ybar[0] - thetaseq[j])){
				F_the[0] = F_the[0] + 1.0/10000.0;	
			}
			if(sd_ystark[k] < (Ybar[0] - thetaseq[j])){
				F_ystark[0] = F_ystark[0] + 1.0/10000.0;	
			}
			if(sd_ystark_exs[k] < (Ybar[0] - thetaseq[j])){
				F_ystark_exs[0] = F_ystark_exs[0] + 1.0/10000.0;	
			}
		}
		for(int k = 0; k < 10000; k++){
			if(H[k] <= (1.0-std::abs(2.0*F_the[0] - 1))){
				plausestheta[j] = plausestheta[j] + 1.0/10000.0;	
			}
			if(H[k] <= (1.0-std::abs(2.0*F_ystark[0] - 1))){
				plausesystar[j] = plausesystar[j] + 1.0/10000.0;	
			}
			if(H[k] <= (1.0-std::abs(2.0*F_ystark_exs[0] - 1))){
				plausesystarexs[j] = plausesystarexs[j] + 1.0/10000.0;	
			}
		}
	}
	
	NumericVector maxplaus(1, 0.0);
	NumericVector maxplausstar(1, 0.0);
	NumericVector maxplausstarexs(1, 0.0);
	for(int k = 0; k < m_the; k++){
		maxplaus[0] = std::max(maxplaus[0], plausestheta[k]);
		maxplausstar[0] = std::max(maxplausstar[0], plausesystar[k]);
		maxplausstarexs[0] = std::max(maxplausstarexs[0], plausesystarexs[k]);
	}
	
	for(int k = 0; k < m_the; k++){
		plausestheta[k] = plausestheta[k]/maxplaus[0];
		plausesystar[k] = plausesystar[k]/maxplausstar[0];
		plausesystarexs[k] = plausesystarexs[k]/maxplausstarexs[0];
	}
	
	result = Rcpp::List::create(Rcpp::Named("plausestheta") = plausestheta, Rcpp::Named("plausesnew") = plausesystar, Rcpp::Named("plausesexs") = plausesystarexs);
	return result;
	
}
						   
Rcpp::List plaus_unbalanced_marginal(NumericVector thetaseq, NumericVector n, NumericVector n_i, NumericVector S, NumericVector lambda, NumericVector r, NumericVector Ybar, NumericVector numk, NumericVector sa2, NumericVector se2){
	
	List result;
	int m_the = thetaseq.length();
	int dn_i = n_i.length();

	NumericVector sumn_i2(1, 0.0);
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
	
	NumericVector H(10000,0.0);
	H = Rcpp::runif(10000,0.0,1.0); 
	
	
	NumericVector plausestheta(m_the, 0.0);
	NumericVector plausesystar(m_the, 0.0);
	NumericVector plausesystarexs(m_the, 0.0);
	NumericVector sd_the(10000, 0.0);
	NumericVector sd_ystark(10000, 0.0);
	NumericVector sd_ystark_exs(10000, 0.0);
	NumericVector F_the(1, 0.0);
	NumericVector F_ystark(1, 0.0);
	NumericVector F_ystark_exs(1, 0.0);
	NumericVector Z(1, 0.0);
	
	for(int j = 0; j < 10000; j++){
		Z[0] = R::rnorm(0.0, 1.0);
		sd_the[j] = Z[0]*std::sqrt(sa2[j]*(1+(1/(n[0]*n[0]))*sumn_i2[0]) + se2[j]/n[0]);
		sd_ystark[j] = Z[0]*std::sqrt(sa2[j]*(1+(1/(n[0]*n[0]))*sumn_i2[0]) + se2[j]*(1.0/numk[0] + 1.0/n[0]));		
		sd_ystark_exs[j] = Z[0]*std::sqrt(sa2[j]*(1-2.0*n_i[dn_i-1]/n[0]+(1/(n[0]*n[0]))*sumn_i2[0]) + se2[j]*(1.0/numk[0] + 1.0/n[0]));	
	}
	
	
	for(int j = 0; j < m_the; j++){
		F_the[0] = 0.0;F_ystark[0] = 0.0;F_ystark_exs[0] = 0.0;
		for(int k = 0; k < 10000; k++){
			if(sd_the[k] < Ybar[0] - thetaseq[j]){
				F_the[0] = F_the[0] + 1.0/10000.0;	
			}
			if(sd_ystark[k] < Ybar[0] - thetaseq[j]){
				F_ystark[0] = F_ystark[0] + 1.0/10000.0;	
			}
			if(sd_ystark_exs[k] < Ybar[0] - thetaseq[j]){
				F_ystark_exs[0] = F_ystark_exs[0] + 1.0/10000.0;	
			}
		}
		for(int k = 0; k < 10000; k++){
			if(H[k] <= (1.0-std::abs(2.0*F_the[0] - 1))){
				plausestheta[j] = plausestheta[j] + 1.0/10000.0;	
			}
			if(H[k] <= (1.0-std::abs(2.0*F_ystark[0] - 1))){
				plausesystar[j] = plausesystar[j] + 1.0/10000.0;	
			}
			if(H[k] <= (1.0-std::abs(2.0*F_ystark_exs[0] - 1))){
				plausesystarexs[j] = plausesystarexs[j] + 1.0/10000.0;	
			}
		}
	}
	
	result = Rcpp::List::create(Rcpp::Named("plausestheta") = plausestheta, Rcpp::Named("plausesnew") = plausesystar, Rcpp::Named("plausesexs") = plausesystarexs);
	return result;
	
}						   


Rcpp::List plaus_marginal_lmer(NumericVector thetaseq, NumericVector total_sigma, NumericVector xBy){

	List result;
	int m_the = thetaseq.length();
	NumericVector plausestheta(m_the,0.0);
	NumericVector F_the(1,0.0);

	for(int j = 0; j < m_the; j++){
		F_the[0] = 0.0;
		for(int k = 0; k < 10000; k++){
			if(total_sigma[k] < (xBy[0] - thetaseq[j])){
				F_the[0] = F_the[0] + 1.0/10000.0;	
			}
		}
		F_the[0] = (1.0-std::abs(2.0*F_the[0] - 1.0));
		plausestheta[j] = F_the[0];	
	}
	
	result = Rcpp::List::create(Rcpp::Named("plausestheta") = plausestheta);
	return result;
	
}						   




Rcpp::List randsetspredlmer(NumericMatrix S, NumericVector dimS, NumericVector U, NumericMatrix C1, NumericMatrix C2, NumericVector By, NumericVector x, NumericVector ztz) {
	
	List result;
	
	arma::colvec xa = as<arma::colvec>(x);
	arma::colvec Cxa; Cxa.zeros(2);
	arma::colvec Bya = as<arma::colvec>(By);
	double xBy = dot(xa, Bya);
	arma::mat Csigma2;  Csigma2.zeros(2,2);
	NumericVector zeroes4(4,0.0);
	NumericMatrix Csigma(2,2,zeroes4.begin());
	NumericVector total_sigma(1,0.0);
	NumericVector Uu(1, 0.0); NumericVector Ul(1, 0.0);
	NumericVector zeroes(20000,0.0);
	NumericMatrix randsetpred = NumericMatrix(10000, 2, zeroes.begin());
	NumericVector thetas(10000, 0.0);
	NumericVector Z(1, 0.0);
		
	for(int j=0; j < 10000; j++){
		Z[0] = R::rnorm(0.0,1.0);
		for(int r = 0; r < 2; r++){
			for(int s = 0; s < 2; s++){
				Csigma(r,s) = S(j,1)*C1(r,s) + S(j,0)*C2(r,s);
			}
		}
		Csigma2 = as<arma::mat>(Csigma);
		Cxa = Csigma2*xa;
		total_sigma[0] = dot(xa, Cxa) + S(j,0)*ztz[0];
		total_sigma[0] = std::sqrt(total_sigma[0]);
		thetas[j] = Z[0]*total_sigma[0] + xBy;
	}

	std::sort(thetas.begin(), thetas.end());
	
	for(int j=0; j < 10000; j++){
		Ul[0] = 0.5-std::fabs(U[j]-0.5);
		Uu[0] = 1.0-Ul[0];
		randsetpred(j,0) = thetas[std::round(10000*Ul[0])];
		randsetpred(j,1) = thetas[std::round(10000*Uu[0])];
	}
	
	result = Rcpp::List::create(Rcpp::Named("randsetpred") = randsetpred);
	
	return result;
}



Rcpp::List randsetspred(NumericMatrix S, NumericVector dimS, NumericVector nsize, NumericVector n_i, NumericVector dimn_i, NumericVector k, NumericVector U, NumericVector Ybar) {
	
	List result;
	int M = int(dimS[0]);
	int n = int(nsize[0]);
	int index = int(1);
	int dn_i = int(dimn_i[0]);
	NumericVector sumn_i2(1,0.0);
	NumericVector Qsampsw(10000,0.0);
	NumericVector Qsampsn(10000,0.0);
	NumericVector QsampsT(10000,0.0);
	NumericVector Qwl(1,0.0);
	NumericVector Qwu(1,0.0);
	NumericVector Qnl(1,0.0);
	NumericVector Qnu(1,0.0);
	NumericVector QTl(1,0.0);
	NumericVector QTu(1,0.0);
	NumericVector Ul(1,0.0);
	NumericVector Uu(1,0.0);
	NumericVector Z = NumericVector(1, 0.0); 
	NumericVector zeroes = NumericVector(10000*6, 0.0); 
        NumericMatrix randsetpred = NumericMatrix(10000, 6, zeroes.begin());

	
	/*for(int j=0; j < M; j++){
		Sa[j] = S(j,0); Se[j] = S(j,1);	
	}
	std::sort(Sa.begin(), Sa.end()); std::sort(Se.begin(), Se.end());
	*/
	
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
	
	
	/*if(Sa[0] < 0.0){	
	for(int j=0; j < M; j++){
		U3[0] = 0.0;
		U3[0] = std::max(R::runif(0.0,1.0), U3[0]);U3[0] = std::max(R::runif(0.0,1.0), U3[0]);U3[0] = std::max(R::runif(0.0,1.0), U3[0]);
		U3l[0] = 0.5 - std::abs(0.5 - U3[0]);
		U3u[0] = 1.0 -U3l[0];
		siga[0] = std::max(Sa[round(U3u[0]*M)],0.0);
		sige[0] = std::max(Se[round(U3u[0]*M)],0.0);
		Zint[0] = R::qnorm(U3l[0],0.0,1.0,1,0); Zint[1] = R::qnorm(U3u[0],0.0,1.0,1,0);
		Qwl[0] = Zint[0]*std::sqrt(siga[0]*(1-(2*n_i[dn_i-1]/n)+(1/(n*n))*sumn_i2[0])+sige[0]*((1/n)+1/(k[0])));
		Qwu[0] = Zint[1]*std::sqrt(siga[0]*(1-(2*n_i[dn_i-1]/n)+(1/(n*n))*sumn_i2[0])+sige[0]*((1/n)+1/(k[0])));
		Qnl[0] = Zint[0]*std::sqrt(siga[0]*(1+(1/(n*n))*sumn_i2[0])+sige[0]*((1/n)+1/(k[0])));
		Qnu[0] = Zint[1]*std::sqrt(siga[0]*(1+(1/(n*n))*sumn_i2[0])+sige[0]*((1/n)+1/(k[0])));
		QTl[0] = Zint[0]*std::sqrt(siga[0]*(1+(1/(n*n))*sumn_i2[0])+sige[0]*((1/n)));
		QTu[0] = Zint[1]*std::sqrt(siga[0]*(1+(1/(n*n))*sumn_i2[0])+sige[0]*((1/n)));
		randsetpred(j,0) = Ybar[0]+Qwl[0];
		randsetpred(j,1) = Ybar[0]+Qwu[0];
		randsetpred(j,2) = Ybar[0]+Qnl[0];
		randsetpred(j,3) = Ybar[0]+Qnu[0];
		randsetpred(j,4) = Ybar[0]+QTl[0];
		randsetpred(j,5) = Ybar[0]+QTu[0];
	}
	}else {
*/		
	
	for(int j=0; j < 10000; j++){
		Z[0] = R::rnorm(0.0,1.0);
		index = round(M*R::runif(0.0,1.0));
		Qsampsw[j] = Z[0]*std::sqrt(S(index,0)*(1-(2*n_i[dn_i-1]/n)+(1/(n*n))*sumn_i2[0])+S(index,1)*((1/n)+1/(k[0])));
		Qsampsn[j] = Z[0]*std::sqrt(S(index,0)*(1+(1/(n*n))*sumn_i2[0])+S(index,1)*((1/n)+1/(k[0])));
		QsampsT[j] = Z[0]*std::sqrt(S(index,0)*(1+(1/(n*n))*sumn_i2[0])+S(index,1)*(1/n));
	}


	std::sort(Qsampsw.begin(), Qsampsw.end());
	std::sort(Qsampsn.begin(), Qsampsn.end());
	std::sort(QsampsT.begin(), QsampsT.end());
	
	for(int j=0; j < 10000; j++){
		Ul[0] = 0.5-std::fabs(U[j]-0.5);
		Uu[0] = 1.0-Ul[0];
		Qwl[0] = Qsampsw[std::round(10000*Ul[0])];
		Qwu[0] = Qsampsw[std::round(10000*Uu[0])];
		Qnl[0] = Qsampsn[std::round(10000*Ul[0])];
		Qnu[0] = Qsampsn[std::round(10000*Uu[0])];
		QTl[0] = QsampsT[std::round(10000*Ul[0])];
		QTu[0] = QsampsT[std::round(10000*Uu[0])];
		randsetpred(j,0) = Ybar[0]+Qwl[0];
		randsetpred(j,1) = Ybar[0]+Qwu[0];
		randsetpred(j,2) = Ybar[0]+Qnl[0];
		randsetpred(j,3) = Ybar[0]+Qnu[0];
		randsetpred(j,4) = Ybar[0]+QTl[0];
		randsetpred(j,5) = Ybar[0]+QTu[0];

	}
//	}

result = Rcpp::List::create(Rcpp::Named("randsetpred") = randsetpred);

	return result;
	
	
}


Rcpp::List genIM(NumericVector Y, NumericMatrix Z, NumericVector museq, NumericVector saseq, NumericVector seseq, NumericVector M) {
	
	List result;
	int n = Y.length();
	int m = round(M[0]);
	int s_par = museq.length(); 
	
	double data_liks[s_par][s_par][s_par];
	double max_data_liks = -10000000.0;
	NumericVector max_sim_liks(m,-100000000.0);
	
	NumericVector U(n,0.0); 
	NumericVector lik(1, 0.0);
	arma::vec z; z.zeros(n);
	arma::vec ym; ym.zeros(n);arma::vec ym2; ym2.zeros(n);arma::vec ymsim; ymsim.zeros(n);
	arma::mat ZZ = as<arma::mat>(Z); 
	arma::mat Sigma; Sigma.zeros(n,n);
	arma::mat chSigma; arma::mat tmp; arma::mat rss; arma::mat tmpsim; arma::mat rsssim; 
	arma::mat I_n; I_n.zeros(n,n);
	arma::mat Ud; Ud.zeros(n,m); NumericVector templik(m,0.0);arma::mat siglik; siglik.zeros(s_par,s_par); 
	arma::vec ztz; ztz.zeros(m);
	
	
	for(int i = 0; i < n; i++){
		I_n(i,i) = 1.0;
	}
	
	for(int j = 0; j < s_par; j++){
		for(int q = 0; q < n; q++){
			ym(q) = Y[q] - museq[j];	
		}
		for(int k = 0; k < s_par; k++){
			for(int t = 0; t < s_par; t++){
				for(int q = 0; q < n; q++){
					for(int r = 0; r<n; r++){
						Sigma(q,r) = ZZ(q,r)*saseq[k] + I_n(q,r)*seseq[t];	
					}
				}
				chSigma = arma::chol(Sigma);
				tmp = solve(trimatl(chSigma.t()), ym);
				rss = dot(tmp,tmp);
				lik[0] = 0.0;
				for(int q = 0; q < n; q++){
					lik[0] = lik[0] - log(chSigma(q,q));
				}
				siglik(k,t) = lik[0];
				lik[0] = lik[0] - 0.5 * n * log(2 * M_PI) - 0.5 * rss(0,0);
				data_liks[j][k][t] = lik[0];
				max_data_liks = std::max(max_data_liks, lik[0]);
			}
		}
	}
	
	
	double dataratios[s_par][s_par][s_par];
	for(int j = 0; j < s_par; j++){
		for(int k = 0; k < s_par; k++){
			for(int t = 0; t < s_par; t++){
				dataratios[j][k][t] = std::exp(data_liks[j][k][t]-max_data_liks);
			}
		}
	}
	
	// Loop for sim data
	
	for(int q = 0; q < m; q++){
		U = Rcpp::rnorm(n,0.0,1.0);
		for(int s = 0; s < n; s++){
			Ud(s,q) = U[s];
			ztz(q) = ztz(q) + U[s]*U[s];
		}			
	}
	
	double nums[s_par][s_par][s_par]; double dens[s_par][s_par][s_par]; double maxdens = -10000000.0;
	double simratios[s_par][s_par][s_par]; double plauses[s_par][s_par][s_par]; arma::vec Uu;Uu.zeros(n);
	
	for(int j = 0; j < s_par; j++){
		for(int k = 0; k < s_par; k++){
			for(int t = 0; t < s_par; t++){
				plauses[j][k][t] = 0.0;
			}
		}
	}
	
	
	for(int q = 0; q < m; q++){
		for(int s = 0; s < n; s++){
			Uu(s) = Ud(s,q);
		}
		for(int p = 0; p < s_par; p++ ){
			for(int k = 0; k < s_par; k++){
				for(int t = 0; t < s_par; t++){
					nums[p][k][t] = siglik(k,t)- 0.5 * n * log(2 * M_PI) - 0.5 * ztz(q);
					for(int s = 0; s < n; s++){
						for(int r = 0; r<n; r++){
							Sigma(s,r) = ZZ(s,r)*saseq[k] + I_n(s,r)*seseq[t];	
						}
					}
					chSigma = arma::chol(Sigma);
					ym = chSigma.t()*Uu;
					for(int r = 0; r < n; r++){
						ym(r) = ym(r) + museq[p];	
					}
					maxdens = -1000000000.0;
						for(int j = 0; j < s_par; j++){
							for(int l = 0; l < s_par; l++){
								for(int v = 0; v < s_par; v++){
									for(int s = 0; s < n; s++){
										for(int r = 0; r<n; r++){
											Sigma(s,r) = ZZ(s,r)*saseq[l] + I_n(s,r)*seseq[v];	
										}
									}
									for(int r = 0; r < n; r++){
										ym2(r) = ym(r) - museq[j];	
									}
									tmp = solve(trimatl(chSigma.t()), ym2);
									rss = dot(tmp,tmp);
									dens[j][l][v] = siglik(l,v) - 0.5 * n * log(2 * M_PI) - 0.5 * rss(0,0);
									maxdens = std::max(maxdens, dens[j][l][v]);
								}
							}
						}
					simratios[p][k][t] = std::exp(nums[p][k][t]-maxdens);
				}
			}
		}
		for(int j = 0; j < s_par; j++){
			for(int k = 0; k < s_par; k++){
				for(int t = 0; t < s_par; t++){
					if(simratios[j][k][t] <= dataratios[j][k][t]){
						plauses[j][k][t] = plauses[j][k][t] + 1.0/m; 
					}
				}
			}
		}		
	}
	
	NumericVector zeros(s_par*s_par,0.0);  NumericMatrix plauses_musa(s_par, s_par, zeros.begin());
	
	for(int j = 0; j < s_par; j++){
		for(int k = 0; k < s_par; k++){
			for(int t = 0; t < s_par; t++){
				plauses_musa(j,k) = std::max(plauses_musa(j,k), plauses[j][k][t]);
			}
		}
	}
	
	
	result = Rcpp::List::create(Rcpp::Named("plauses") = plauses_musa);
	
	return result;
	
	
}




Rcpp::List randsetspreddens(NumericVector sigsampdens, NumericVector dimS, NumericVector nsize, NumericVector n_i, NumericVector dimn_i, NumericVector k, NumericVector Ybar, NumericVector predgrid, NumericVector dim_predgrid, NumericVector localpt, NumericVector logdenslocalpt) {
	
	List result;
	int M = int(dimS[0]);
	int n = int(nsize[0]);
	int dn_i = int(dimn_i[0]);
	int L = int(dim_predgrid[0]);
	NumericVector sumn_i2(1,0.0);
	NumericVector solnw(L, 0.0);
	NumericVector solnn(L, 0.0);
	NumericVector solnT(L, 0.0);
	NumericVector denssolnw(L, 0.0);
	NumericVector denssolnn(L, 0.0);
	NumericVector denssolnT(L, 0.0);
	NumericVector denssamp(M, 0.0);
	NumericVector Z(1, 0.0);
	NumericVector localplausesw(L, 0.0);NumericVector localplausesn(L, 0.0);NumericVector localplausesT(L, 0.0);
	
	for(int j=0; j<dn_i; j++){
		sumn_i2[0] = sumn_i2[0] + n_i[j]*n_i[j];	
	}
	
	for(int j=0; j<L; j++){
		solnw[j] = (predgrid[j]-Ybar[0])/std::sqrt(localpt[0]*(1-(2*n_i[dn_i-1]/n)+(1/(n*n))*sumn_i2[0])+localpt[1]*((1/n)+1/(k[0])));
		denssolnw[j] = logdenslocalpt[0] + R::dnorm(solnw[j], 0.0,1.0,1);
		solnn[j] = (predgrid[j]-Ybar[0])/std::sqrt(localpt[0]*(1+(1/(n*n))*sumn_i2[0])+localpt[1]*((1/n)+1/(k[0])));
		denssolnn[j] = logdenslocalpt[0] + R::dnorm(solnn[j], 0.0,1.0,1);
		solnT[j] = (predgrid[j]-Ybar[0])/std::sqrt(localpt[0]*(1+(1/(n*n))*sumn_i2[0])+localpt[1]*(1/n));
		denssolnT[j] = logdenslocalpt[0] + R::dnorm(solnT[j], 0.0,1.0,1);
	}
	
	for(int j=0; j<M; j++){
		Z[0] = R::rnorm(0.0,1.0);
		denssamp[j] = sigsampdens[j] + R::dnorm(Z[0], 0.0, 1.0, 1);
		for(int k = 0; k<L; k++){
			if(denssolnw[k] > denssamp[j]){
				localplausesw[k] = localplausesw[k]+(1.0/M); 	
			}
			if(denssolnn[k] > denssamp[j]){
				localplausesn[k] = localplausesn[k]+(1.0/M); 	
			}
			if(denssolnT[k] > denssamp[j]){
				localplausesT[k] = localplausesT[k]+(1.0/M); 	
			}
		}
	}
	
	

result = Rcpp::List::create(Rcpp::Named("localplausesw") = localplausesw,Rcpp::Named("localplausesn") = localplausesn, Rcpp::Named("localplausesT") = localplausesT);

	return result;
	
	
}








//NumericVector(*f)(NumericVector x, NumericVector uu, NumericVector vv, NumericVector yy, NumericVector zz)
Rcpp::NumericVector zeroin(NumericVector ax, NumericVector bx, NumericVector u, NumericVector v, NumericVector y, NumericVector z, Function f , NumericVector tol) {
    // code here
/*
double zeroin(ax,bx,f,tol)		An estimate to the root	
Xdouble ax;				Left border | of the range	
Xdouble bx;  				Right border| the root is seeked
Xdouble (*f)(double x);			Function under investigation	
Xdouble tol;				Acceptable tolerance		
	*/

NumericVector a(1,0.0); NumericVector b(1,0.0); NumericVector c(1,0.0); NumericVector fa(1,0.0); NumericVector fb(1,0.0); NumericVector fc(1,0.0);


  a[0] = ax[0];  b[0] = bx[0];  fa = f(a, u, v, y, z);  fb = f(b, u, v, y, z);
  c[0] = a[0];   fc[0] = fa[0];

  for(;;)		/* Main iteration loop	*/
  {
    double prev_step = b[0]-a[0];		/* Distance from the last but one*/
					/* to the last approximation	*/
    double tol_act;			/* Actual tolerance		*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
  					/* sion operations is delayed   */
 					/* until the last moment	*/
    double new_step;      		/* Step at this iteration       */
    double dbl_eps = std::numeric_limits<double>::epsilon();
	  
    if( fabs(fc[0]) < fabs(fb[0]) )
    {                         		/* Swap data for b to be the 	*/
	a[0] = b[0];  b[0] = c[0];  c[0] = a[0];          /* best approximation		*/
	fa[0]=fb[0];  fb[0]=fc[0];  fc[0]=fa[0];
    }
    tol_act = 2*dbl_eps*fabs(b[0]) + tol[0]/2;
    new_step = (c[0]-b[0])/2;

    if( fabs(new_step) <= tol_act || fb[0] == (double)0 )
    {
      return b;				/* Acceptable approx. is found	*/
    }

    			/* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	&& fabs(fa[0]) > fabs(fb[0]) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
	register double t1,cb,t2;
	cb = c[0]-b[0];
	if( a[0]==c[0] )			/* If we have only two distinct	*/
	{				/* points linear interpolation 	*/
	  t1 = fb[0]/fa[0];			/* can only be applied		*/
	  p = cb*t1;
	  q = 1.0 - t1;
 	}
	else				/* Quadric inverse interpolation*/
	{
	  q = fa[0]/fc[0];  t1 = fb[0]/fc[0];  t2 = fb[0]/fa[0];
	  p = t2 * ( cb*q*(q-t1) - (b[0]-a[0])*(t1-1.0) );
	  q = (q-1.0) * (t1-1.0) * (t2-1.0);
	}
	if( p>(double)0 )
	{/* p was calculated with the op-*/
	  q = -q;	
	}/* posite sign; make p positive	*/
	else	
	{/* and assign possible minus to	*/
	  p = -p;
	}/* q				*/

	if( p < (0.75*cb*q-fabs(tol_act*q)/2)	/* If b+p/q falls in [b,c]*/
	    && p < fabs(prev_step*q/2) )
	{/* and isn't too large	*/
	  new_step = p/q;
	}/* it is accepted	*/
					/* If p/q is too large then the	*/
					/* bissection procedure can 	*/
					/* reduce [b,c] range to more	*/
					/* extent			*/
    }

    if( fabs(new_step) < tol_act )
    {/* Adjust the step to be not less*/
      if( new_step > (double)0 )
      {/* than tolerance		*/
	new_step = tol_act;
      }
      else
      {
	new_step = -tol_act;
      }
    }

    a[0] = b[0];  fa[0] = fb[0];			/* Save the previous approx.	*/
    b[0] += new_step;  fb = f(b, u, v, y, z);	/* Do step to a new approxim.	*/
    if( (fb[0] > 0 && fc[0] > 0) || (fb[0] < 0 && fc[0] < 0) )
    {                 			/* Adjust c for it to have a sign*/
      c[0] = a[0];  fc[0] = fa[0];                  /* opposite to that of b	*/
    }
  }

}

Rcpp::NumericVector root_function(NumericVector x, NumericVector Sampsj, NumericVector SL, NumericVector aL, NumericVector lambdaL) {
	
	int L = int(aL[0]);
	NumericVector sigsolnsj(1,0.0);
	NumericVector f(1,0.0);
	
	sigsolnsj[0] = std::exp(std::log(SL[L-1])-Sampsj[1]);	
	for(int k = 0; k < (L-1); k++){
		f[0] += (std::log(lambdaL[k]*x[0]+sigsolnsj[0])-std::log(SL[k]));
	}
	f[0] += Sampsj[0];
	
	return(f);
}



Rcpp::List sigmaSolve(NumericMatrix Samps, NumericVector SL, NumericVector aL, NumericVector aM, NumericVector lambdaL) {

	Rcpp::Function zeroin("zeroin");
	Rcpp::Function root_function("root_function");
	List result;
	int L = int(aL[0]);
	int M = int(aM[0]);
	NumericVector tol(1,0.0001);
	NumericVector sigsolnsj(1,0.0);
	NumericVector solnj(1,0.0);
	NumericVector soln(1,99.0);
	NumericVector zeroes = NumericVector(M*2, 0.0); 
        NumericMatrix solution = NumericMatrix(M, 2, zeroes.begin());
	NumericVector u(1,0.00001);
	NumericVector l(1,20000.0);
	NumericVector fu(1,0.0);
	NumericVector fl(1,0.0);
	NumericVector Sampsj(2,0.0);
	
	for(int k = 0; k < (M-1); k++){
		Sampsj[0] = Samps(k,0);Sampsj[1] = Samps(k,1);
		sigsolnsj[0] = std::exp(std::log(SL[L-1])-Sampsj[1]);	
		fl = root_function(l, Sampsj, SL, aL, lambdaL);
		fu = root_function(u, Sampsj, SL, aL, lambdaL);
		if(fl[0]*fu[0] < 0.0) 
		{
			soln = zeroin(l, u, Sampsj, SL, aL, lambdaL, root_function, tol);
		}
		solution(k,0) = soln[0]; solution(k,1) = sigsolnsj[0];
	}
	
	result = Rcpp::List::create(Rcpp::Named("solution") = solution);

	return result;
	
	
}
	
	

