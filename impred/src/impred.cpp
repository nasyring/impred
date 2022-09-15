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

Rcpp::List IMTS_mh_sampler(NumericVector lU0, NumericVector V0, NumericVector H0, NumericMatrix Minv, NumericVector rL, NumericVector thetau0s, NumericVector propsd1, NumericVector propsd2){

	List result;
	int L = H0.length() + 2;
	int T = thetau0s.length();
	NumericVector rL0(L+1, 0.0); rL0[0] = 1.0; rL0[L] = rL[L-1];
	NumericVector prop1(1,0.0);
	for(int j = 0; j < (L-1); j++){
		rL0[j+1] = rL[j];
		prop1[0] = prop1[0] + std::log(rL[j]);	
	}
	arma::mat M = as<arma::mat>(Minv); M = M.t();
	NumericVector allLU(L,0.0);
	allLU[0] = lU0[0]; allLU[1] = V0[0];
	for(int i = 0; i < (L-2); i++){
		allLU[i+2] = H0[i];	
	}
	arma::vec lUn = as<arma::vec>(allLU);
	arma::vec lU = M*lUn;
	NumericVector lf(1,0.0);
	NumericVector lf1(1,0.0);NumericVector lf2(1,0.0);NumericVector lf3(1,0.0);
	for(int i = 0; i < L; i++){
		lf1[0] = lf1[0] + 0.5*lU[i]*rL0[i];
		lf2[0] = lf2[0] + 0.5*rL0[i];
		lf3[0] = lf3[0] + rL0[i]*std::exp(lU[i])/rL0[L];
	}
	lf2[0] = lf2[0] + 0.5*rL0[L];
	lf[0] = lf1[0] - lf2[0]*(0.5+0.5*lf3[0]);
	
	NumericVector lfold(1, lf[0]); NumericVector lfnew(1, 0.0);
	NumericVector uold(1, lU0[0]); NumericVector vold(1, V0[0]); 	
	NumericVector unew(1, 0.0); NumericVector vnew(1, 0.0);
	NumericVector unif1(1, 0.0); NumericVector unif2(1, 0.0);
	NumericVector logdens(5000, 0.0);  
	
	for(int m = 0; m < 5000; m++){
		unew[0] = R::rnorm(uold[0],propsd1[0]);	
		vnew[0] = R::rnorm(vold[0],propsd2[0]);

		allLU[0] = unew[0]; allLU[1] = vnew[0];
		lUn = as<arma::vec>(allLU);
		lU = M*lUn;
		lf1[0]=0.0; lf2[0]=0.0; lf3[0]=0.0;
		for(int i = 0; i < L; i++){
			lf1[0] = lf1[0] + 0.5*lU[i]*rL0[i];
			lf2[0] = lf2[0] + 0.5*rL0[i];
			lf3[0] = lf3[0] + rL0[i]*std::exp(lU[i])/rL0[L];
		}
		lf2[0] = lf2[0] + 0.5*rL0[L];
		lfnew[0] = lf1[0] - lf2[0]*(0.5+0.5*lf3[0]);		
		
		unif1[0] = std::exp(lfnew[0] - lfold[0]);
		unif2[0] = R::runif(0.0,1.0);
		if(unif1[0] >= unif2[0]){
			logdens[m] = lfnew[0];
			lfold[0] = lfnew[0];
			uold[0] = unew[0];
			vold[0] = vnew[0];			
		}else {
			logdens[m] = lfold[0];
		}
		
	}
	
	
	NumericVector lfseq(T,0.0);
	allLU[1] = V0[0];
	for(int t = 0; t < T; t++){
		allLU[0] = thetau0s[t]; 
		lUn = as<arma::vec>(allLU);
		lU = M*lUn;
		lf1[0]=0.0; lf2[0]=0.0; lf3[0]=0.0;
		for(int i = 0; i < L; i++){
			lf1[0] = lf1[0] + 0.5*lU[i]*rL0[i];
			lf2[0] = lf2[0] + 0.5*rL0[i];
			lf3[0] = lf3[0] + rL0[i]*std::exp(lU[i])/rL0[L];
		}
		lf2[0] = lf2[0] + 0.5*rL0[L];
		lfseq[t] = lf1[0] - lf2[0]*(0.5+0.5*lf3[0]);
	}
	
	
	result = Rcpp::List::create(Rcpp::Named("logdenssamps") = logdens, Rcpp::Named("logdens0") = lf, Rcpp::Named("logdensthetas") = lfseq);
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


