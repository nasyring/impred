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

inline double GibbsMCMC(RVector<double> nn, RMatrix<double> data, RVector<double> nnp, RMatrix<double> priordata, RVector<double> priorweight, RMatrix<double> thetaboot,
	RVector<double> bootmean0, RVector<double> bootmean1, RMatrix<double> databoot, RVector<double> scheduleLen, RMatrix<double> priorSched,
	RVector<double> alpha, RVector<double> M_samp, RVector<double> w, std::size_t i) {
   	

	int M = int(M_samp[0]);
	int n = int(nn[0]);
	int np = int(nnp[0]);
	int sN = int(scheduleLen[0]);
   	NumericVector prop0(1,0.0);
   	NumericVector prop1(1,0.0);
	double cov_ind = 0.0;
	NumericVector YIboot(1,0.0);
   	NumericVector datamin(1,0.0);
	NumericVector datamax(1,0.0);
   	NumericVector theta0old(1,0.0);
	NumericVector theta0new(1,0.0);
	NumericVector theta1old(1,0.0);
	NumericVector theta1new(1,0.0);
	NumericVector loglikdiff(1,0.0);
	NumericVector r(1,0.0);
	NumericVector uu(1,0.0);
	NumericVector vv(1,0.0);
	NumericVector postsamples0(M,0.0);
	NumericVector postsamples1(M,0.0);
	NumericVector logpost(M,0.0);
	NumericVector l0(1,0.0);
	NumericVector l1(1,0.0);
	NumericVector u0(1,0.0);
	NumericVector u1(1,0.0);
	NumericVector n0(1,0.0);
	NumericVector n1(1,0.0);
	NumericVector n2(1,0.0);
	NumericVector F0_c0old(1,0.0);
	NumericVector F0_c0new(1,0.0);
	NumericVector F1_c0old(1,0.0);
	NumericVector F1_c0new(1,0.0);
	NumericVector F1_c1old(1,0.0);
	NumericVector F1_c1new(1,0.0);
	NumericVector F2_c1old(1,0.0);
	NumericVector F2_c1new(1,0.0);
	NumericVector n0p(1,0.0);
	NumericVector n1p(1,0.0);
	NumericVector n2p(1,0.0);
	NumericVector F0_c0oldp(1,0.0);
	NumericVector F0_c0newp(1,0.0);
	NumericVector F1_c0oldp(1,0.0);
	NumericVector F1_c0newp(1,0.0);
	NumericVector F1_c1oldp(1,0.0);
	NumericVector F1_c1newp(1,0.0);
	NumericVector F2_c1oldp(1,0.0);
	NumericVector F2_c1newp(1,0.0);
	theta0old(0) = thetaboot(i,0);
	theta1old(0) = thetaboot(i,1);
	NumericVector YI(M,0.0);
	NumericVector YIl(1,0.0);
	NumericVector YIu(1,0.0);
	datamin(0) = databoot(0,2*i+1);
	datamax(0) = databoot(0,2*i+1);
	
	for(int j=0; j<sN; j++){
		if(w[0]<=priorSched(j,0)){
			prop0(0) = priorSched(j,1);	
		}
		if(w[0]<=priorSched(j,2)){
			prop1(0) = priorSched(j,3);	
		}		
	}

	
	for(int k=0; k<n; k++){
		if(databoot(k,2*i)==1){
			n0(0) = n0(0) + 1;
			if(databoot(k,2*i+1)<=theta0old(0)){
				F0_c0old(0) = 	F0_c0old(0)+1;
			}	
		}
		else if(databoot(k,2*i)==2){
			n1(0) = n1(0) + 1;
			if(databoot(k,2*i+1)<=theta0old(0)){
				F1_c0old(0) = 	F1_c0old(0)+1;
			}
			if(databoot(k,2*i+1)<=theta1old(0)){
				F1_c1old(0) = 	F1_c1old(0)+1;
			}	
		}
		else {
			n2(0) = n2(0) + 1;
			if(databoot(k,2*i+1)<=theta1old(0)){
				F2_c1old(0) = 	F2_c1old(0)+1;
			}
		}
		if(datamin(0) > databoot(k,2*i+1)){
			datamin(0) = databoot(k,2*i+1);	
		}
		if(datamax(0) < databoot(k,2*i+1)){
			datamax(0) = databoot(k,2*i+1);	
		}
		YIboot(0) = YIboot(0)+thetaboot(k,2);
	}
	YIboot(0) = YIboot(0)/n;
	F0_c0old(0) = 	F0_c0old(0)/n0(0);
	F1_c0old(0) = 	F1_c0old(0)/n1(0);
	F1_c1old(0) = 	F1_c1old(0)/n1(0);
	F2_c1old(0) = 	F2_c1old(0)/n2(0);
	
	for(int k=0; k<np; k++){
		if(priordata(k,0)==1){
			n0p(0) = n0p(0) + 1;
			if(priordata(k,1)<=theta0old(0)){
				F0_c0oldp(0) = 	F0_c0oldp(0)+1;
			}	
		}
		else if(priordata(k,0)==2){
			n1p(0) = n1p(0) + 1;
			if(priordata(k,1)<=theta0old(0)){
				F1_c0oldp(0) = 	F1_c0oldp(0)+1;
			}
			if(priordata(k,1)<=theta1old(0)){
				F1_c1oldp(0) = 	F1_c1oldp(0)+1;
			}	
		}
		else {
			n2p(0) = n2p(0) + 1;
			if(priordata(k,1)<=theta1old(0)){
				F2_c1oldp(0) = 	F2_c1oldp(0)+1;
			}
		}
	}
	F0_c0oldp(0) = 	F0_c0oldp(0)/n0p(0);
	F1_c0oldp(0) = 	F1_c0oldp(0)/n1p(0);
	F1_c1oldp(0) = 	F1_c1oldp(0)/n1p(0);
	F2_c1oldp(0) = 	F2_c1oldp(0)/n2p(0);

	
	for(int j=0; j<M; j++) {
		theta0new(0) = R::rnorm(theta0old(0), prop0(0));
		theta1new(0) = R::rnorm(theta1old(0), prop1(0));
	
		if( (theta0new(0)>theta1new(0)) ||  (theta0new(0)<datamin(0)) || (theta1new(0)>datamax(0))   ){
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}else {
		
		loglikdiff(0) = 0.0;
		for(int k=0; k<n; k++){
			if(databoot(k,2*i)==1){
				if(databoot(k,2*i+1)<=theta0new(0)){
					F0_c0new(0) = 	F0_c0new(0)+1.0;
				}	
			}
			else if(databoot(k,2*i)==2){
				if(databoot(k,2*i+1)<=theta0new(0)){
					F1_c0new(0) = 	F1_c0new(0)+1.0;
				}
				if(databoot(k,2*i+1)<=theta1new(0)){
					F1_c1new(0) = 	F1_c1new(0)+1.0;
				}	
			}
			else {
				if(databoot(k,2*i+1)<=theta1new(0)){
					F2_c1new(0) = 	F2_c1new(0)+1.0;
				}
			}
		}
		F0_c0new(0) = 	F0_c0new(0)/n0(0);
		F1_c0new(0) = 	F1_c0new(0)/n1(0);
		F1_c1new(0) = 	F1_c1new(0)/n1(0);
		F2_c1new(0) = 	F2_c1new(0)/n2(0);
		
		for(int k=0; k<np; k++){
			if(priordata(k,0)==1){
				if(priordata(k,1)<=theta0new(0)){
					F0_c0newp(0) = 	F0_c0newp(0)+1.0;
				}	
			}
			else if(priordata(k,0)==2){
				if(priordata(k,1)<=theta0new(0)){
					F1_c0newp(0) = 	F1_c0newp(0)+1.0;
				}
				if(priordata(k,1)<=theta1new(0)){
					F1_c1newp(0) = 	F1_c1newp(0)+1.0;
				}	
			}
			else {
				if(priordata(k,1)<=theta1new(0)){
					F2_c1newp(0) = 	F2_c1newp(0)+1.0;
				}
			}
		}
		F0_c0newp(0) = 	F0_c0newp(0)/n0p(0);
		F1_c0newp(0) = 	F1_c0newp(0)/n1p(0);
		F1_c1newp(0) = 	F1_c1newp(0)/n1p(0);
		F2_c1newp(0) = 	F2_c1newp(0)/n2p(0);
	
		loglikdiff(0) = w[0]*n*((F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))-(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0)))  +   w[0]*priorweight[0]*np*((F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0))-(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0)));
		loglikdiff(0) = fmin(std::exp(loglikdiff(0)), 1.0);
		uu[0] = R::runif(0.0,1.0);
		if(uu(0) <= loglikdiff(0)) {
			postsamples0(j) = theta0new(0);
			postsamples1(j) = theta1new(0);
			logpost(j) = w[0]*n*(F0_c0new(0)+F1_c1new(0)-F1_c0new(0)-F2_c1new(0))+w[0]*priorweight[0]*np*(F0_c0newp(0)+F1_c1newp(0)-F1_c0newp(0)-F2_c1newp(0));
			theta0old(0) = theta0new(0);
			theta1old(0) = theta1new(0);
			F0_c0old(0)=F0_c0new(0);
			F1_c0old(0)=F1_c0new(0);
			F1_c1old(0)=F1_c1new(0);
			F2_c1old(0)=F2_c1new(0);
			F0_c0oldp(0)=F0_c0newp(0);
			F1_c0oldp(0)=F1_c0newp(0);
			F1_c1oldp(0)=F1_c1newp(0);
			F2_c1oldp(0)=F2_c1newp(0);
			YI(j) = F0_c0new(0) + F1_c1new(0) - F1_c0new(0) - F2_c1new(0);
		}
		else {
			postsamples0(j) = theta0old(0);
			postsamples1(j) = theta1old(0);
			logpost(j) = w[0]*n*(F0_c0old(0)+F1_c1old(0)-F1_c0old(0)-F2_c1old(0))+w[0]*priorweight[0]*np*(F0_c0oldp(0)+F1_c1oldp(0)-F1_c0oldp(0)-F2_c1oldp(0));
			YI(j) = F0_c0old(0) + F1_c1old(0) - F1_c0old(0) - F2_c1old(0);
		}
		}
	}	

	std::sort(postsamples0.begin(), postsamples0.end());
	std::sort(postsamples1.begin(), postsamples1.end());
	std::sort(YI.begin(), YI.end());
	l0[0] = postsamples0(M*.025-1);
	u0[0] = postsamples0(M*.975-1);
	l1[0] = postsamples1(M*.025-1);
	u1[0] = postsamples1(M*.975-1);
	YIl[0] = YI[M*.025-1];
	YIu[0] = YI[M*.975-1];
	/*if ( (YIl[0] < YIboot[0]) && (YIu[0] > YIboot[0]) ){
			cov_ind = 1.0;
	} else {cov_ind = 0.0;}*/
	if ( (l0[0] < bootmean0[0]) && (u0[0] > bootmean0[0]) ){
			cov_ind = 1.0;
	} else {cov_ind = 0.0;}
	return cov_ind;

	
}






struct GPCYI_yi_mcmc_parallel : public Worker {

	const RVector<double> nn;
	const RMatrix<double> data;
	const RVector<double> nnp;
	const RMatrix<double> priordata;
	const RVector<double> priorweight;
	const RMatrix<double> thetaboot;
	const RVector<double> bootmean0;
	const RVector<double> bootmean1;
	const RMatrix<double> databoot;
	const RVector<double> scheduleLen;
	const RMatrix<double> priorSched;
	const RVector<double> alpha;
	const RVector<double> M_samp;
	const RVector<double> B_resamp;
	const RVector<double> w;
	RVector<double> cover;

   // initialize with source and destination
   GPCYI_yi_mcmc_parallel(const NumericVector nn, const NumericMatrix data, const NumericVector nnp, const NumericMatrix priordata, const NumericVector priorweight, const NumericMatrix thetaboot,
	const NumericVector bootmean0, const NumericVector bootmean1, const NumericMatrix databoot, const NumericVector scheduleLen, const NumericMatrix priorSched,
	const NumericVector alpha, const NumericVector M_samp, const NumericVector B_resamp,
	const NumericVector w, NumericVector cover) 
			: nn(nn), data(data), nnp(nnp), priordata(priordata), priorweight(priorweight),  thetaboot(thetaboot), bootmean0(bootmean0), bootmean1(bootmean1), databoot(databoot), scheduleLen(scheduleLen), priorSched(priorSched), alpha(alpha), M_samp(M_samp), B_resamp(B_resamp), w(w), cover(cover) {}   

   // operator
void operator()(std::size_t begin, std::size_t end) {
		for (std::size_t i = begin; i < end; i++) {
			cover[i] = GibbsMCMC(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, databoot, scheduleLen, priorSched, alpha, M_samp, w, i);	
		}
	}
};

// [[Rcpp::export]]
NumericVector rcpp_parallel_yi(NumericVector nn, NumericMatrix data, NumericVector nnp, NumericMatrix priordata, NumericVector priorweight, NumericMatrix thetaboot, NumericVector bootmean0,
	NumericVector bootmean1, NumericMatrix databoot, NumericVector scheduleLen, NumericMatrix priorSched, NumericVector alpha, NumericVector M_samp, NumericVector B_resamp,
	NumericVector w) {
	
   int B = int(B_resamp[0]);
   // allocate the matrix we will return
   NumericVector cover(B,2.0); 

   // create the worker
   GPCYI_yi_mcmc_parallel gpcWorker(nn, data, nnp, priordata, priorweight, thetaboot, bootmean0, bootmean1, databoot, scheduleLen, priorSched, alpha, M_samp, B_resamp, w, cover);
     
   // call it with parallelFor
   
   parallelFor(0, B, gpcWorker);

   return cover;
}

// [[Rcpp::export]]
Rcpp::List GPCYI_yi_parallel(SEXP & nn, SEXP & data, SEXP & nnp, SEXP & priordata, SEXP & priorweight, SEXP & theta_boot, SEXP & data_boot, SEXP & scheduleLen, SEXP & priorSched, SEXP & alpha, SEXP & M_samp, SEXP & B_resamp) {

RNGScope scp;
Rcpp::Function _GPCYI_rcpp_parallel_yi("rcpp_parallel_yi");
List result;
List finalsample;
double eps 			= 0.01; 
NumericVector nn_ = Rcpp::as<NumericVector>(nn);
NumericMatrix data_ = Rcpp::as<NumericMatrix>(data);
NumericVector nnp_ = Rcpp::as<NumericVector>(nnp);
NumericMatrix priordata_ = Rcpp::as<NumericMatrix>(priordata);
NumericVector priorweight_ = Rcpp::as<NumericVector>(priorweight);
NumericMatrix thetaboot_ = Rcpp::as<NumericMatrix>(theta_boot);
NumericVector bootmean0(1,0.0);
NumericVector bootmean1(1,0.0);
NumericMatrix databoot_ = Rcpp::as<NumericMatrix>(data_boot);
NumericVector scheduleLen_ = Rcpp::as<NumericVector>(scheduleLen);
NumericMatrix priorSched_ = Rcpp::as<NumericMatrix>(priorSched);
NumericVector alpha_ = Rcpp::as<NumericVector>(alpha);
NumericVector M_samp_ = Rcpp::as<NumericVector>(M_samp);
NumericVector B_resamp_ = Rcpp::as<NumericVector>(B_resamp);
NumericVector w(1,0.5);
double diff;
bool go 			= TRUE;
int t				=1; 
double sumcover;
int B = int(B_resamp_[0]);
NumericVector cover;
	
for (int i=0; i<B; i++) {
	bootmean0[0] = bootmean0[0] + thetaboot_(i,0);
	bootmean1[0] = bootmean1[0] + thetaboot_(i,1);
}
bootmean0 = bootmean0/B;
bootmean1 = bootmean1/B;

while(go){	
cover = _GPCYI_rcpp_parallel_yi(nn_, data_, nnp_, priordata_, priorweight_, thetaboot_, bootmean0, bootmean1, databoot_, alpha_, M_samp_, B_resamp_, w);
sumcover = 0.0;
for(int s = 0; s<B; s++){sumcover = sumcover + cover(s);}
diff = (sumcover/B) - (1.0-alpha_[0]);
if(((abs(diff)<= eps)&&(diff>=0)) || t>16) {
   go = FALSE;
} else {
   t = t+1;
   w[0] = fmax(w[0] + (pow(1+t,-0.51)*diff),0.1);
} 
}

// Final sample

NumericVector M_final; M_final[0] = 2*M_samp_[0];
finalsample = GibbsMCMC2(nn_, data_, nnp_, priordata_, priorweight_, thetaboot_, bootmean0, bootmean1, scheduleLen_, priorSched_, alpha_, M_final, w);
	
result = Rcpp::List::create(Rcpp::Named("w") = w,Rcpp::Named("t") = t,Rcpp::Named("diff") = diff, Rcpp::Named("list_cis") = finalsample);
	
return result;
}




// [[Rcpp::export]]

inline bool compare(std::array<double, 6> a, std::array<double, 6> b){
    return (a[5] < b[5]);
}




