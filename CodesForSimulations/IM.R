##################################
##
##
##
##
##
##
##
##
##################################

############ Dependent Packages

# matlib, lme4, merDeriv

############ Functions

######## Functions for setting up one-way random effects ANOVA (variance components model)

## make.Z
## make.Z takes a vector n_i as input and outputs a design matrix Z.  The length of n_i is the number of groups.  
## The value of each element of n_i is the number of replications within that group. 
## e.g. make.Z(c(5,5,5,5)) creats a design matrix for a balanced ANOVA with 4 groups each having 5 replicates
 
make.Z <- function(n_i){
	I <- length(n_i)
	n <- sum(n_i)
	Z <- matrix(NA,n,I)
	for(i in 1: I){
		if(i==1){
			z <- c(rep(1,n_i[1]), rep(0, sum(n_i[-1])))
		}else if(i != I){
			z <- c(rep(0, sum(n_i[1:(i-1)])), rep(1,n_i[i]), rep(0, sum(n_i[(i+1):I])))
		}else {
			z <- c(rep(0, sum(n_i[1:(I-1)])), rep(1, n_i[I]))
		}
		Z[,i]<-z
	}
	return(Z)
}


## sample.data
## sample.data takes as argument a vector n_i used for making a design matrix Z (see make.Z),
## and values for the group variance (sigma_a2), the common variance (sigma_2), and the mean (mu).
## It returns a list with the response (Y) and the design matrix (Z).

sample.data <- function(n_i, sigma_a2, sigma_2, mu){
	I <- length(n_i)
	n <- sum(n_i)
	Z <- make.Z(n_i)
	alpha <- matrix(rnorm(I, 0, sigma_a2),I,1)
	eps <- matrix(rnorm(n, 0, sigma_2),n,1)
	Y <- matrix(rep(mu,n),n,1) + Z%*%alpha + eps
	return(list(Y=Y, Z=Z))
}

sample.Y <- function(Z, sigma_a2, sigma_2, mu){
	I <- ncol(Z)
	n <- sum(Z)
	alpha <- matrix(rnorm(I, 0, sigma_a2),I,1)
	eps <- matrix(rnorm(n, 0, sigma_2),n,1)
	Y <- matrix(rep(mu,n),n,1) + Z%*%alpha + eps
	return(Y)
}

######## Function for Common computations for associations given a data set

## aov.statistics
## aov.statistics takes as argument the output of sample.data, which is response vector and design matrix
## and outputs as a list the eignevalues and multiplicities of t(P)GP, along with the sums of squares S; see
## equation (2) of the paper.  Note, these statistics are not dependent on fixed variance components

aov.statistics <- function(aov.data){
	Y <- aov.data$Y
	Z <- aov.data$Z
	n <- length(Y)

	J <- diag(n)-(1/n)*matrix(1,n,n)
	K.eig <- eigen(J)
	K <- K.eig$vectors[,1:(n-1)]%*%diag(sqrt(K.eig$values[1:(n-1)]))
	G <- t(K)%*%Z%*%t(Z)%*%K
	Ge <- eigen(G)
	lambda <- round(Ge$values,8)
	tab.l <- table(lambda)
	L <- dim(tab.l)
	lambda.L <- rep(NA,L)
	r.L <- rep(NA,L)
	for(k in 1:L){
		lambda.L[k] <- as.numeric(names(tab.l)[L-k+1])
		r.L[k] <- tab.l[L-k+1]
	}
	P.inv <- Ge$vectors
	P <- solve(P.inv)
	P <- t(P)
	S.L <- rep(NA,L)
	for(k in 1:L){
		if(k == 1){
			start <- 1
			end <- r.L[1]
		}else if(k != L){
			start <- sum(r.L[1:(k-1)])+1
			end <- sum(r.L[1:k])
		}else {
			start <- sum(r.L[1:(L-1)])+1
			end <- sum(r.L)
		}
		P.k <- P[,start:end]
		S.L[k] <- t(Y)%*%K%*%P.k%*%t(P.k)%*%t(K)%*%Y
	}
	return(list(lambda.L = lambda.L, r.L = r.L, S.L = S.L))
}



######## Functions for random set for mean of k future responses


## function: rand.sets
## inputs:  U - an M x 1 matrix iid random unif(0,1), e.g. runif(M) 
##		k - how many future responses to predict mean for
##	     aov.data - output of the function sample.data
##	     aov.statistics - output of the function aov.statistics
##         sig02 - a 2-dim vector of fixed (positive) values for the two variance components in order
##                 of sigma_a2 and sigma_2
## output: an M x 4-dim matrix, each row gives two random sets (two intervals) for the mean of k future responses from the same group.
##		The first interval corresponds to a "within-experiment" prediction for the mean of k new responses from an existing group
##		which is assumed to be the last group in the aov.data structure.  The second interval corresponds to a "out-of-experiment"
##		prediction for mean of k new responses from a new group.


rand.sets <- function(sig02, U, k, aov.data, aov.statistics, predgrid){
	Y <- aov.data$Y
	Z <- aov.data$Z
	n_i <- colSums(Z)
	n <- length(Y)	
	lambda.L <- aov.statistics$lambda.L
	r.L <- aov.statistics$r.L 
	S.L <- aov.statistics$S.L
	sigma_a02 <- sig02[1]
	sigma_02 <- sig02[2]
	
	L <- length(lambda.L)
	if(L==2){
		samps <-cbind(rchisq(10000,r.L[1]), rchisq(10000,r.L[2]))
		s2 <- S.L[2]/samps[,2]
		sigma.solns <-cbind((S.L[1]/samps[,1] - s2)/lambda.L[1], s2)
	}else {
		W0 <- matrix(NA, L, 2)
			for(l in 1:(L-1)){
				W0[l,] <- -c(lambda.L[l]/(lambda.L[l]*sigma_a02 + sigma_02), 1/(lambda.L[l]*sigma_a02 + sigma_02))
			}
		W0[L,] <- -c(0, 1/sigma_02)
		W0t <- t(W0)
		Pi_soln <- echelon(W0t, matrix(0,2,1))
		Pi.mat <- function(Pi_soln){
			L <- ncol(Pi_soln)-1
			Pi <- matrix(0, L-2,L)
			for(l in 1:(L-2)){
				v <- matrix(runif(L-2),L-2,1)
				lin.comb1 <- matrix(Pi_soln[1,3:L],1,L-2)%*%v
				lin.comb2 <- matrix(Pi_soln[2,3:L],1,L-2)%*%v
				Pi[l,]<-c(-lin.comb1, -lin.comb2, v)
			}
			return(Pi)
		}
		Pi <- Pi.mat(Pi_soln)
		A <- rbind(Pi, c(rep(1, L-1), 0), c(rep(0, L-1), 1))
		H <- Pi%*%matrix(log(S.L/(lambda.L*sigma_a02+sigma_02)),L, 1)
		A.inv <- solve(A)
		samps <- randsetsMCMC(H,A.inv,r.L,dim(H),10000)
		logdens <- samps$logdens
		samps <- cbind(samps$samples1, samps$samples2)
		#### Predict and combine steps
		#cpp.sig.solns <- sigmaSolve(samps, S.L, L, 10000, lambda.L)  # c++ version of this predict/combine step, not faster than R
		sigma.solve <- function(samps){
			M <- nrow(samps)
			sig.solns <- exp(log(S.L[L])-samps[,2])
			siga.soln <- function(arg){
				tau1<-arg[1]
				sig<-arg[2]
				g <- function(x) sum(log(lambda.L[1:(L-1)]*x+sig))+tau1 - sum(log(S.L[1:(L-1)]))
				if(g(.00001)*g(4000)>=0){
					out <- NA
				}else {
					out <- uniroot(g, c(0.000001,20000))$root
				}
				return(out)
			}
			siga.solns <- apply(matrix(cbind(samps[,1], sig.solns), M,2),1,siga.soln)
			sigs.solns <- matrix(cbind(siga.solns, sig.solns),M,2)
			sigs.solns <- sigs.solns[complete.cases(sigs.solns), ]
			return(sigs.solns)
		}
		sigma.solns <- sigma.solve(samps)
	}
	rand.sets.pred <- randsetspred(sigma.solns, dim(sigma.solns), n, n_i, length(n_i), k, U, mean(Y))$randsetpred
	Omega1 <- sum(log(S.L[1:(L-1)])) - sum(log(lambda.L[1:(L-1)]*sigma_a02 + sigma_02))
	Omega2 <- log(S.L[L]) - log(sigma_02)
	matprod <- A.inv%*%matrix(c(H,Omega1, Omega2), L,1)
	dens <- sum(0.5*r.L*matprod-0.5*exp(matprod))
	densplaus <- randsetspreddens(logdens, length(logdens), n, n_i, length(n_i), k, mean(Y), predgrid, length(predgrid), c(sigma_a02, sigma_02), dens)
	return(rand.sets.pred)
}

# same function but entirely in R code (no c++), not used
rand.sets.R <- function(sig02, U, k, aov.data, aov.statistics){
	Y <- aov.data$Y
	Z <- aov.data$Z
	n_i <- colSums(Z)
	n <- length(Y)	
	lambda.L <- aov.statistics$lambda.L
	r.L <- aov.statistics$r.L 
	S.L <- aov.statistics$S.L
	sigma_a02 <- sig02[1]
	sigma_02 <- sig02[2]
	
	L <- length(lambda.L)
	W0 <- matrix(NA, L, 2)
	for(l in 1:(L-1)){
		W0[l,] <- -c(lambda.L[l]/(lambda.L[l]*sigma_a02 + sigma_02), 1/(lambda.L[l]*sigma_a02 + sigma_02))
	}
	W0[L,] <- -c(0, 1/sigma_02)
	W0t <- t(W0)
	Pi_soln <- echelon(W0t, matrix(0,2,1))
	Pi.mat <- function(Pi_soln){
		L <- ncol(Pi_soln)-1
		Pi <- matrix(0, L-2,L)
		for(l in 1:(L-2)){
			v <- matrix(runif(L-2),L-2,1)
			lin.comb1 <- matrix(Pi_soln[1,3:L],1,L-2)%*%v
			lin.comb2 <- matrix(Pi_soln[2,3:L],1,L-2)%*%v
			Pi[l,]<-c(-lin.comb1, -lin.comb2, v)
		}
		return(Pi)
	}
	Pi <- Pi.mat(Pi_soln)
	A <- rbind(Pi, c(rep(1, L-1), 0), c(rep(0, L-1), 1))
	H <- Pi%*%matrix(log(S.L/(lambda.L*sigma_a02+sigma_02)),L, 1)
	log.den_logchisq.df <- function(v, df){
		return(0.5*df*v-0.5*exp(v))
	}
	joint.log.dens <- function(u){
		U<-rbind(H, u[1], u[2])
		A.inv <- solve(A)
		arg <- A.inv%*%U
		log.joint <- 0
		for(j in 1:length(U)) log.joint <- log.joint + log.den_logchisq.df(arg[j], r.L[j])
		return(log.joint)
	}
	mh_j <- function(j, dprop, rprop, u) {  # Metropolis Hastings steps for each component of theta (MH within Gibbs)
		u_prop <- u
		uj <- rprop(u[j],j)
		u_prop[j] <- uj
		r <- joint.log.dens(u_prop) - joint.log.dens(u) + log(dprop(uj , u[j], j)) - log(dprop(u[j] , uj, j))
   		R <- min(exp(r), 1)
    		if(runif(1) <= R) {
      		u <- u_prop
    		} 
     		return(u)
	}
	dprop <- function(x,theta,i){   # MH proposal distribution density
		if(i==1){
			s <- 10
		}else if(i == 2){
			s <- 0.5
		}else {
			s <- 0.5
		}
    		return(dnorm(x, mean = theta,sd = s))
	}
	rprop <- function(theta,i){     # MH proposal distribution sampling
		if(i==1){
			s <- 10
		}else if(i == 2){
			s <- 0.5
		}else {
			s <- 0.5
		}
    		return(rnorm(1,mean = theta, sd = s))
	}
	M<- 10000
	samps <- matrix(NA, M, 2)
	acc <- c(0,0)
	u <- c(1,1)
	for( s in 1:(M+100) ){
		for(j in 1:2){
			uold <- u
			u <- mh_j(j, dprop, rprop, u) 
			if(j ==1){
			  acc[1] <- acc[1] + ifelse(u[1]==uold[1],0,1)
			} else {
			  acc[2] <- acc[2] + ifelse(u[2]==uold[2],0,1)
			}
		}
		if(s>=101) samps[s-100,]<-u
	}
	#### Predict and combine steps
	sigma.solve <- function(samps){
		M <- nrow(samps)
		sig.solns <- exp(log(S.L[L])-samps[,2])
		siga.soln <- function(arg){
			tau1<-arg[1]
			sig<-arg[2]
			g <- function(x) sum(log(lambda.L[1:(L-1)]*x+sig))+tau1 - sum(log(S.L[1:(L-1)]))
			if(g(.00001)*g(4000)>=0){
				out <- NA
			}else {
				out <- uniroot(g, c(0.000001,20000))$root
			}
			return(out)
		}
		siga.solns <- apply(matrix(cbind(samps[,1], sig.solns), M,2),1,siga.soln)
		sigs.solns <- matrix(cbind(siga.solns, sig.solns),M,2)
		sigs.solns <- sigs.solns[complete.cases(sigs.solns), ]
		return(sigs.solns)
	}
	sigma.solns <- sigma.solve(samps)
	rand.norms <- rnorm(nrow(sigma.solns))
	Q.samps <- function(Zs, sigmas){
		triplet.solve <- function(vec){
			z<-vec[1]
			siga<-vec[2]
			sig<-vec[3]
			return(z*sqrt(siga*(1+(1/(n^2))*sum(n_i^2))+sig*(1/n+1/k)))
		}
		all.triples <- matrix(cbind(Zs, sigmas),nrow(sigma.solns),3)
		return(sort(apply(all.triples, 1, triplet.solve)))
	}
	Qs <- Q.samps(rand.norms, sigma.solns)
	rand.set.pred.U <- function(u){
		U.l <- 0.5-abs(u-0.5)
		U.h <- 0.5+abs(u-0.5)
		rand.set.Q <- quantile(Qs, probs = c(U.l, U.h))
		rand.set.pred <- mean(Y) + rand.set.Q
		return(rand.set.pred)
	}
	rand.sets.pred <- apply(U,1,rand.set.pred.U)
	return(t(rand.sets.pred))
}






######## Functions for computing Plausibility 

## function: pl.0
## inputs: vals - an M x 1 matrix of values (a grid) for which to compute plausibilities 
##	     rand.sets - output of the function rand.sets
## output: a list containing 'plauses' an M x 1 vector of proportions, plausibilites for each value in vals
##		and 'vals', the input

pl.0 <- function(vals, randsets){
	M <- nrow(randsets)
	pl.0.val <- function(val, randsets){
		plaus <- sum((val > randsets[,1])*(val < randsets[,2]))/M
		return(plaus)
	}
	pl.0.val.T <- function(val, randsets){
		randsets1<-randsets[,1];randsets2<-randsets[,2]
		none.na<- (1-apply(matrix(cbind(is.na(randsets1), is.na(randsets2)), M, 2),1,any))
		M <-	sum(none.na)
		plaus <- sum(((val > randsets1)*(val < randsets2))[none.na==1])/M
		return(plaus)
	}
	return(list(plauses.w = apply(vals,1,pl.0.val, randsets = randsets[,1:2]), plauses.n = apply(vals,1,pl.0.val, randsets = randsets[,3:4]),
		plauses.T = apply(vals,1,pl.0.val.T, randsets = randsets[,5:6]), vals = vals))
}

plot.plaus <- function(pl.0){
	return(plot(pl.0$vals, pl.0$plauses.n, type = 'l'))
}



## function: apply.plaus
## inputs: sig02 - a N x 2-dim matrix of pairs of fixed (positive) values for the two variance components in order
##                 of sigma_a2 and sigma_2
##	     U - an M x 1 matrix iid random unif(0,1), e.g. runif(M)
##	     k - how many future responses to predict mean for
##	     aov.data - output of the function aov.data
##	     aov.statistics - output of the function aov.statistics
##	     vals - an L x 1 matrix of values (a grid) for which to compute plausibilities
##	     within - integer 1 or 0 indicating if prediction is desired for "within-experiment" group or a new group
##			  if 1 it is assumed the desired group is the last group (Ith group) according to the structure in sample.data
##         
## output: an N x L-dim matrix, each row is a local plausibility contour corresponding to a pair of 
## 		sigma_a2 and sigma_2.  Taking the columnwise maxima of output gives a (fused) plausibility
##		contour for the mean of the next k responses.


apply.plaus <- function(sig02, U, k, ex.data, ex.statistics, vals){
	N <- nrow(sig02)
	L <- length(vals)
	all.plaus.w <- matrix(NA, N, L)
	all.plaus.n <- matrix(NA, N, L)
	all.plaus.T <- matrix(NA, N, L)
	for(j in 1:N){
		pl.j <- pl.0(vals, rand.sets(sig02[j,], U, k, ex.data, ex.statistics))
		all.plaus.w[j,] <- pl.j$plauses.w
		all.plaus.n[j,] <- pl.j$plauses.n
		all.plaus.T[j,] <- pl.j$plauses.T
	}
	fused.plaus.w <- apply(all.plaus.w, 2, max)
	fused.plaus.n <- apply(all.plaus.n, 2, max)
	fused.plaus.T <- apply(all.plaus.T, 2, max)
	return(list(all.plaus.w = all.plaus.w, all.plaus.n = all.plaus.n, all.plaus.T = all.plaus.T, fused.plaus.w = fused.plaus.w, fused.plaus.n = fused.plaus.n, fused.plaus.T = fused.plaus.T))
}




############
############	TESTING FUNCTIONS
############
############



# testing cpp codes

#par(mfrow = c(1,2))
#ex.data <- sample.data(n_i = c(4,4,4,6,32,17,8), sigma_a2 = 2, sigma_2 = 1, mu = 4)
#ex.statistics <- aov.statistics(ex.data)
#t1 <- proc.time()
#ex.rand.sets <- rand.sets(sig02 = c(0.5,0.5), matrix(runif(10000),10000,1), 10, ex.data, ex.statistics)
#ex.grid <- matrix(seq(from = -20, to = 25, length.out = 1000),1000,1)
#ex.plauses <- pl.0(ex.grid, ex.rand.sets)
#proc.time()-t1
#plot.plaus(ex.plauses)

#par(mfrow = c(1,2))
#plot(ex.grid, ex.plauses$plauses.w, type = 'l')
#plot(ex.grid, ex.plauses$plauses.n, type = 'l')


#plaus.grid <- seq(from = min(ex.data$Y)-6, to = max(ex.data$Y)+6, length.out = 500)
#sa <- seq(.01, 4, length.out = 10)
#s0 <- seq(.01, 4, length.out = 10)
#fusing.grid <- as.matrix(expand.grid(x = sa, y = s0))
#t1 <- proc.time()
#fused <- apply.plaus(fusing.grid, matrix(runif(10000),10000,1), 1, ex.data, ex.statistics, matrix(plaus.grid, length(plaus.grid),1))
#proc.time()-t1


#plot(plaus.grid, fused$fused.plaus.w, type = 'l')
#for(j in 1:100){
#lines(plaus.grid, fused$all.plaus.w[j,], col = 'red')
#}
#lines(plaus.grid, fused$fused.plaus.w, lwd = 2, col = 'blue')

#plot(plaus.grid, fused$fused.plaus.n, type = 'l')
#for(j in 1:100){
#lines(plaus.grid, fused$all.plaus.n[j,], col = 'red')
#}
#lines(plaus.grid, fused$fused.plaus.n, lwd = 2, col = 'blue')







