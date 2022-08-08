##################################
##
##	IM.R
##	Helper functions for IM computations and simulations 
##	
##
##################################

############ Dependent Packages

# matlib, lme4, merDeriv
library(devtools)
install_github("nasyring/impred", subdir = "impred")
library(impred)
library(matlib)
library(lme4)
library(merDeriv)


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
## Appendix A of the paper.  Note, these statistics are not dependent on fixed variance components

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


## lmer.statistics
## lmer.statistics - analogous to aov.statistics for general two-stage model
lmer.statistics <- function(data){
	Y <- data$Y
	X <- data$X
	Z <- data$Z
	
	n <- length(Y)
	G <- Z%*%t(Z)
	XX <- t(X)%*%X
	XXinv <- solve(XX)
	B <- XXinv%*%t(X)
	By <- B%*%Y
	In <- diag(n)
	K1 <- In - X%*%B
	Keig <- eigen(K1)
	Ke <- diag(sqrt(round(Keig$values,8)))
	Kv <- Keig$vectors
	K <- Kv%*%Ke
	H <- t(K)%*%G%*%K
	C1 <- B%*%t(B)
	C2 <- B%*%G%*%t(B)
	He <- eigen(H)
	lambda <- round(He$values,8)
	tab.l <- table(lambda)
	L <- dim(tab.l)
	lambda.L <- rep(NA,L)
	r.L <- rep(NA,L)
	for(k in 1:L){
		lambda.L[k] <- as.numeric(names(tab.l)[L-k+1])
		r.L[k] <- tab.l[L-k+1]
	}
	P.inv <- He$vectors
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

	return(list(lambda.L = lambda.L, r.L = r.L, S.L = S.L, C1=C1, C2=C2, By=By))

}




## plaus_aov_unbalanced_full
## Calls rcpp function "plaus_unbalanced_aov_full" to compute IM for random intercept model
## data is a data set in the format of sample.data
## theta is a sequence of theta prediction values at which to compute plausibility pointwise
## ratio.grid is a sequence in (0,1) of rho values 
plaus_aov_unbalanced_full <- function(data, theta, ratio.grid){
	Y <- data$Y
	Z <- data$Z
	n <- length(Y)
	n_i <- colSums(Z)
	
	Ybar <- mean(Y)
	stat <- aov.statistics(data)
	S <- stat$S.L
	lambda <- stat$lambda.L
	r <- stat$r.L
	L <- length(S)
	T <- length(theta)
	plauses.t <- rep(0, T)
	plauses.n <- rep(0, T)
	plauses.e <- rep(0, T)

	M <- length(ratio.grid)
	plauses.t.all <- matrix(0, T, M)

	
	for(m in 1:M){
		ratio <- (1-ratio.grid[m])/ratio.grid[m]
		plauses <- plaus_unbalanced_aov_full(theta, Ybar, S, lambda, r, n, n_i, ratio) 
		plauses.t.all[,m] <- plauses$plauses.theta
		plauses.t <- apply(cbind(plauses$plauses.theta, plauses.t),1,max)	
		plauses.n <- apply(cbind(plauses$plauses.new, plauses.n),1,max)	
		plauses.e <- apply(cbind(plauses$plauses.exs, plauses.e),1,max)	
	}
	
	return(list(plauses.theta = plauses.t, plauses.new = plauses.n, plauses.exs = plauses.e, plauses.all = plauses.t.all))

}

## plaus_two_stage_lme_full
## Calls rcpp function "plaus_two_stage_full" to compute IM for two stage model
## data is a data set in the format of sample.data
## theta is a sequence of theta prediction values at which to compute plausibility pointwise
## x is a covariate vector, i.e. theta has expectation t(x)\beta
## ratio.grid is a sequence in (0,1) of rho values 
plaus_two_stage_lme_full <- function(data, theta, x, ratio.grid){
	Y <- data$Y
	Z <- data$Z
	X <- data$X
	n <- length(Y)
	stat <- lmer.statistics(data)
	x <- matrix(x, length(x),1)
	xBy <- t(x)%*%stat$By
	S <- stat$S.L
	L <- length(S)
	lambda <- stat$lambda.L
	r <- stat$r.L
	C1 <- stat$C1
	C2 <- stat$C2	

	M <- length(ratio.grid)
	plauses.all.t <- matrix(0, length(theta), M)
	plauses.all.n <- matrix(0, length(theta), M)
	for(m in 1:M){
		ratio <- ratio.grid[m]
		csigma <- as.numeric(t(x)%*%(C1 * 1.0 + C2 * ratio)%*%x + ratio)
		plauses.im <- plaus_two_stage_full(theta, xBy, S, lambda, r, csigma, ratio)
		plauses.all.t[,m] <- plauses.im$plauses.theta
		plauses.all.n[,m] <- plauses.im$plauses.new
	}
	plauses.theta <- apply(plauses.all.t, 1, max)
	plauses.new <- apply(plauses.all.n, 1, max)

	return(list(plauses.theta = plauses.theta, plauses.new = plauses.new, plauses.all.theta = plauses.all.t, plauses.all.new = plauses.all.n))

}


## t.satterthwaite.aov
## Computes prediction intervals using the Satterthwaite method to approximate the distribution of the RHS in (9) (see paper)
t.satterthwaite.aov <- function(data){
	Y <- data$Y
	Z <- data$Z
	n <- length(Y)
	n_i <- colSums(Z)
	
	Ybar <- mean(Y)
	stat <- aov.statistics(data)
	S <- stat$S.L
	lambda <- stat$lambda.L
	r <- stat$r.L
	L <- length(S)

	groups <- rep(0, sum(n_i))
	group = 1
	for(j in 1:length(groups)){
		group = min(which(j <= descum))
		groups[j] =  group
	}

	model <- lmer(Y ~ 1+(1|groups))
	summ <- summary(model)
	se2.hat <- (summ$sigma^2)
	sa2.hat <- summ$varcor[1]$groups[1]

	std.err2 <- (se2.hat / n) + (sa2.hat *(1+sum(n_i^2)/(n^2)))
	std.err3 <- (se2.hat * (1+ 1/n)) + (sa2.hat *(1+sum(n_i^2)/(n^2)))
	std.err4 <- (se2.hat * (1+ 1/n)) + (sa2.hat *(1+sum(n_i^2)/(n^2)) - 2*n_i[L]/n)
	std.err <- (sa2.hat * sum(r*lambda) + se2.hat*sum(r))
	nu.hat <- ((sa2.hat^2)*(sum(r*lambda)^2) + 2*(sa2.hat*se2.hat)*sum(r*lambda)*sum(r) + ((se2.hat^2)*(sum(r)^2)))/(sa2.hat*sum(r*(lambda^2)) + 2*sa2.hat*se2.hat*sum(r*lambda) + se2.hat*sum(r))
	
	me.95 <- qt(0.975, nu.hat)*sqrt(std.err2*sum(S)/std.err)
	pi.95 <- c(Ybar - me.95, Ybar + me.95)

	me.90 <- qt(0.95, nu.hat)*sqrt(std.err2*sum(S)/std.err)
	pi.90 <- c(Ybar - me.90, Ybar + me.90)

	me.80 <- qt(0.90, nu.hat)*sqrt(std.err2*sum(S)/std.err)
	pi.80 <- c(Ybar - me.80, Ybar + me.80)

	me.95 <- qt(0.975, nu.hat)*sqrt(std.err3*sum(S)/std.err)
	pi.new.95 <- c(Ybar - me.95, Ybar + me.95)

	me.90 <- qt(0.95, nu.hat)*sqrt(std.err3*sum(S)/std.err)
	pi.new.90 <- c(Ybar - me.90, Ybar + me.90)

	me.80 <- qt(0.90, nu.hat)*sqrt(std.err3*sum(S)/std.err)
	pi.new.80 <- c(Ybar - me.80, Ybar + me.80)

	me.95 <- qt(0.975, nu.hat)*sqrt(std.err4*sum(S)/std.err)
	pi.exs.95 <- c(Ybar - me.95, Ybar + me.95)

	me.90 <- qt(0.95, nu.hat)*sqrt(std.err4*sum(S)/std.err)
	pi.exs.90 <- c(Ybar - me.90, Ybar + me.90)

	me.80 <- qt(0.90, nu.hat)*sqrt(std.err4*sum(S)/std.err)
	pi.exs.80 <- c(Ybar - me.80, Ybar + me.80)

	return(list(pi.95 = pi.95, pi.90 = pi.90, pi.80 = pi.80, pi.new.95 = pi.new.95, pi.new.90 = pi.new.90, pi.new.80 = pi.new.80, pi.exs.95 = pi.exs.95, pi.exs.90 = pi.exs.90, pi.exs.80 = pi.exs.80))

}




## t.gen.satterthwaite.aov
## Computes prediction intervals using the generalized Satterthwaite method (see Francq et al., 2019, cited in paper)
t.gen.satterthwaite.aov <- function(data){
	Y <- data$Y
	Z <- data$Z
	n <- length(Y)
	n_i <- colSums(Z)
	
	Ybar <- mean(Y)
	stat <- aov.statistics(data)
	S <- stat$S.L
	lambda <- stat$lambda.L
	r <- stat$r.L
	L <- length(S)

	groups <- rep(0, sum(n_i))
	group = 1
	for(j in 1:length(groups)){
		group = min(which(j <= descum))
		groups[j] =  group
	}

	model <- lmer(Y ~ 1+(1|groups))
	summ <- summary(model)
	se2.hat <- (summ$sigma^2)
	sa2.hat <- summ$varcor[1]$groups[1]
	model <- lmer(Y ~ 1+(1|groups), REML = FALSE)
	mat <-  vcov(model, information = 'expected', full = TRUE)

	std.err2 <- (se2.hat / n) + (sa2.hat *(1+sum(n_i^2)/(n^2)))
	std.err3 <- (se2.hat * (1+ 1/n)) + (sa2.hat *(1+sum(n_i^2)/(n^2)))
	std.err4 <- (se2.hat * (1+ 1/n)) + (sa2.hat *(1+sum(n_i^2)/(n^2)) - 2*n_i[L]/n)
	std.err <- (sa2.hat * sum(r*lambda) + se2.hat*sum(r))
	nu.hat <- 2*((se2.hat + sa2.hat)^2)/sum(mat[-1,-1])
	
	me.95 <- qt(0.975, nu.hat)*sqrt(std.err2*sum(S)/std.err)
	pi.95 <- c(Ybar - me.95, Ybar + me.95)

	me.90 <- qt(0.95, nu.hat)*sqrt(std.err2*sum(S)/std.err)
	pi.90 <- c(Ybar - me.90, Ybar + me.90)

	me.80 <- qt(0.90, nu.hat)*sqrt(std.err2*sum(S)/std.err)
	pi.80 <- c(Ybar - me.80, Ybar + me.80)

	me.95 <- qt(0.975, nu.hat)*sqrt(std.err3*sum(S)/std.err)
	pi.new.95 <- c(Ybar - me.95, Ybar + me.95)

	me.90 <- qt(0.95, nu.hat)*sqrt(std.err3*sum(S)/std.err)
	pi.new.90 <- c(Ybar - me.90, Ybar + me.90)

	me.80 <- qt(0.90, nu.hat)*sqrt(std.err3*sum(S)/std.err)
	pi.new.80 <- c(Ybar - me.80, Ybar + me.80)

	me.95 <- qt(0.975, nu.hat)*sqrt(std.err4*sum(S)/std.err)
	pi.exs.95 <- c(Ybar - me.95, Ybar + me.95)

	me.90 <- qt(0.95, nu.hat)*sqrt(std.err4*sum(S)/std.err)
	pi.exs.90 <- c(Ybar - me.90, Ybar + me.90)

	me.80 <- qt(0.90, nu.hat)*sqrt(std.err4*sum(S)/std.err)
	pi.exs.80 <- c(Ybar - me.80, Ybar + me.80)

	return(list(pi.95 = pi.95, pi.90 = pi.90, pi.80 = pi.80, pi.new.95 = pi.new.95, pi.new.90 = pi.new.90, pi.new.80 = pi.new.80, pi.exs.95 = pi.exs.95, pi.exs.90 = pi.exs.90, pi.exs.80 = pi.exs.80))

}


## pi.bootMer.aov
## Computes prediction intervals for random intercept model using parametric bootstrap as implemented in bootMer and the lme4 package
pi.bootMer.aov <- function(data){
	Y <- data$Y
	Z <- data$Z
	n <- length(Y)
	n_i <- colSums(Z)
	L <- length(n_i)
	c1.theta <- 1+(sum(n_i^2)/(n^2))
	c1.new <- 1+(sum(n_i^2)/(n^2))
	c1.exs <- 1+(sum(n_i^2)/(n^2))-(2*n_i[L]/n)
	c2.theta <- 1/n
	c2.new <- (1/n)+1
	c2.exs <- (1/n)+1

	groups <- rep(0, sum(n_i))
	group = 1
	for(j in 1:length(groups)){
		group = min(which(j <= descum))
		groups[j] =  group
	}

	model <- lmer(Y ~ 1+(1|groups))

	fun <- function(model){
			summ <- summary(model)
			beta.hat <- summ$coefficients[1]
			se2.hat <- (summ$sigma^2)
			sa2.hat <- summ$varcor[1]$groups[1]
			z <- rnorm(1)
			return(c(z*sqrt(c1.theta*sa2.hat+c2.theta*se2.hat)+beta.hat, z*sqrt(c1.new*sa2.hat+c2.new*se2.hat)+beta.hat, z*sqrt(c1.exs*sa2.hat+c2.exs*se2.hat)+beta.hat, z*sqrt(sa2.hat)+beta.hat))	
	}	

	booted <- bootMer(model, fun, type = 'parametric', nsim = 100)
	pi.95 <- as.numeric(quantile(booted$t[,1], c(0.025, 0.975)))
	pi.90 <- as.numeric(quantile(booted$t[,1], c(0.05, 0.95)))
	pi.80 <- as.numeric(quantile(booted$t[,1], c(0.1, 0.9)))

	pi.new.95 <- as.numeric(quantile(booted$t[,2], c(0.025, 0.975)))
	pi.new.90 <- as.numeric(quantile(booted$t[,2], c(0.05, 0.95)))
	pi.new.80 <- as.numeric(quantile(booted$t[,2], c(0.1, 0.9)))

	pi.exs.95 <- as.numeric(quantile(booted$t[,3], c(0.025, 0.975)))
	pi.exs.90 <- as.numeric(quantile(booted$t[,3], c(0.05, 0.95)))
	pi.exs.80 <- as.numeric(quantile(booted$t[,3], c(0.1, 0.9)))

	pi.95.a <- as.numeric(quantile(booted$t[,4], c(0.025, 0.975)))
	pi.90.a <- as.numeric(quantile(booted$t[,4], c(0.05, 0.95)))
	pi.80.a <- as.numeric(quantile(booted$t[,4], c(0.1, 0.9)))

	return(list(pi.95 = pi.95, pi.90 = pi.90, pi.80 = pi.80, pi.new.95 = pi.new.95, pi.new.90 = pi.new.90, pi.new.80 = pi.new.80, pi.exs.95 = pi.exs.95, pi.exs.90 = pi.exs.90, pi.exs.80 = pi.exs.80, pi.95.a = pi.95.a, pi.90.a = pi.90.a, pi.80.a = pi.80.a))

}




