


#source("IM.R")


######################################################################################### SIMULATING DATA

### 
#
#	Notes:  The following simulation settings are used in the paper
#		1. mu = 0, sigma_a^2 = 0.1, sigma^2 = 1.0, balanced design 5 groups of 6 
#		2. mu = 0, sigma_a^2 = 0.5, sigma^2 = 0.5, balanced design 5 groups of 6
#		3. mu = 0, sigma_a^2 = 1.0, sigma^2 = 0.1, balanced design 5 groups of 6
#		4. mu = 0, sigma_a^2 = 0.1, sigma^2 = 1.0, balanced design 10 groups of 12
#		5. mu = 0, sigma_a^2 = 0.5, sigma^2 = 0.5, balanced design 10 groups of 12
#		6. mu = 0, sigma_a^2 = 1.0, sigma^2 = 0.1, balanced design 10 groups of 12
#		7. mu = 0, sigma_a^2 = 0.1, sigma^2 = 1.0, unbalanced design 5 groups of (4,4,4,6,12)
#		8. mu = 0, sigma_a^2 = 0.5, sigma^2 = 0.5, unbalanced design 5 groups of (4,4,4,6,12)
#		9. mu = 0, sigma_a^2 = 1.0, sigma^2 = 0.1, unbalanced design 5 groups of (4,4,4,6,12)	
#	     10. mu = 0, sigma_a^2 = 0.1, sigma^2 = 1.0, unbalanced design 10 groups of (4,4,7,11,13,16,16,16,16,17)
#	     11. mu = 0, sigma_a^2 = 0.5, sigma^2 = 0.5, unbalanced design 10 groups of (4,4,7,11,13,16,16,16,16,17)
#	     12. mu = 0, sigma_a^2 = 1.0, sigma^2 = 0.1, unbalanced design 10 groups of (4,4,7,11,13,16,16,16,16,17)
#
###

# setting the seed will not result in perfect replication of simulation results
# because it does not affect random number generation within c++ functions used by impred package 
set.seed(852021)
# Number of data sets to be simulated
K <- 2000
# Design matrix for ANOVA model including two observations to be used as "true" predicted values
# one additional value is within the last group and one additional value is in a new group
des <- c(6,6,6,6,6)
I <- length(des)
des.pred <- c(6,6,6,6,7,1)
Z <- make.Z(des)
Z.1 <- make.Z(des.pred) 
n.1 <- sum(des.pred)
# number of observations is two less than number of generated responses
n <- n.1-2  
# design matrix of the observed data
Z <- Z.1[1:n,1:(ncol(Z.1)-1)]	
# across/between group variance
sigma2_a <- 0.1
# across/between response variance
sigma2 <- 1.0
# mean parameter	
mu <- 0	
# matrix of simulated responses and matrix of simulated new group means
Y.mat <- matrix(NA,n.1,K)	
theta.new <- matrix(NA,K,1)
for(i in 1:K){
	Y.mat[,i] <- sample.Y(Z.1, sqrt(sigma2_a), sqrt(sigma2), mu)
	theta.new[i,] <- rnorm(1, mu, sqrt(sigma2_a))
}
# matrix of responses without held-out new values to be predicted
Y <- Y.mat[1:n,]	
# matrix of values representing "true" predictions, first column are within existing group, last column are new group
Y.new <- t(Y.mat[(n+1):n.1,])	


########################################################################################## CONSTRUCTING PREDICTIONS

#### MLE estimation functions used for Higgins-type plug-in Student t prediction intervals

neg.log.lik <- function(par, y, Z){
	n<-length(y)
	mu <- par[1]; sigma2_a <- exp(par[2]); sigma2 <- exp(par[3])
	Sigma <- sigma2_a*(Z%*%t(Z))+diag(sigma2, n)
	Sigma.inv <- solve(Sigma)
	M<-matrix(y-rep(mu,n), n, 1)
	return(log(abs(det(Sigma)))+(t(M)%*%Sigma.inv%*%M))
}

MLE <- function(y, Z){
	par <- c(mu, sigma2_a, sigma2)
	return(optim(par, neg.log.lik, y=y, Z=Z))
} 


#### Functions for conformal prediction

pred.lvl.new <- (1/(n+1))*floor(0.90*(n+1))
pred.alpha.new <- (1-pred.lvl.new)/2

n.I <- des[I]
pred.lvl.exs <- (1/(n.I+1))*floor(0.90*(n.I+1))
pred.alpha.exs <- (1-pred.lvl.exs)/2

pred.lvl.theta <- (1/(I+1))*floor(0.90*(I+1))
pred.alpha.theta <- (1-pred.lvl.theta)/2


nonconformity <- function(y, x) abs(mean(y)-x)
nonconformity.i <- function(i, y){
	x<-y[i]
	y<-y[-i]
	return(nonconformity(y,x))
}
nonconformity.all <- function(y){
	n<-length(y)
	return(apply(matrix(1:n,n,1),1,nonconformity, y=y))
}
nonconformity.plaus <- function(x, y){
	yx<-c(y,x)
	n<-length(yx)
	nonconformities <- nonconformity.all(yx)
	nonconformity.x <- nonconformity(y,x)
	plaus <- (1/n)*sum(nonconformities>=nonconformity.x)
}
nonconformity.plaus.grid <- function(grid, y){
	return(apply(matrix(grid, length(grid), 1), 1, nonconformity.plaus, y=y))
}




#### Looping through K simulated data sets
results <- list( oracle = matrix(NA, K, 18) , stdt = matrix(NA, K, 18), im = matrix(NA, K, 18), conf = matrix(NA, K, 6)  )
for(k in 1:K){

	#### The data and quantities to be predicted for the kth run
	response <- Y[,k]
	response.last.group <- response[(n-des[I]+1):n]
	response.group.averages <- (t(Z)%*%matrix(response,n,1))/des
	mean.response <- mean(response)
	mean.response.last.group <- mean(response.last.group)
	mean.averages <- mean(response.group.averages)
	Y.star.exs <- Y.new[k,1]
	Y.star.new <- Y.new[k,2]
	theta.star <- theta.new[k,1]

	#### Predictions using "Oracle" --- true model parameters
	oracle.std.dev.new <- sqrt((sigma2_a)*(1+(1/(n^2))*sum(des^2)) + (sigma2)*(1/n+1/1)) 
	oracle.std.dev.exs <- sqrt((sigma2_a)*(1+(1/(n^2))*sum(des^2)) + (sigma2)*(1/n+1/1)) 

	oracle.pdi.new.80 <-  c(mu + qnorm(0.10)*oracle.std.dev.new, mu + qnorm(0.90)*oracle.std.dev.new)
	oracle.pdi.new.90 <-  c(mu + qnorm(0.05)*oracle.std.dev.new, mu + qnorm(0.95)*oracle.std.dev.new)
	oracle.pdi.new.95 <-  c(mu + qnorm(0.025)*oracle.std.dev.new, mu + qnorm(0.975)*oracle.std.dev.new)
	
	results$oracle[k,1] <- ifelse(oracle.pdi.new.80[1]<=Y.star.new & oracle.pdi.new.80[2]>=Y.star.new, 1, 0)
	results$oracle[k,2] <- oracle.pdi.new.80[2]-oracle.pdi.new.80[1]
	results$oracle[k,3] <- ifelse(oracle.pdi.new.90[1]<=Y.star.new & oracle.pdi.new.90[2]>=Y.star.new, 1, 0)
	results$oracle[k,4] <- oracle.pdi.new.90[2]-oracle.pdi.new.90[1]
	results$oracle[k,5] <- ifelse(oracle.pdi.new.95[1]<=Y.star.new & oracle.pdi.new.95[2]>=Y.star.new, 1, 0)
	results$oracle[k,6] <- oracle.pdi.new.95[2]-oracle.pdi.new.95[1]

	oracle.pdi.exs.80 <-  c(mu + qnorm(0.10)*oracle.std.dev.exs, mu + qnorm(0.90)*oracle.std.dev.exs)
	oracle.pdi.exs.90 <-  c(mu + qnorm(0.05)*oracle.std.dev.exs, mu + qnorm(0.95)*oracle.std.dev.exs)
	oracle.pdi.exs.95 <-  c(mu + qnorm(0.025)*oracle.std.dev.exs, mu + qnorm(0.975)*oracle.std.dev.exs)

	results$oracle[k,7] <- ifelse(oracle.pdi.exs.80[1]<=Y.star.exs & oracle.pdi.exs.80[2]>=Y.star.exs, 1, 0)
	results$oracle[k,8] <- oracle.pdi.exs.80[2]-oracle.pdi.exs.80[1]
	results$oracle[k,9] <- ifelse(oracle.pdi.exs.90[1]<=Y.star.exs & oracle.pdi.exs.90[2]>=Y.star.exs, 1, 0)
	results$oracle[k,10] <- oracle.pdi.exs.90[2]-oracle.pdi.exs.90[1]
	results$oracle[k,11] <- ifelse(oracle.pdi.exs.95[1]<=Y.star.exs & oracle.pdi.exs.95[2]>=Y.star.exs, 1, 0)
	results$oracle[k,12] <- oracle.pdi.exs.95[2]-oracle.pdi.exs.95[1]

	oracle.pdi.theta.80 <-  c(mu + qnorm(0.10)*sqrt(sigma2_a), mu + qnorm(0.90)*sqrt(sigma2_a))
	oracle.pdi.theta.90 <-  c(mu + qnorm(0.05)*sqrt(sigma2_a), mu + qnorm(0.95)*sqrt(sigma2_a))
	oracle.pdi.theta.95 <-  c(mu + qnorm(0.025)*sqrt(sigma2_a), mu + qnorm(0.975)*sqrt(sigma2_a))

	results$oracle[k,13] <- ifelse(oracle.pdi.theta.80[1]<=theta.star & oracle.pdi.theta.80[2]>=theta.star, 1, 0)
	results$oracle[k,14] <- oracle.pdi.theta.80[2]-oracle.pdi.theta.80[1]
	results$oracle[k,15] <- ifelse(oracle.pdi.theta.90[1]<=theta.star & oracle.pdi.theta.90[2]>=theta.star, 1, 0)
	results$oracle[k,16] <- oracle.pdi.theta.90[2]-oracle.pdi.theta.90[1]
	results$oracle[k,17] <- ifelse(oracle.pdi.theta.95[1]<=theta.star & oracle.pdi.theta.95[2]>=theta.star, 1, 0)
	results$oracle[k,18] <- oracle.pdi.theta.95[2]-oracle.pdi.theta.95[1]


	#### Predictions using Higgins-type plug-in Student t prediction intervals 

	mles <- MLE(response,Z)
	stdt.pdi.new.80 <- qt(c(0.10, 0.90), df = n-1)*sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n+1))+mles$par[1]
	stdt.pdi.new.90 <- qt(c(0.05, 0.95), df = n-1)*sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n+1))+mles$par[1]
	stdt.pdi.new.95 <- qt(c(0.025, 0.975), df = n-1)*sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n+1))+mles$par[1]	

	results$stdt[k,1] <- ifelse(stdt.pdi.new.80[1]<=Y.star.new & stdt.pdi.new.80[2]>=Y.star.new, 1, 0)
	results$stdt[k,2] <- stdt.pdi.new.80[2]-stdt.pdi.new.80[1]
	results$stdt[k,3] <- ifelse(stdt.pdi.new.90[1]<=Y.star.new & stdt.pdi.new.90[2]>=Y.star.new, 1, 0)
	results$stdt[k,4] <- stdt.pdi.new.90[2]-stdt.pdi.new.90[1]
	results$stdt[k,5] <- ifelse(stdt.pdi.new.95[1]<=Y.star.new & stdt.pdi.new.95[2]>=Y.star.new, 1, 0)
	results$stdt[k,6] <- stdt.pdi.new.95[2]-stdt.pdi.new.95[1]

	stdt.pdi.exs.80 <- qt(c(0.100, 0.900), df = n-1)*sqrt(exp(mles$par[2])*(1-2*(I/n)+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n+1))+mles$par[1]
	stdt.pdi.exs.90 <- qt(c(0.050, 0.950), df = n-1)*sqrt(exp(mles$par[2])*(1-2*(I/n)+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n+1))+mles$par[1]
	stdt.pdi.exs.95 <- qt(c(0.025, 0.975), df = n-1)*sqrt(exp(mles$par[2])*(1-2*(I/n)+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n+1))+mles$par[1]

	results$stdt[k,7] <- ifelse(stdt.pdi.exs.80[1]<=Y.star.exs & stdt.pdi.exs.80[2]>=Y.star.exs, 1, 0)
	results$stdt[k,8] <- stdt.pdi.exs.80[2]-stdt.pdi.exs.80[1]
	results$stdt[k,9] <- ifelse(stdt.pdi.exs.90[1]<=Y.star.exs & stdt.pdi.exs.90[2]>=Y.star.exs, 1, 0)
	results$stdt[k,10] <- stdt.pdi.exs.90[2]-stdt.pdi.exs.90[1]
	results$stdt[k,11] <- ifelse(stdt.pdi.exs.95[1]<=Y.star.exs & stdt.pdi.exs.95[2]>=Y.star.exs, 1, 0)
	results$stdt[k,12] <- stdt.pdi.exs.95[2]-stdt.pdi.exs.95[1]

	stdt.pdi.theta.80 <- qt(c(0.100, 0.900), df = n-1)*sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n))+mles$par[1]
	stdt.pdi.theta.90 <- qt(c(0.050, 0.950), df = n-1)*sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n))+mles$par[1]
	stdt.pdi.theta.95 <- qt(c(0.025, 0.975), df = n-1)*sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(des^2))+exp(mles$par[3])*(1/n))+mles$par[1]

	results$stdt[k,13] <- ifelse(stdt.pdi.theta.80[1]<=theta.star & stdt.pdi.theta.80[2]>=theta.star, 1, 0)
	results$stdt[k,14] <- stdt.pdi.theta.80[2]-stdt.pdi.theta.80[1]
	results$stdt[k,15] <- ifelse(stdt.pdi.theta.90[1]<=theta.star & stdt.pdi.theta.90[2]>=theta.star, 1, 0)
	results$stdt[k,16] <- stdt.pdi.theta.90[2]-stdt.pdi.theta.90[1]
	results$stdt[k,17] <- ifelse(stdt.pdi.theta.95[1]<=theta.star & stdt.pdi.theta.95[2]>=theta.star, 1, 0)
	results$stdt[k,18] <- stdt.pdi.theta.95[2]-stdt.pdi.theta.95[1]


	#### Predictions using Inferential Model

	ex.data <- list(Y=response, Z=Z)
	ex.statistics <- aov.statistics(ex.data)
	ex.grid.size <- 500
	ex.grid <- matrix(seq(from = min(ex.data$Y)-6, to = max(ex.data$Y)+6, length.out = ex.grid.size),ex.grid.size,1)
	sig.grid = as.matrix(expand.grid(seq(from = 0.01, to = 7, length.out = 14), seq(from = 0.01, to = 7, length.out = 14)))
	ex.grid.plauses <-  apply.plaus(sig.grid, matrix(runif(10000),10000,1), 1, ex.data, ex.statistics, ex.grid)
	
	plaus.max.new <- which.max(ex.grid.plauses$fused.plaus.n)	
	im.pdi.new.80 <- c(ex.grid[which.min(abs(0.100-ex.grid.plauses$fused.plaus.n[1:plaus.max.new]))], ex.grid[which.min(abs(0.900-ex.grid.plauses$fused.plaus.n[(plaus.max.new+1):ex.grid.size]))])
	im.pdi.new.90 <- c(ex.grid[which.min(abs(0.050-ex.grid.plauses$fused.plaus.n[1:plaus.max.new]))], ex.grid[which.min(abs(0.950-ex.grid.plauses$fused.plaus.n[(plaus.max.new+1):ex.grid.size]))])
	im.pdi.new.95 <- c(ex.grid[which.min(abs(0.025-ex.grid.plauses$fused.plaus.n[1:plaus.max.new]))], ex.grid[which.min(abs(0.975-ex.grid.plauses$fused.plaus.n[(plaus.max.new+1):ex.grid.size]))])

	results$im[k,1] <- ifelse(im.pdi.new.80[1]<=Y.star.new & im.pdi.new.80[2]>=Y.star.new, 1, 0)
	results$im[k,2] <- im.pdi.new.80[2]-im.pdi.new.80[1]
	results$im[k,3] <- ifelse(im.pdi.new.90[1]<=Y.star.new & im.pdi.new.90[2]>=Y.star.new, 1, 0)
	results$im[k,4] <- im.pdi.new.90[2]-im.pdi.new.90[1]
	results$im[k,5] <- ifelse(im.pdi.new.95[1]<=Y.star.new & im.pdi.new.95[2]>=Y.star.new, 1, 0)
	results$im[k,6] <- im.pdi.new.95[2]-im.pdi.new.95[1]

	plaus.max.exs <- which.max(ex.grid.plauses$fused.plaus.w)	
	im.pdi.exs.80 <- c(ex.grid[which.min(abs(0.100-ex.grid.plauses$fused.plaus.w[1:plaus.max.exs]))], ex.grid[which.min(abs(0.900-ex.grid.plauses$fused.plaus.w[(plaus.max.exs+1):ex.grid.size]))])
	im.pdi.exs.90 <- c(ex.grid[which.min(abs(0.050-ex.grid.plauses$fused.plaus.w[1:plaus.max.exs]))], ex.grid[which.min(abs(0.950-ex.grid.plauses$fused.plaus.w[(plaus.max.exs+1):ex.grid.size]))])
	im.pdi.exs.95 <- c(ex.grid[which.min(abs(0.025-ex.grid.plauses$fused.plaus.w[1:plaus.max.exs]))], ex.grid[which.min(abs(0.975-ex.grid.plauses$fused.plaus.w[(plaus.max.exs+1):ex.grid.size]))])

	results$im[k,7] <- ifelse(im.pdi.exs.80[1]<=Y.star.exs & im.pdi.exs.80[2]>=Y.star.exs, 1, 0)
	results$im[k,8] <- im.pdi.exs.80[2]-im.pdi.exs.80[1]
	results$im[k,9] <- ifelse(im.pdi.exs.90[1]<=Y.star.exs & im.pdi.exs.90[2]>=Y.star.exs, 1, 0)
	results$im[k,10] <- im.pdi.exs.90[2]-im.pdi.exs.90[1]
	results$im[k,11] <- ifelse(im.pdi.exs.95[1]<=Y.star.exs & im.pdi.exs.95[2]>=Y.star.exs, 1, 0)
	results$im[k,12] <- im.pdi.exs.95[2]-im.pdi.exs.95[1]

	plaus.max.theta <- which.max(ex.grid.plauses$fused.plaus.T)	
	im.pdi.theta.80 <- c(ex.grid[which.min(abs(0.100-ex.grid.plauses$fused.plaus.T[1:plaus.max.theta]))], ex.grid[which.min(abs(0.900-ex.grid.plauses$fused.plaus.T[(plaus.max.theta+1):ex.grid.size]))])
	im.pdi.theta.90 <- c(ex.grid[which.min(abs(0.050-ex.grid.plauses$fused.plaus.T[1:plaus.max.theta]))], ex.grid[which.min(abs(0.950-ex.grid.plauses$fused.plaus.T[(plaus.max.theta+1):ex.grid.size]))])
	im.pdi.theta.95 <- c(ex.grid[which.min(abs(0.025-ex.grid.plauses$fused.plaus.T[1:plaus.max.theta]))], ex.grid[which.min(abs(0.975-ex.grid.plauses$fused.plaus.T[(plaus.max.theta+1):ex.grid.size]))])

	results$im[k,13] <- ifelse(im.pdi.theta.80[1]<=theta.star & im.pdi.theta.80[2]>=theta.star, 1, 0)
	results$im[k,14] <- im.pdi.theta.80[2]-im.pdi.theta.80[1]
	results$im[k,15] <- ifelse(im.pdi.theta.90[1]<=theta.star & im.pdi.theta.90[2]>=theta.star, 1, 0)
	results$im[k,16] <- im.pdi.theta.90[2]-im.pdi.theta.90[1]
	results$im[k,17] <- ifelse(im.pdi.theta.95[1]<=theta.star & im.pdi.theta.95[2]>=theta.star, 1, 0)
	results$im[k,18] <- im.pdi.theta.95[2]-im.pdi.theta.95[1]


	#### Conformal Prediction

	grid.new.size <- 480
	grid.new <- seq(from = mean.response	 - 6, to = mean.response + 6, length.out = grid.new.size)
	conformal.plaus.new <- nonconformity.plaus.grid(grid.new, response)
	conf.plaus.max.new <- which.max(conformal.plaus.new)	 
	conf.pdi.new <- c(grid.new[which.min(abs(pred.alpha.new-conformal.plaus.new[1:conf.plaus.max.new]))], grid.new[which.min(abs((1-pred.alpha.new)-conformal.plaus.new[(conf.plaus.max.new+1):grid.new.size]))])

	results$conf[k,1] <- ifelse(conf.pdi.new[1]<=Y.star.new & conf.pdi.new[2]>=Y.star.new, 1, 0)
	results$conf[k,2] <- conf.pdi.new[2]-conf.pdi.new[1]

	grid.exs.size <- 480
	grid.exs <- seq(from = mean.response.last.group	 - 6, to = mean.response.last.group + 6, length.out = grid.exs.size)
	conformal.plaus.exs <- nonconformity.plaus.grid(grid.exs, response.last.group)
	conf.plaus.max.exs <- which.max(conformal.plaus.exs)	
	conf.pdi.exs <- c(grid.exs[which.min(abs(pred.alpha.exs-conformal.plaus.exs[1:conf.plaus.max.exs]))], grid.exs[which.min(abs((1-pred.alpha.exs)-conformal.plaus.exs[(conf.plaus.max.exs+1):grid.exs.size]))])

	results$conf[k,3] <- ifelse(conf.pdi.exs[1]<=Y.star.exs & conf.pdi.exs[2]>=Y.star.exs, 1, 0)
	results$conf[k,4] <- conf.pdi.exs[2]-conf.pdi.exs[1]

	grid.theta.size <- 480
	grid.theta <- seq(from = mean.averages	 - 6, to = mean.averages + 6, length.out = grid.theta.size)
	conformal.plaus.theta <- nonconformity.plaus.grid(grid.theta, response.group.averages)
	conf.plaus.max.theta <- which.max(conformal.plaus.theta)	
	conf.pdi.theta <- c(grid.theta[which.min(abs(pred.alpha.theta-conformal.plaus.theta[1:conf.plaus.max.theta]))], grid.theta[which.min(abs((1-pred.alpha.theta)-conformal.plaus.theta[(conf.plaus.max.theta+1):grid.theta.size]))])

	results$conf[k,5] <- ifelse(conf.pdi.theta[1]<=theta.star & conf.pdi.theta[2]>=theta.star, 1, 0)
	results$conf[k,6] <- conf.pdi.theta[2]-conf.pdi.theta[1]


	if(k%%10==0) print(k)
}
