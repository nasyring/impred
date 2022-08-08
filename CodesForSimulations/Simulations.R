
# Packages needed
library(impred)
library(predintma)
library(lme4)
library(rstanarm)
#need to load the functions in the IM.R script
#source("IM.R")  

######################################################################################### SIMULATING DATA

###  IT IS NECESSARY TO MANUALLY CHANGE THE SIMULATION SETTINGS BELOW
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
K <- 1000
# Design matrix for ANOVA model including two observations to be used as "true" predicted values
# one additional value is within the last group and one additional value is in a new group
des <- c(6,6,6,6,6)
descum <- cumsum(des)
groups <- rep(0, sum(des))
group = 1
for(j in 1:length(groups)){
	group = min(which(j <= descum))
	groups[j] =  group
}
I <- length(des)
des.pred <- c(6,6,6,6,7,1)
Z <- make.Z(des)
Z.1 <- make.Z(des.pred) 
n.1 <- sum(des.pred)
# number of observations is two less than number of generated responses
n <- n.1-2  
# design matrix of the observed data
Z <- Z.1[1:n,1:(ncol(Z.1)-1)]	
# between group variance
sigma2_a <- 0.01
# within group variance
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


#### Looping through K simulated data sets
results <- list( oracle = matrix(NA, K, 18) , stdt = matrix(NA, K, 18), im = matrix(NA, K, 18), boot = matrix(NA, K, 18), para.boot = matrix(NA, K, 18), bayes.cauchy = matrix(NA, K, 18)  )
for(k in 1:K){
	then = proc.time()

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
	 oracle.std.dev.exs <- sqrt((sigma2_a)*(1-2*des[length(des)]/n+(1/(n^2))*sum(des^2)) + (sigma2)*(1/n+1/1))
	 oracle.std.dev.theta <- sqrt((sigma2_a)*(1+(1/(n^2))*sum(des^2)) + (sigma2)*(1/n))
	
	 oracle.pdi.new.80 <-  c(mean.response + qnorm(0.10)*oracle.std.dev.new, mean.response + qnorm(0.90)*oracle.std.dev.new)
	 oracle.pdi.new.90 <-  c(mean.response + qnorm(0.05)*oracle.std.dev.new, mean.response + qnorm(0.95)*oracle.std.dev.new)
	 oracle.pdi.new.95 <-  c(mean.response + qnorm(0.025)*oracle.std.dev.new, mean.response + qnorm(0.975)*oracle.std.dev.new)
	
	 results$oracle[k,1] <- ifelse(oracle.pdi.new.80[1]<=Y.star.new & oracle.pdi.new.80[2]>=Y.star.new, 1, 0)
	 results$oracle[k,2] <- oracle.pdi.new.80[2]-oracle.pdi.new.80[1]
	 results$oracle[k,3] <- ifelse(oracle.pdi.new.90[1]<=Y.star.new & oracle.pdi.new.90[2]>=Y.star.new, 1, 0)
	 results$oracle[k,4] <- oracle.pdi.new.90[2]-oracle.pdi.new.90[1]
	 results$oracle[k,5] <- ifelse(oracle.pdi.new.95[1]<=Y.star.new & oracle.pdi.new.95[2]>=Y.star.new, 1, 0)
	 results$oracle[k,6] <- oracle.pdi.new.95[2]-oracle.pdi.new.95[1]
	
	 oracle.pdi.exs.80 <-  c(mean.response + qnorm(0.10)*oracle.std.dev.exs, mean.response + qnorm(0.90)*oracle.std.dev.exs)
	 oracle.pdi.exs.90 <-  c(mean.response + qnorm(0.05)*oracle.std.dev.exs, mean.response + qnorm(0.95)*oracle.std.dev.exs)
	 oracle.pdi.exs.95 <-  c(mean.response + qnorm(0.025)*oracle.std.dev.exs, mean.response + qnorm(0.975)*oracle.std.dev.exs)
	
	 results$oracle[k,7] <- ifelse(oracle.pdi.exs.80[1]<=Y.star.exs & oracle.pdi.exs.80[2]>=Y.star.exs, 1, 0)
	 results$oracle[k,8] <- oracle.pdi.exs.80[2]-oracle.pdi.exs.80[1]
	 results$oracle[k,9] <- ifelse(oracle.pdi.exs.90[1]<=Y.star.exs & oracle.pdi.exs.90[2]>=Y.star.exs, 1, 0)
	 results$oracle[k,10] <- oracle.pdi.exs.90[2]-oracle.pdi.exs.90[1]
	 results$oracle[k,11] <- ifelse(oracle.pdi.exs.95[1]<=Y.star.exs & oracle.pdi.exs.95[2]>=Y.star.exs, 1, 0)
	 results$oracle[k,12] <- oracle.pdi.exs.95[2]-oracle.pdi.exs.95[1]
	
	 oracle.pdi.theta.80 <-  c(mean.response + qnorm(0.10)*oracle.std.dev.theta, mean.response + qnorm(0.90)*oracle.std.dev.theta)
	 oracle.pdi.theta.90 <-  c(mean.response + qnorm(0.05)*oracle.std.dev.theta, mean.response + qnorm(0.95)*oracle.std.dev.theta)
	 oracle.pdi.theta.95 <-  c(mean.response + qnorm(0.025)*oracle.std.dev.theta, mean.response + qnorm(0.975)*oracle.std.dev.theta)
	
	 results$oracle[k,13] <- ifelse(oracle.pdi.theta.80[1]<=theta.star & oracle.pdi.theta.80[2]>=theta.star, 1, 0)

	 results$oracle[k,14] <- oracle.pdi.theta.80[2]-oracle.pdi.theta.80[1]
	 results$oracle[k,15] <- ifelse(oracle.pdi.theta.90[1]<=theta.star & oracle.pdi.theta.90[2]>=theta.star, 1, 0)
	 results$oracle[k,16] <- oracle.pdi.theta.90[2]-oracle.pdi.theta.90[1]
	 results$oracle[k,17] <- ifelse(oracle.pdi.theta.95[1]<=theta.star & oracle.pdi.theta.95[2]>=theta.star, 1, 0)
	 results$oracle[k,18] <- oracle.pdi.theta.95[2]-oracle.pdi.theta.95[1]



	#### Predictions using Inferential Model
	
	 ex.data <- list(Y=response, Z=Z)
	 ex.statistics <- aov.statistics(ex.data)
	 ex.grid.size <- 600
	 ex.grid <- matrix(seq(from = min(ex.data$Y)-8, to = max(ex.data$Y)+8, length.out = ex.grid.size),ex.grid.size,1)

	 if(all(des - mean(des) == 0)){  # balanced anova
		model <- lmer(response ~ 1+ (1|groups))
		fun <- function(model){
			summ <- summary(model)
			return(1/(1+ex.statistics$r.L[1] *(summ$varcor[1]$groups[1]/(summ$sigma^2))))
		}
		booted.ratio <- bootMer(model, fun, 100)
		boot.mean <- mean(booted.ratio$t)
		boot.sd <- sd(booted.ratio$t)
		var.grid <- seq(from = max(boot.mean-boot.sd, 0.01), to = min(boot.mean+boot.sd, 0.99), length.out = 3)
	 	plauses <- plaus_balanced_aov(ex.grid, mean(response), ex.statistics$S.L, ex.statistics$lambda.L, ex.statistics$r.L, n, des, var.grid)
	 }else {   # unbalanced anova design
		model <- lmer(response ~ 1+ (1|groups))
		fun <- function(model){
			summ <- summary(model)
			return((summ$sigma^2)/((summ$sigma^2)+summ$varcor[1]$groups[1]))
		}
		booted.ratio <- bootMer(model, fun, 100)
		boot.mean <- mean(booted.ratio$t)
		boot.sd <- sd(booted.ratio$t)
		ratio.grid <- seq(from = max(boot.mean-boot.sd, 0.01), to = min(boot.mean+boot.sd, 0.99), length.out = 3)
		plauses <- plaus_aov_unbalanced_full(ex.data, ex.grid, ratio.grid)
	 }

	 plaus.max.theta <- which.max(plauses$plauses.theta)
	 im.pdi.theta.80 <- c(ex.grid[which.min(abs(0.200-plauses$plauses.theta[1:plaus.max.theta]))], ex.grid[plaus.max.theta+which.min(abs(0.200-plauses$plauses.theta[(plaus.max.theta+1):ex.grid.size]))])
	 im.pdi.theta.90 <- c(ex.grid[which.min(abs(0.1-plauses$plauses.theta[1:plaus.max.theta]))], ex.grid[plaus.max.theta+which.min(abs(0.1-plauses$plauses.theta[(plaus.max.theta+1):ex.grid.size]))])
	 im.pdi.theta.95 <- c(ex.grid[which.min(abs(0.05-plauses$plauses.theta[1:plaus.max.theta]))], ex.grid[plaus.max.theta+which.min(abs(0.05-plauses$plauses.theta[(plaus.max.theta+1):ex.grid.size]))])


	 results$im[k,13] <- ifelse(im.pdi.theta.80[1]<=theta.star & im.pdi.theta.80[2]>=theta.star, 1, 0)
	 results$im[k,14] <- im.pdi.theta.80[2]-im.pdi.theta.80[1]
	 results$im[k,15] <- ifelse(im.pdi.theta.90[1]<=theta.star & im.pdi.theta.90[2]>=theta.star, 1, 0)
	 results$im[k,16] <- im.pdi.theta.90[2]-im.pdi.theta.90[1]
	 results$im[k,17] <- ifelse(im.pdi.theta.95[1]<=theta.star & im.pdi.theta.95[2]>=theta.star, 1, 0)
	 results$im[k,18] <- im.pdi.theta.95[2]-im.pdi.theta.95[1]
	

	 plaus.max.new <- which.max(plauses$plauses.new)
	 im.pdi.new.80 <- c(ex.grid[which.min(abs(0.200-plauses$plauses.new[1:plaus.max.new]))], ex.grid[plaus.max.new+which.min(abs(0.200-plauses$plauses.new[(plaus.max.new+1):ex.grid.size]))])
	 im.pdi.new.90 <- c(ex.grid[which.min(abs(0.1-plauses$plauses.new[1:plaus.max.new]))], ex.grid[plaus.max.new+which.min(abs(0.1-plauses$plauses.new[(plaus.max.new+1):ex.grid.size]))])
	 im.pdi.new.95 <- c(ex.grid[which.min(abs(0.05-plauses$plauses.new[1:plaus.max.new]))], ex.grid[plaus.max.new+which.min(abs(0.05-plauses$plauses.new[(plaus.max.new+1):ex.grid.size]))])

	 results$im[k,1] <- ifelse(im.pdi.new.80[1]<=Y.star.new & im.pdi.new.80[2]>=Y.star.new, 1, 0)
	 results$im[k,2] <- im.pdi.new.80[2]-im.pdi.new.80[1]
	 results$im[k,3] <- ifelse(im.pdi.new.90[1]<=Y.star.new & im.pdi.new.90[2]>=Y.star.new, 1, 0)
	 results$im[k,4] <- im.pdi.new.90[2]-im.pdi.new.90[1]
	 results$im[k,5] <- ifelse(im.pdi.new.95[1]<=Y.star.new & im.pdi.new.95[2]>=Y.star.new, 1, 0)
	 results$im[k,6] <- im.pdi.new.95[2]-im.pdi.new.95[1]


	 plaus.max.exs <- which.max(plauses$plauses.exs)
	 im.pdi.exs.80 <- c(ex.grid[which.min(abs(0.200-plauses$plauses.exs[1:plaus.max.exs]))], ex.grid[plaus.max.exs+which.min(abs(0.200-plauses$plauses.exs[(plaus.max.exs+1):ex.grid.size]))])
	 im.pdi.exs.90 <- c(ex.grid[which.min(abs(0.1-plauses$plauses.exs[1:plaus.max.exs]))], ex.grid[plaus.max.exs+which.min(abs(0.1-plauses$plauses.exs[(plaus.max.exs+1):ex.grid.size]))])
	 im.pdi.exs.95 <- c(ex.grid[which.min(abs(0.05-plauses$plauses.exs[1:plaus.max.exs]))], ex.grid[plaus.max.exs+which.min(abs(0.05-plauses$plauses.exs[(plaus.max.exs+1):ex.grid.size]))])

	 results$im[k,7] <- ifelse(im.pdi.exs.80[1]<=Y.star.exs & im.pdi.exs.80[2]>=Y.star.exs, 1, 0)
	 results$im[k,8] <- im.pdi.exs.80[2]-im.pdi.exs.80[1]
	 results$im[k,9] <- ifelse(im.pdi.exs.90[1]<=Y.star.exs & im.pdi.exs.90[2]>=Y.star.exs, 1, 0)
	 results$im[k,10] <- im.pdi.exs.90[2]-im.pdi.exs.90[1]
	 results$im[k,11] <- ifelse(im.pdi.exs.95[1]<=Y.star.exs & im.pdi.exs.95[2]>=Y.star.exs, 1, 0)
	 results$im[k,12] <- im.pdi.exs.95[2]-im.pdi.exs.95[1]

	#### Student's t Wald-type plug-in prediction intervals

	summ <- summary(model)
	summ$varcor[1]$groups[1]
	summ$sigma^2
  stdt.pdi.new.80 <- qt(c(0.10, 0.90), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n+1))+mean.response
  stdt.pdi.new.90 <- qt(c(0.05, 0.95), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n+1))+mean.response
  stdt.pdi.new.95 <- qt(c(0.025, 0.975), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n+1))+mean.response

  results$stdt[k,1] <- ifelse(stdt.pdi.new.80[1]<=Y.star.new & stdt.pdi.new.80[2]>=Y.star.new, 1, 0)
  results$stdt[k,2] <- stdt.pdi.new.80[2]-stdt.pdi.new.80[1]
  results$stdt[k,3] <- ifelse(stdt.pdi.new.90[1]<=Y.star.new & stdt.pdi.new.90[2]>=Y.star.new, 1, 0)
  results$stdt[k,4] <- stdt.pdi.new.90[2]-stdt.pdi.new.90[1]
  results$stdt[k,5] <- ifelse(stdt.pdi.new.95[1]<=Y.star.new & stdt.pdi.new.95[2]>=Y.star.new, 1, 0)
  results$stdt[k,6] <- stdt.pdi.new.95[2]-stdt.pdi.new.95[1]

  stdt.pdi.exs.80 <- qt(c(0.100, 0.900), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1-2*(I/n)+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n+1))+mean.response
  stdt.pdi.exs.90 <- qt(c(0.050, 0.950), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1-2*(I/n)+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n+1))+mean.response
  stdt.pdi.exs.95 <- qt(c(0.025, 0.975), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1-2*(I/n)+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n+1))+mean.response

  results$stdt[k,7] <- ifelse(stdt.pdi.exs.80[1]<=Y.star.exs & stdt.pdi.exs.80[2]>=Y.star.exs, 1, 0)
  results$stdt[k,8] <- stdt.pdi.exs.80[2]-stdt.pdi.exs.80[1]
  results$stdt[k,9] <- ifelse(stdt.pdi.exs.90[1]<=Y.star.exs & stdt.pdi.exs.90[2]>=Y.star.exs, 1, 0)
  results$stdt[k,10] <- stdt.pdi.exs.90[2]-stdt.pdi.exs.90[1]
  results$stdt[k,11] <- ifelse(stdt.pdi.exs.95[1]<=Y.star.exs & stdt.pdi.exs.95[2]>=Y.star.exs, 1, 0)
  results$stdt[k,12] <- stdt.pdi.exs.95[2]-stdt.pdi.exs.95[1]

  stdt.pdi.theta.80 <- qt(c(0.100, 0.900), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n))+mean.response
  stdt.pdi.theta.90 <- qt(c(0.050, 0.950), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n))+mean.response
  stdt.pdi.theta.95 <- qt(c(0.025, 0.975), df = n-1)*sqrt(summ$varcor[1]$groups[1]*(1+(1/(n^2))*sum(des^2))+(summ$sigma^2)*(1/n))+mean.response

  results$stdt[k,13] <- ifelse(stdt.pdi.theta.80[1]<=theta.star & stdt.pdi.theta.80[2]>=theta.star, 1, 0)
  results$stdt[k,14] <- stdt.pdi.theta.80[2]-stdt.pdi.theta.80[1]
  results$stdt[k,15] <- ifelse(stdt.pdi.theta.90[1]<=theta.star & stdt.pdi.theta.90[2]>=theta.star, 1, 0)
  results$stdt[k,16] <- stdt.pdi.theta.90[2]-stdt.pdi.theta.90[1]
  results$stdt[k,17] <- ifelse(stdt.pdi.theta.95[1]<=theta.star & stdt.pdi.theta.95[2]>=theta.star, 1, 0)
  results$stdt[k,18] <- stdt.pdi.theta.95[2]-stdt.pdi.theta.95[1]


	#### Parametric Bootstrap predictions 
	
	 ex.data <- list(Y=response, Z=Z)
	 t.pis <- pi.bootMer.aov(ex.data)

	 results$para.boot[k,13] <- ifelse(t.pis$pi.80[1]<=theta.star & t.pis$pi.80[2]>=theta.star, 1, 0)
	 results$para.boot[k,14] <- t.pis$pi.80[2]-t.pis$pi.80[1]
	 results$para.boot[k,15] <- ifelse(t.pis$pi.90[1]<=theta.star & t.pis$pi.90[2]>=theta.star, 1, 0)
	 results$para.boot[k,16] <- t.pis$pi.90[2]-t.pis$pi.90[1]
	 results$para.boot[k,17] <- ifelse(t.pis$pi.95[1]<=theta.star & t.pis$pi.95[2]>=theta.star, 1, 0)
	 results$para.boot[k,18] <- t.pis$pi.95[2]-t.pis$pi.95[1]
	
	 results$para.boot[k,1] <- ifelse(t.pis$pi.new.80[1]<=Y.star.new & t.pis$pi.new.80[2]>=Y.star.new, 1, 0)
	 results$para.boot[k,2] <- t.pis$pi.new.80[2]-t.pis$pi.new.80[1]
	 results$para.boot[k,3] <- ifelse(t.pis$pi.new.90[1]<=Y.star.new & t.pis$pi.new.90[2]>=Y.star.new, 1, 0)
	 results$para.boot[k,4] <- t.pis$pi.new.90[2]-t.pis$pi.new.90[1]
	 results$para.boot[k,5] <- ifelse(t.pis$pi.new.95[1]<=Y.star.new & t.pis$pi.new.95[2]>=Y.star.new, 1, 0)
	 results$para.boot[k,6] <- t.pis$pi.new.95[2]-t.pis$pi.new.95[1]
	
	 results$para.boot[k,7] <- ifelse(t.pis$pi.exs.80[1]<=Y.star.exs & t.pis$pi.exs.80[2]>=Y.star.exs, 1, 0)
	 results$para.boot[k,8] <- t.pis$pi.exs.80[2]-t.pis$pi.exs.80[1]
	 results$para.boot[k,9] <- ifelse(t.pis$pi.exs.90[1]<=Y.star.exs & t.pis$pi.exs.90[2]>=Y.star.exs, 1, 0)
	 results$para.boot[k,10] <- t.pis$pi.exs.90[2]-t.pis$pi.exs.90[1]
	 results$para.boot[k,11] <- ifelse(t.pis$pi.exs.95[1]<=Y.star.exs & t.pis$pi.exs.95[2]>=Y.star.exs, 1, 0)
	 results$para.boot[k,12] <- t.pis$pi.exs.95[2]-t.pis$pi.exs.95[1]

	

		#### Bayes - half cauchy priors for variance components
  
		pr1 <-  brms::prior(cauchy(0, 1), class = "sd") + 
 		brms::prior(cauchy(0, 1), class = "sigma")
 		bayes.cauchy.posterior <- brm(response ~ 1 + (1 | group), data = my.data, prior = pr1)
  
		bayes.cauchy.theta <- posterior_epred(bayes.cauchy.posterior, newdata = data.frame(group = I+1), re.form = NULL,allow_new_levels = TRUE)	
 		bayes.cauchy.pdi.theta.80 <- quantile(bayes.cauchy.theta, probs = c(0.1, 0.8))
 		bayes.cauchy.pdi.theta.90 <- quantile(bayes.cauchy.theta, probs = c(0.05, 0.95))
		bayes.cauchy.pdi.theta.95 <- quantile(bayes.cauchy.theta, probs = c(0.025, 0.975))
  
 		bayes.cauchy.new <- posterior_predict(bayes.cauchy.posterior, newdata = data.frame(group = I+1), re.form = NULL,allow_new_levels = TRUE)	
 		bayes.cauchy.pdi.new.80 <- quantile(bayes.cauchy.new, probs = c(0.1, 0.8))
 		bayes.cauchy.pdi.new.90 <- quantile(bayes.cauchy.new, probs = c(0.05, 0.95))
 		bayes.cauchy.pdi.new.95 <- quantile(bayes.cauchy.new, probs = c(0.025, 0.975))
  
 		bayes.cauchy.exs <- posterior_predict(bayes.cauchy.posterior, newdata = data.frame(group = I), re.form = NULL)	
		bayes.cauchy.pdi.exs.80 <- quantile(bayes.cauchy.exs, probs = c(0.1, 0.8))
 		bayes.cauchy.pdi.exs.90 <- quantile(bayes.cauchy.exs, probs = c(0.05, 0.95))
 		bayes.cauchy.pdi.exs.95 <- quantile(bayes.cauchy.exs, probs = c(0.025, 0.975))
  
 		results$bayes.cauchy[k,1] <- ifelse(bayes.cauchy.pdi.new.80[1]<=Y.star.new & bayes.cauchy.pdi.new.80[2]>=Y.star.new, 1, 0)
 		results$bayes.cauchy[k,2] <- bayes.cauchy.pdi.new.80[2]-bayes.cauchy.pdi.new.80[1]
 		results$bayes.cauchy[k,3] <- ifelse(bayes.cauchy.pdi.new.90[1]<=Y.star.new & bayes.cauchy.pdi.new.90[2]>=Y.star.new, 1, 0)
 		results$bayes.cauchy[k,4] <- bayes.cauchy.pdi.new.90[2]-bayes.cauchy.pdi.new.90[1]
 		results$bayes.cauchy[k,5] <- ifelse(bayes.cauchy.pdi.new.95[1]<=Y.star.new & bayes.cauchy.pdi.new.95[2]>=Y.star.new, 1, 0)
 		results$bayes.cauchy[k,6] <- bayes.cauchy.pdi.new.95[2]-bayes.cauchy.pdi.new.95[1]
  
 		results$bayes.cauchy[k,7] <- ifelse(bayes.cauchy.pdi.exs.80[1]<=Y.star.exs & bayes.cauchy.pdi.exs.80[2]>=Y.star.exs, 1, 0)
 		results$bayes.cauchy[k,8] <- bayes.cauchy.pdi.exs.80[2]-bayes.cauchy.pdi.exs.80[1]
		results$bayes.cauchy[k,9] <- ifelse(bayes.cauchy.pdi.exs.90[1]<=Y.star.exs & bayes.cauchy.pdi.exs.90[2]>=Y.star.exs, 1, 0)
 		results$bayes.cauchy[k,10] <- bayes.cauchy.pdi.exs.90[2]-bayes.cauchy.pdi.exs.90[1]
 		results$bayes.cauchy[k,11] <- ifelse(bayes.cauchy.pdi.exs.95[1]<=Y.star.exs & bayes.cauchy.pdi.exs.95[2]>=Y.star.exs, 1, 0)
 		results$bayes.cauchy[k,12] <- bayes.cauchy.pdi.exs.95[2]-bayes.cauchy.pdi.exs.95[1]
  
 		results$bayes.cauchy[k,13] <- ifelse(bayes.cauchy.pdi.theta.80[1]<=theta.star & bayes.cauchy.pdi.theta.80[2]>=theta.star, 1, 0)
 		results$bayes.cauchy[k,14] <- bayes.cauchy.pdi.theta.80[2]-bayes.cauchy.pdi.theta.80[1]
		results$bayes.cauchy[k,15] <- ifelse(bayes.cauchy.pdi.theta.90[1]<=theta.star & bayes.cauchy.pdi.theta.90[2]>=theta.star, 1, 0)
 		results$bayes.cauchy[k,16] <- bayes.cauchy.pdi.theta.90[2]-bayes.cauchy.pdi.theta.90[1]
 		results$bayes.cauchy[k,17] <- ifelse(bayes.cauchy.pdi.theta.95[1]<=theta.star & bayes.cauchy.pdi.theta.95[2]>=theta.star, 1, 0)
 		results$bayes.cauchy[k,18] <- bayes.cauchy.pdi.theta.95[2]-bayes.cauchy.pdi.theta.95[1]
 
  #### Non-parametric bootstrap
  
  
  boot.response.predictions <- rep(NA, 1000)
  for(b in 1:1000){
    boot.responses <- sample(response, n, replace = TRUE)
    boot.response.predictions[b] <- sample(boot.responses, 1)
  }
  
  boot.pdi.new.80 <- quantile(boot.response.predictions, c(0.1, 0.9))
  boot.pdi.new.90 <- quantile(boot.response.predictions, c(0.05, 0.95))
  boot.pdi.new.95 <- quantile(boot.response.predictions, c(0.025, 0.975))
  
  results$boot[k,1] <- ifelse(boot.pdi.new.80[1]<=Y.star.new & boot.pdi.new.80[2]>=Y.star.new, 1, 0)
  results$boot[k,2] <- boot.pdi.new.80[2]-boot.pdi.new.80[1]
  results$boot[k,3] <- ifelse(boot.pdi.new.90[1]<=Y.star.new & boot.pdi.new.90[2]>=Y.star.new, 1, 0)
  results$boot[k,4] <- boot.pdi.new.90[2]-boot.pdi.new.90[1]
  results$boot[k,5] <- ifelse(boot.pdi.new.95[1]<=Y.star.new & boot.pdi.new.95[2]>=Y.star.new, 1, 0)
  results$boot[k,6] <- boot.pdi.new.95[2]-boot.pdi.new.95[1]
  
  boot.response.last.group.predictions <- rep(NA, 1000)
  for(b in 1:1000){
    boot.group.responses <- sample(response.last.group, des[I], replace = TRUE)
    boot.response.last.group.predictions[b] <- sample(boot.group.responses, 1)
  }
  
  boot.pdi.exs.80 <- quantile(boot.response.last.group.predictions, c(0.1, 0.9))
  boot.pdi.exs.90 <- quantile(boot.response.last.group.predictions, c(0.05, 0.95))
  boot.pdi.exs.95 <- quantile(boot.response.last.group.predictions, c(0.025, 0.975))
  
  results$boot[k,7] <- ifelse(boot.pdi.exs.80[1]<=Y.star.exs & boot.pdi.exs.80[2]>=Y.star.exs, 1, 0)
  results$boot[k,8] <- boot.pdi.exs.80[2]-boot.pdi.exs.80[1]
  results$boot[k,9] <- ifelse(boot.pdi.exs.90[1]<=Y.star.exs & boot.pdi.exs.90[2]>=Y.star.exs, 1, 0)
  results$boot[k,10] <- boot.pdi.exs.90[2]-boot.pdi.exs.90[1]
  results$boot[k,11] <- ifelse(boot.pdi.exs.95[1]<=Y.star.exs & boot.pdi.exs.95[2]>=Y.star.exs, 1, 0)
  results$boot[k,12] <- boot.pdi.exs.95[2]-boot.pdi.exs.95[1]
  
  my.data = data.frame(response, group)
  pdi.boot <- pred_int_boot(formula = response ~ 1 + (1|group), data = my.data, level = c(0.8,0.9,0.95), R = 5000)
  boot.pdi.theta.80 <- c(pdi.boot[2], pdi.boot[5])
  boot.pdi.theta.90 <- c(pdi.boot[3], pdi.boot[6])
  boot.pdi.theta.95 <- c(pdi.boot[4], pdi.boot[7])
  
  results$boot[k,13] <- ifelse(boot.pdi.theta.80[1]<=theta.star & boot.pdi.theta.80[2]>=theta.star, 1, 0)
  results$boot[k,14] <- boot.pdi.theta.80[2]-boot.pdi.theta.80[1]
  results$boot[k,15] <- ifelse(boot.pdi.theta.90[1]<=theta.star & boot.pdi.theta.90[2]>=theta.star, 1, 0)
  results$boot[k,16] <- boot.pdi.theta.90[2]-boot.pdi.theta.90[1]
  results$boot[k,17] <- ifelse(boot.pdi.theta.95[1]<=theta.star & boot.pdi.theta.95[2]>=theta.star, 1, 0)
  results$boot[k,18] <- boot.pdi.theta.95[2]-boot.pdi.theta.95[1]
  
  



	oof <- proc.time() - then
	print(c(k, oof[1]))
}


colMeans(results$im)


