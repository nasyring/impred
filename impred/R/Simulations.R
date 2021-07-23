### STEP 1. load functions from the IM.R script

source("IM.R")

### STEP 2. simulate data sets

set.seed(54321) # for reproducibility, however, the use of c++ means we cannot perfectly reproduce runs because we cannot set the seed used in MCMC sampling inside c++
K <- 2000
Z.1 <- make.Z(c(13,13,13,13,13,13,13,13,13,14,1)) # this creates a design matrix for the model with the associated group structure.  This particular setup has 10 groups of 13 observations, plus two new observations we want to predict, one in an existing group and one in a new group
n.1 <- sum(c(13,13,13,13,13,13,13,13,13,14,1))
n <- n.1-2  # 130 total data points, and two unobserved points we will predict
Z <- Z.1[1:n,1:(ncol(Z.1)-1)]	# the design matrix taking out the new observations
sigma_a2 <- 0.1	# within group standard deviation
sigma_2 <- 2	# across group standard deviation
mu <- 6	# mean parameter
Y.mat <- matrix(0,n.1,K)	# generating matrix of data plus two new observations
for(i in 1:K){
	Y.mat[,i] <- sample.Y(Z.1, sigma_a2, sigma_2, mu)
}
Y <- Y.mat[1:n,]	# keeping just the n observations
Y.new <- t(Y.mat[(n+1):n.1,])	# picking out just the 2 new obs for prediction


### STEP 3. compute fused plausibility curves for each data set

plaus.results <- matrix(0,500,2*K)	# we have a total of K data sets of size n each along with 2 new obs to predict for each of the K data sets

# looping through all K sets and computing IM predictive plausibility curves for our two predictions: one within group and one new group
for(k in 1:K){
	time <- proc.time()
	ex.data <- list(Y=Y[,k], Z=Z)
	ex.statistics <- aov.statistics(ex.data)
	ex.grid <- matrix(seq(from = min(ex.data$Y)-6, to = max(ex.data$Y)+6, length.out = 500),500,1)
	sig.grid = as.matrix(expand.grid(seq(from = 0.01, to = 7, length.out = 14), seq(from = 0.01, to = 7, length.out = 14)))
	ex.grid.plauses <-  apply.plaus(sig.grid, matrix(runif(10000),10000,1), 1, ex.data, ex.statistics, ex.grid)
	plaus.results[,(2*k-1)] <- ex.grid.plauses$fused.plaus.w
	plaus.results[,(2*k)] <- ex.grid.plauses$fused.plaus.n
	time2 <- proc.time() - time
	print(c(k, time2))
}


### STEP 4. check coverage of prediction intervals

# looping through all K data sets and computing prediction intervals at 95%, 90%, and 80% from our IMs for our 2 predictions each

coverage<-matrix(0,K,6)
length<-coverage
for(k in 1:K){
	# within prediction
	plaus.max <- which.max(plaus.results[,(2*k-1)])
	lower.plaus <- which.min(abs(.05-plaus.results[1:plaus.max,(2*k-1)]))
	upper.plaus <- which.min(abs(.05-plaus.results[(plaus.max+1):500,(2*k-1)])) +plaus.max 
	pi.95 <- c(ex.grid[lower.plaus], ex.grid[upper.plaus])
	coverage[k,1] <- ifelse((pi.95[1]<=Y.new[k,1]) & (pi.95[2]>=Y.new[k,1]), 1, 0)
	length[k,1] <- pi.95[2]-pi.95[1]
	lower.plaus <- which.min(abs(.1-plaus.results[1:plaus.max,(2*k-1)]))
	upper.plaus <- which.min(abs(.1-plaus.results[(plaus.max+1):500,(2*k-1)])) +plaus.max 
	pi.90 <- c(ex.grid[lower.plaus], ex.grid[upper.plaus])
	coverage[k,2] <- ifelse((pi.90[1]<=Y.new[k,1]) & (pi.90[2]>=Y.new[k,1]), 1, 0)
	length[k,2] <- pi.90[2]-pi.90[1]	
	lower.plaus <- which.min(abs(.2-plaus.results[1:plaus.max,(2*k-1)]))
	upper.plaus <- which.min(abs(.2-plaus.results[(plaus.max+1):500,(2*k-1)]))+plaus.max 
	pi.80 <- c(ex.grid[lower.plaus], ex.grid[upper.plaus])
	coverage[k,3] <- ifelse((pi.80[1]<=Y.new[k,1]) & (pi.80[2]>=Y.new[k,1]), 1, 0)
	length[k,3] <- pi.80[2]-pi.80[1]

	# new prediction
	plaus.max <- which.max(plaus.results[,2*k])
	lower.plaus <- which.min(abs(.05-plaus.results[1:plaus.max,2*k]))
	upper.plaus <- which.min(abs(.05-plaus.results[(plaus.max+1):500,2*k])) +plaus.max 
	pi.95 <- c(ex.grid[lower.plaus], ex.grid[upper.plaus])
	coverage[k,4] <- ifelse((pi.95[1]<=Y.new[k,2]) & (pi.95[2]>=Y.new[k,2]), 1, 0)
	length[k,4] <- pi.95[2]-pi.95[1]
	lower.plaus <- which.min(abs(.1-plaus.results[1:plaus.max,2*k]))
	upper.plaus <- which.min(abs(.1-plaus.results[(plaus.max+1):500,2*k])) +plaus.max 
	pi.90 <- c(ex.grid[lower.plaus], ex.grid[upper.plaus])
	coverage[k,5] <- ifelse((pi.90[1]<=Y.new[k,2]) & (pi.90[2]>=Y.new[k,2]), 1, 0)
	length[k,5] <- pi.90[2]-pi.90[1]	
	lower.plaus <- which.min(abs(.2-plaus.results[1:plaus.max,2*k]))
	upper.plaus <- which.min(abs(.2-plaus.results[(plaus.max+1):500,2*k]))+plaus.max 
	pi.80 <- c(ex.grid[lower.plaus], ex.grid[upper.plaus])
	coverage[k,6] <- ifelse((pi.80[1]<=Y.new[k,2]) & (pi.80[2]>=Y.new[k,2]), 1, 0)
	length[k,6] <- pi.80[2]-pi.80[1]
}


#### STEP 5: Compare to bootstrap - using Fernando's bootstrap function in predintma

group.des <- c(13,13,13,13,13,13,13,13,13,13)
k <- length(group.des)
group <- c()
for(j in 1:k){
	group <- c(group, rep(j,group.des[j])) 
}
coverage.boot.w <- matrix(NA, K,3)
length.boot.w <- matrix(NA, K,3)
coverage.boot.n <- matrix(NA, K,3)
length.boot.n <- matrix(NA, K,3)


# first set of bootstrap prediction intervals is for the new observation in an existing group
for(i in 1:K){
	response <- c(Y[,i],Y.new[i,1])
	my.data = data.frame(response, group)
	pdi.80 <- pred_int_boot(formula = response ~ 1 + (1|group), data = my.data, level = 0.8, R = 5000)
	pdi.90 <- pred_int_boot(formula = response ~ 1 + (1|group), data = my.data, level = 0.9, R = 5000)
	pdi.95 <- pred_int_boot(formula = response ~ 1 + (1|group), data = my.data, level = 0.95, R = 5000)
	coverage.boot.w[i,3] <- ifelse((pdi.80[2]<=Y.new[i,2]) & (pdi.80[3]>=Y.new[i,2]), 1, 0)
	length.boot.w[i,3] <- pdi.80[3]-pdi.80[2]
	coverage.boot.w[i,2] <- ifelse((pdi.90[2]<=Y.new[i,2]) & (pdi.90[3]>=Y.new[i,2]), 1, 0)
	length.boot.w[i,2] <- pdi.90[3]-pdi.90[2]
	coverage.boot.w[i,1] <- ifelse((pdi.95[2]<=Y.new[i,2]) & (pdi.95[3]>=Y.new[i,2]), 1, 0)
	length.boot.w[i,1] <- pdi.95[3]-pdi.95[2]
	if(i%%100 == 0) print(i)
}

colMeans(coverage.boot.w, na.rm = TRUE)
colMeans(length.boot.w)

> colMeans(coverage.boot[1:500,])
[1] 0.552 0.484 0.372
> colMeans(length.boot[1:500,])
[1] 2.918202 2.467598 1.929513

> colMeans(coverage.boot, na.rm = TRUE)
[1] 0.5262631 0.4532266 0.3621811
> colMeans(length.boot, na.rm = TRUE)
[1] 2.926005 2.473079 1.934625



# second set of bootstrap prediction intervals is for the new observation in an existing group
for(i in 1:K){
	response <- c(Y[,i],Y.new[i,2])
	my.data = data.frame(response, group)
	pdi.80 <- pred_int_boot(formula = response ~ 1 + (1|group), data = my.data, level = 0.8, R = 5000)
	pdi.90 <- pred_int_boot(formula = response ~ 1 + (1|group), data = my.data, level = 0.9, R = 5000)
	pdi.95 <- pred_int_boot(formula = response ~ 1 + (1|group), data = my.data, level = 0.95, R = 5000)
	coverage.boot.n[i,3] <- ifelse((pdi.80[2]<=Y.new[i,2]) & (pdi.80[3]>=Y.new[i,2]), 1, 0)
	length.boot.n[i,3] <- pdi.80[3]-pdi.80[2]
	coverage.boot.n[i,2] <- ifelse((pdi.90[2]<=Y.new[i,2]) & (pdi.90[3]>=Y.new[i,2]), 1, 0)
	length.boot.n[i,2] <- pdi.90[3]-pdi.90[2]
	coverage.boot.n[i,1] <- ifelse((pdi.95[2]<=Y.new[i,2]) & (pdi.95[3]>=Y.new[i,2]), 1, 0)
	length.boot.n[i,1] <- pdi.95[3]-pdi.95[2]
	if(i%%100 == 0) print(i)
}

colMeans(coverage.boot.n, na.rm = TRUE)
colMeans(length.boot.n)
