###############################################################################
#### Step 1: Load packages, functions and data
###############################################################################

#library(devtools)
#install_github("nasyring/impred", subdir = "impred")
library(impred)
library(lme4)
library(predintma)
library(tidyverse)
library(matlib)
library(ggplot2)
library(brms)
library(rstanarm)

# Load IM.R from github/nasyring/impred/CodesForSimulations
#source("IM.R")

# Load data from github/nasyring/impred/CodesForSimulations
# load('G:\\prediction\\soysff.rda')
# load('C:\\Users\\nsyring\\Desktop\\prediction\\soysff.rda')
names(soysff)[names(soysff) == "Trial_ID"] <- "groups"


###############################################################################
#### Step 2: Fit model using lmer
###############################################################################

# fit random effects anova
model <- fit.rand.fx <- lmer(soysff$lrr ~ (1 | groups), data = soysff, REML=TRUE)
#boundary (singular) fit: see ?isSingular
groups <- soysff$groups

# alternatively compute MLEs
neg.log.lik <- function(par, y, Z){
	n<-length(y)
	mu <- par[1]; sigma2_a <- exp(par[2]); sigma2 <- exp(par[3]) # exponentiate variance components so they aren't <=0
	Sigma <- sigma2_a*(Z%*%t(Z))+diag(sigma2, n)
	Sigma.inv <- solve(Sigma)
	eigen.S <- eigen(Sigma)
	M<-matrix(y-rep(mu,n), n, 1)
	return(sum(log(eigen.S$values))+(t(M)%*%Sigma.inv%*%M))
}

MLE <- function(param, y, Z){
	return(optim(param, neg.log.lik, y=y, Z=Z))
} 

tab.Z<-as.numeric(table(soysff$groups))
n_i <- tab.Z
I<-length(tab.Z)
n<-sum(tab.Z)
Z <- make.Z(tab.Z)
Y <- soysff$lrr
mu <- mean(Y)
sigma2_a <- sigma2 <- log(var(Y))
param <- c(mu, sigma2_a, sigma2)
mles<-MLE(param, Y,Z)
mles
#$par
#[1]  0.02810272 -21.35012751  -5.53543359
# estimates are 0.028, exp(-21.35012751 ) = 5.34e-10 and exp( -5.53543359) = 3.94e-03
#


#############################################################################
#### Step 3: Plot data
#############################################################################

# Plot the data 
bxpltDF <- data.frame(Y = Y, Group = soysff$groups) %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(median.group = median(Y)) %>%
  dplyr::ungroup()
bxpltDF <- bxpltDF %>% dplyr::mutate(grp.rnk = min_rank(median.group))

boxplot(Y~grp.rnk, data = bxpltDF, horizontal = TRUE, axes = FALSE, ylab = 'Farm', xlab = 'Yield', ylim = c(-0.3,0.3))
axis(1, at=c(-0.3, -0.2,-0.1,0,0.1,0.2,0.3))

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 14))

# ggplot version
ggplot(data = bxpltDF, mapping = aes(group = grp.rnk, y = Y)) + 
  geom_boxplot() +
  ylab("Fungicide Effect") +
  theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()) +
  coord_flip()+My_Theme

tsum <- function(x, var.names = c("y","trial"), 
                 method = c("ici","lm","min-max"),
                 level = 0.95, order = TRUE,
                 add.id = TRUE,
                 formula = NULL){
  
  if(class(x) != "data.frame") stop("only for data.frames")
  
  method <- match.arg(method)
  if(!missing(formula)) var.names <- all.vars(as.formula(formula))
  
  if(method == "ici"){
    ## Calcualte mean by trial
    frm <- as.formula(paste0(var.names[1],"~",var.names[2]))
    ans <- aggregate(formula = frm, data = x, FUN = mean)
    ans$lb <- aggregate(formula = frm, data = x, FUN = lb_fun, level = level)[,2]
    ans$ub <- aggregate(formula = frm, data = x, FUN = ub_fun, level = level)[,2]
    ans$n <- aggregate(formula = frm, data = x, FUN = length)[,2]
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(method == "lm"){
    ## Calcualte intervals using lm
    frm <- as.formula(paste0(var.names[1],"~",var.names[2]))
    nt <- aggregate(formula = frm, data = x, FUN = length)[,2]
    fit <- lm(formula = frm, data = x)
    ndat <- data.frame(unique(x[,var.names[2]]))
    names(ndat) <- var.names[2]
    prd <- predict(fit, newdata = ndat, interval = "confidence")
    ans <- cbind(ndat, prd)
    ans$n <- nt
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(method == "min-max"){
    ## Calcualte mean by trial
    frm <- as.formula(paste0(var.names[1],"~",var.names[2]))
    ans <- aggregate(formula = frm, data = x, FUN = mean)
    ans$lb <- aggregate(formula = frm, data = x, FUN = min)[,2]
    ans$ub <- aggregate(formula = frm, data = x, FUN = max)[,2]
    ans$n <- aggregate(formula = frm, data = x, FUN = length)[,2]
    names(ans) <- c(var.names[2], "m","lb","ub","n")
  }
  
  if(order & add.id){
    ans2 <- ans[order(ans[,"m"]),]
    ans2$id <- 1:nrow(ans2)
    return(ans2)
  }else{
    return(ans)
  }
}

lb_fun <- function(x, level = 0.95){
  n.x <- length(x)
  m.x <- mean(x)
  se.x <- sd(x)/sqrt(n.x)
  qnt <- 1 - (1 - level)/2
  tv <- qt(qnt, n.x-1)
  lbv <- m.x -  tv * se.x
  return(lbv)
}

ub_fun <- function(x, level = 0.95){
  n.x <- length(x)
  m.x <- mean(x)
  se.x <- sd(x)/sqrt(n.x)
  qnt <- 1 - (1 - level)/2
  tv <- qt(qnt, n.x-1)
  ubv <- m.x +  tv * se.x
  return(ubv)
}

includes <- function(x, lb, ub){
  
  if(ub < lb) stop("upper bound should be greater than lower bound")
  tmp <- 0
  for(i in 1:length(x)){
    if(x[i] >= lb & x[i] <= ub){
      tmp <- tmp + 1
    }
  }
  ans <- tmp/length(x)
  ans
}

data("soysff")

ggplot(data = soysff, aes(x = lrr, y = Trial_ID)) + 
  geom_point() + 
  ggtitle("Raw data")

soysff.tsum <- tsum(soysff, var.names = c("lrr", "Trial_ID"))

ggplot(soysff.tsum, aes(x = m, y = id)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = lb, xmax = ub)) +
  geom_vline(xintercept = 0) + 
  ggtitle("Response to fungicide in soybeans")

## New graph
tmp.soysff <- soysff.tsum[, c("Trial_ID", "id", "m")]
tmp.soysff$Trial_ID <- factor(tmp.soysff$Trial_ID, levels = tmp.soysff$Trial_ID)
nsoysff <- merge(tmp.soysff, soysff)

ggplot(data = nsoysff, aes(x = lrr, y = Trial_ID)) + 
  geom_point(alpha = 0.4) + 
  geom_point(aes(x = m), size = 1.1, color = "blue") + 
  geom_vline(xintercept = 0) + 
  ylab("Trial ID") + xlab("log RR") + 
  theme_classic()







############################################################################################
#### Step 4: IM method
############################################################################################


ex.data <- list(Y=Y, Z=Z)
ex.statistics <- aov.statistics(ex.data)
ex.grid.size <- 300
ex.grid <- matrix(seq(from = min(ex.data$Y), to = max(ex.data$Y), length.out = ex.grid.size),ex.grid.size,1)

fun <- function(model){
	summ <- summary(model)
	return((summ$sigma^2)/((summ$sigma^2)+summ$varcor[1]$groups[1]))
}
booted.ratio <- bootMer(model, fun, 100)
boot.mean <- mean(booted.ratio$t)
boot.sd <- sd(booted.ratio$t)
ratio.grid <- seq(from = max(boot.mean-boot.sd, 0.01), to = min(boot.mean+boot.sd, 0.99), length.out = 3)
plauses <- plaus_aov_unbalanced_full(ex.data, ex.grid, (1:99)/100)
plauses2 <- plaus_aov_unbalanced_full(ex.data, ex.grid, ratio.grid)

plot(ex.grid, plauses2$plauses.theta, type = 'l', ylim = c(0,1), xlim = c(-0.05, 0.1), xlab = expression(theta), ylab = expression(pi[n](theta)))
min.plaus <- apply(plauses2$plauses.all,1,min)
polygon(c(ex.grid, rev(ex.grid)), c(min.plaus, rev(plauses2$plauses.theta)),col = 'grey', lty = 0)

lines(ex.grid, plauses$plauses.theta, lwd = 2, lty = 2)
lines(ex.grid, plauses2$plauses.theta, lwd = 2)

plaus.max <- which.max(plauses$plauses.theta)
lower1 <- ex.grid[which.min(abs(0.05 - plauses$plauses.theta[1:plaus.max]))]
upper1 <- ex.grid[plaus.max+which.min(abs(0.05 - plauses$plauses.theta[(plaus.max+1):length(ex.grid)]))]

c(lower1, upper1)
> c(lower1, upper1)
[1] -0.1037250  0.1598329


plaus.max <- which.max(plauses2$plauses.theta)
lower2 <- ex.grid[which.min(abs(0.05 - plauses2$plauses.theta[1:plaus.max]))]
upper2 <- ex.grid[plaus.max+which.min(abs(0.05 - plauses2$plauses.theta[(plaus.max+1):length(ex.grid)]))]

c(lower2, upper2)
> c(lower2, upper2)
[1] -0.002356612  0.058464451

points(c(lower2, upper2), c(0.05, 0.05), pch = 19)
lines(c(lower2,0.05), c(0.05, 0.4), lty = 2 )
lines(c(upper2,0.05), c(0.05, 0.4), lty = 2 )
text(0.075, 0.41, "95% Prediction Interval")
text(0.075, 0.36, "(-0.0024  0.0585)")

points(-0.0006671377, 0.0578, pch = 19)
lines(c(-0.0006671377,-0.01), c(0.0578,0.3), lty = 2 )
text(-0.024, 0.32, expression(pi[n](0) %~~% 0.0578))

plaus.max <- which.max(plauses2$plauses.new)
lower <- ex.grid[which.min(abs(0.05 - plauses2$plauses.new[1:plaus.max]))]
upper <- ex.grid[plaus.max+which.min(abs(0.05 - plauses2$plauses.new[(plaus.max+1):length(ex.grid)]))]
c(lower, upper)
#> c(lower, upper)
#[1] -0.09696715  0.15307499


#####################################################################################################
#### Step 5: Student t intervals
#####################################################################################################


stdt.pi <- mles$par[1] + sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(n_i^2))+exp(mles$par[3])/n)*qt(c(0.025, 0.975), df = I-1)
stdt.pi
#[1] 0.01909582 0.03710961

stdt.pi.new <- mles$par[1] + sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(n_i^2))+exp(mles$par[3])*(1+1/n))*qt(c(0.025, 0.975), df = I-1)
stdt.pi.new
#[1] -0.09959026  0.15579569

#####################################################################################################
#### Step 6: Non-parametric bootstrap
#####################################################################################################

pdi.boot <- pred_int_boot(formula = Y ~ 1 + (1|Group), data = bxpltDF, level = 0.95, R = 10000)
pdi.95 <- c(pdi.boot[2], pdi.boot[3])
pdi.95
#       2.5%       97.5% 
#-0.03705742  0.11904215 

B = 10000
quants.boot <- matrix(NA, B,2)
for(k in 1:B){
	data.boot <- rep(NA,n)
	start <- 1
	end <- n_i[1]
	for(i in 1:I){
		if(i!=1){
			 start <- end+1
			 end <- start-1+n_i[i]
		}
		index <- sample.int(n_i[i], n_i[i], replace = TRUE)
		Y.n_i<-Y[start:end]
		data.boot[start:end]<-Y.n_i[index]
	}
	quants.boot[k,] <- quantile(data.boot, c(0.025, 0.975))
if(k%%100==0) print(k)
}
colMeans(quants.boot)
#> colMeans(quants.boot)
#[1] -0.08203826  0.16248958




#####################################################################################################
#### Step 7: Conformal prediction
#####################################################################################################


pred.lvl.theta <- (1/(I+1))*floor(0.95*(I+1))
pred.alpha.theta <- (1-pred.lvl.theta)

pred.lvl.new <- (1/(n+1))*floor(0.95*(n+1))
pred.alpha.new <- (1-pred.lvl.new)

nonconformity <- function(y, x) abs(mean(y)-x)
nonconformity.i <- function(i, y){
  x<-y[i]
  y<-y[-i]
  return(nonconformity(y,x))
}
nonconformity.all <- function(y){
  n<-length(y)
  return(apply(matrix(1:n,n,1),1,nonconformity.i, y=y))
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

  response.group.averages <- (t(Z)%*%matrix(Y,n,1))/n_i
  mean.averages <- mean(response.group.averages)

grid.theta.size <- 480
  grid.theta <- seq(from = mean.averages	 - 4, to = mean.averages + 4, length.out = grid.theta.size)
  conformal.plaus.theta <- nonconformity.plaus.grid(grid.theta, response.group.averages)
  conf.plaus.max.theta <- which.max(conformal.plaus.theta)	
  conf.pdi.theta <- c(grid.theta[which.min(abs(pred.alpha.theta-conformal.plaus.theta[1:conf.plaus.max.theta]))], grid.theta[conf.plaus.max.theta+which.min(abs(pred.alpha.theta-conformal.plaus.theta[(conf.plaus.max.theta+1):grid.theta.size]))])
  conf.pdi.theta
#[1] -0.04796395  0.10234920


  grid.new.size <- 480
  grid.new <- seq(from = mean(Y)	 - 4, to = mean(Y) + 4, length.out = grid.new.size)
  conformal.plaus.new <- nonconformity.plaus.grid(grid.new, Y)
  conf.plaus.max.new <- which.max(conformal.plaus.new)	 
  conf.pdi.new <- c(grid.new[which.min(abs(pred.alpha.new-conformal.plaus.new[1:conf.plaus.max.new]))], grid.new[conf.plaus.max.new+which.min(abs(pred.alpha.new-conformal.plaus.new[(conf.plaus.max.new+1):grid.new.size]))])
#> conf.pdi.new
#[1] -0.09715695  0.15336497


#####################################################################################################
#### Step 7: Parametric Bootstrap
#####################################################################################################

mu <- mean(Y)
sigma2_a <- sigma2 <- log(var(Y))
param <- c(mu, sigma2_a, sigma2)


B = 1000
theta.star.boot <- rep(NA, B)
Y.star.boot <- rep(NA, B)
estimators <- c(mles$par[1], exp(mles$par[2]), exp(mles$par[3]))
for(k in 1:B){
	data.boot <- rep(NA,n)
	start <- 1
	end <- n_i[1]
	for(i in 1:I){
		if(i!=1){
			 start <- end+1
			 end <- start-1+n_i[i]
		}
		Y.n_i<-rnorm(n_i[i], estimators[1], sqrt(estimators[2]+estimators[3]))
		data.boot[start:end]<-Y.n_i
	}
	mles.boot <- MLE(param, data.boot, Z)
	theta.star.boot[k] <- rnorm(1, mles.boot$par[1], sqrt(exp(mles.boot$par[2])))
	Y.star.boot[k] <- rnorm(1, mles.boot$par[1], sqrt(exp(mles.boot$par[2])+exp(mles.boot$par[3])))	
print(k)
}

pboot.pi<-quantile(theta.star.boot, c(0.025, 0.975))
#       2.5%       97.5% 
#0.006967937 0.049013914 
quantile(Y.star.boot, c(0.025, 0.975))
#       2.5%       97.5% 
#-0.09380917  0.14995249 


#####################################################################################################
#### Step 8: Bayesian 
#####################################################################################################


## This is the Bayes approach, default priors
fm1 <- stan_lmer(lrr ~ 1 + (1 | Trial_ID), data = soysff)

## Prediction intervals from Stan
fm1.pdi0 <- posterior_epred(fm1, newdata = data.frame(Trial_ID = "new-trial"), re.form = NULL)
fm1.pdi <- quantile(fm1.pdi0, probs = c(0.025, 0.975))
#> fm1.pdi
#      2.5%      97.5% 
#0.00203831 0.05428976 

fm1.pdin <- posterior_predict(fm1, newdata = data.frame(Trial_ID = "new-trial"), re.form = NULL)
fm1.pdi.new <- quantile(fm1.pdin, probs = c(0.025, 0.975))
#> fm1.pdi.new
#       2.5%       97.5% 
#-0.09444404  0.15180969 

fm1.pdie <- posterior_epred(fm1, newdata = data.frame(Trial_ID = "ST2015IA385A"), re.form = NULL)
fm1.pdi.exs <- quantile(fm1.pdie, probs = c(0.025, 0.975))

fm1.pdie <- posterior_predict(fm1, newdata = data.frame(Trial_ID = "ST2015IA385A"), re.form = NULL)
fm1.pdi.exs <- quantile(fm1.pdie, probs = c(0.025, 0.975))

fm1.pdie <- posterior_predict(fm1, newdata = data.frame(Trial_ID = "ST2015IA315A"), re.form = NULL)
fm1.pdi.exs <- quantile(fm1.pdie, probs = c(0.025, 0.975))


## or
fm1.pdi2 <- predintma:::pred_int_stanreg(fm1, pmethod = "ntrial")


library(brms)
pr1 <- brms::prior(normal(0,4), class = "Intercept")+
brms::prior(cauchy(0,1), class = "sd")+
brms::prior(cauchy(0,1), class = "sigma")
bayes.default.posterior <- brm(lrr ~ 1 + (1 | Trial_ID), data = soysff, control = list(adapt_delta = 0.95))
bayes.default.theta <- posterior_epred(bayes.default.posterior, newdata = data.frame(Trial_ID = "new-trial"), re.form = NULL,allow_new_levels = TRUE)
bayes.default.new <- posterior_predict(bayes.default.posterior, newdata =data.frame(Trial_ID = "new-trial"), re.form = NULL,allow_new_levels = TRUE)
bayes.pdis.theta <- quantile(bayes.default.theta, c(0.025,0.975))
bayes.pdis.new <- quantile(bayes.default.new, c(0.025,0.975))
bayes.pdis.theta
bayes.pdis.new

> bayes.pdis.theta
       2.5%       97.5% 
0.003970944 0.052723651 
> bayes.pdis.new
       2.5%       97.5% 
-0.09937754  0.15161449 
> 

##########################################################################################
#### Step 9: Plotting data and all prediction intervals
##########################################################################################


## Frequentist approach
fm0 <- lmer(lrr ~ 1 + (1 | Trial_ID), data = soysff)
## ?isSingular means that the estimated group-level variance is zero

ggplot(nsoysff, aes(x = lrr, y = Trial_ID)) + 
       geom_point(alpha = 0.4) + 
       geom_point(aes(x = m), size = 1.1, color = "blue") + 
       geom_vline(xintercept = 0) + 
       geom_point(aes(x = fixef(fm0), y = -1.5), color = "orange") + 
       geom_errorbarh(aes(xmin = 0.019, xmax = 0.037, y = -1.5), color = "orange") +
       geom_text(aes(x = -0.21, y = -1.5), label = "Student's t", color = "orange") +
       geom_point(aes(x = fixef(fm0), y = -2.5), color = "blue") + 
       geom_errorbarh(aes(xmin = 0.007, xmax = 0.049, y = -2.5), color = "blue") +
       geom_text(aes(x = -0.175, y = -2.5), label = "Parametric bootstrap", color = "blue") +
       geom_point(aes(x = fixef(fm0), y = -3.5), color = "purple") + 
       geom_errorbarh(aes(xmin = 0.002, xmax = 0.054, y = -3.5), color = "purple") +
       geom_text(aes(x = -0.215, y = -3.5), label = "Bayesian", color = "purple") + 
       geom_point(aes(x = fixef(fm0), y = -4.5), color = "green") + 
       geom_errorbarh(aes(xmin = -0.0024, xmax = 0.0585, y = -4.5), color = "green") +
       geom_text(aes(x = -0.2375, y = -4.5), label = "IM", color = "green") + 
       geom_point(aes(x = fixef(fm0), y = -5.5), color = "darkgreen") +  
       geom_errorbarh(aes(xmin = -0.037, xmax = 0.119, y = -5.5), color = "darkgreen") +
       geom_text(aes(x = -0.16, y = -5.5), label = "Nonparametric bootstrap ", color = "darkgreen") +
       geom_point(aes(x = fixef(fm0), y = -6.5), color ="red") +  
       geom_errorbarh(aes(xmin = -0.048, xmax = 0.102, y = -6.5), color = "red") +
       geom_text(aes(x = -0.2075, y = -6.5), label = "Conformal ", color = "red") +
       geom_blank(aes(x = fixef(fm0), y = -7.5)) + 
       geom_label(aes(x = -0.2, y = 35), label = "observation   ") + 
	 geom_point(aes(x = -0.16, y = 35), pch = 19, col= 'darkgrey') + 
	 geom_point(aes(x = -0.16, y = 35), pch = 21) + 
       geom_label(aes(x = -0.2, y = 33), label = "trial mean   ", color = "blue") +
       geom_point(aes(x = -0.165, y = 33), color ="blue", pch = 20, cex = 2) + 
       xlab("log RR") + ylab("Trial ID") +
       theme_classic() 

