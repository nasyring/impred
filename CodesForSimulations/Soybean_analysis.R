
# Load packages
library(lme4)
library(predintma)
library(impred)
library(tidyverse)
library(matlib)
# Load IM.R from github/nasyring/impred/CodesForSimulations
#source("IM.R")

# Load data from github/nasyring/impred/CodesForSimulations
# load('soysff.rda')

# fit random effects anova
fit.rand.fx <- lmer(soysff$lrr ~ (1 | Trial_ID), data = soysff, REML=TRUE)
#boundary (singular) fit: see ?isSingular


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

tab.Z<-as.numeric(table(soysff$Trial_ID))
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


# Plot the data 
bxpltDF <- data.frame(Y = Y, Group = soysff$Trial_ID) %>%
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



# IM method
soy.data<-list(Y=Y, Z=Z)
soy.statistics <- aov.statistics(soy.data)
soy.rand.sets <- rand.sets(sig02 = c(0.5,0.5), matrix(runif(10000),10000,1), 1, soy.data, soy.statistics)
soy.grid <- matrix(seq(from = -.15, to = .2, length.out = 1000),1000,1)
soy.plauses <- pl.0(soy.grid, soy.rand.sets)
plot.plaus(soy.plauses)

plot.im.soysff.DF <- data.frame(y = soy.plauses$plauses.T, x = soy.plauses$vals)

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 14),
axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 14))

ggplot(data = plot.im.soysff.DF, mapping = aes(x = x, y = y)) + 
  geom_line()+My_Theme


sig.grid = as.matrix(expand.grid(seq(from = 0.0001, to = 4, length.out = 20), seq(from = 0.0001, to = 4, length.out = 20)))
soy.grid.plauses <-  apply.plaus(sig.grid, matrix(runif(10000),10000,1), 3, soy.data, soy.statistics, soy.grid)

min.plauses <- apply(femiguez.grid.plauses$all.plaus.T,2,min)
fused.plaus.DF <- data.frame(x = soy.grid, y1 = soy.grid.plauses$fused.plaus.T, y2 = min.plauses)

max.plaus <- which.max(soy.grid.plauses$fused.plaus.T)
lower <- soy.grid[which.min(abs(0.05-soy.grid.plauses$fused.plaus.T[1:max.plaus]))]
lower.y <- soy.grid.plauses$fused.plaus.T[which.min(abs(0.05-soy.grid.plauses$fused.plaus.T[1:max.plaus]))]
upper <- soy.grid[max.plaus+which.min(abs(0.05-soy.grid.plauses$fused.plaus.T[(max.plaus+1):length(soy.grid)]))]
upper.y <- soy.grid.plauses$fused.plaus.T[max.plaus+which.min(abs(0.05-soy.grid.plauses$fused.plaus.T[(max.plaus+1):length(soy.grid)]))]
c(lower, upper)
#[1] -0.01161161  0.06791792

max.plaus.new <- which.max(soy.grid.plauses$fused.plaus.n)
lower.n <- soy.grid[which.min(abs(0.05-soy.grid.plauses$fused.plaus.n[1:max.plaus.new]))]
lower.y.n <- soy.grid.plauses$fused.plaus.n[which.min(abs(0.05-soy.grid.plauses$fused.plaus.n[1:max.plaus.new]))]
upper.n <- soy.grid[max.plaus.new+which.min(abs(0.05-soy.grid.plauses$fused.plaus.n[(max.plaus.new+1):length(soy.grid)]))]
upper.y.n <- soy.grid.plauses$fused.plaus.n[max.plaus.new+which.min(abs(0.05-soy.grid.plauses$fused.plaus.n[(max.plaus.new+1):length(soy.grid)]))]
c(lower.n, upper.n)
#[1] -0.05365365  0.11031031



ci.DF <- data.frame(x = c(lower, upper), y = c(lower.y, upper.y))
pl.0.DF <- data.frame(x = soy.grid[429], y=soy.grid.plauses$fused.plaus.T[429])

My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 14),
axis.text.y = element_text(size = 14),
  axis.title.y = element_text(size = 14))

ggplot(data = fused.plaus.DF, mapping = aes(x = x, y = y1)) + 
  geom_ribbon(data=fused.plaus.DF,
                  aes(x=x,ymin=y2,ymax=y1),
                   inherit.aes=FALSE,alpha=0.3,color="grey30") +
  xlim(-0.075, 0.125)+
  ylab("plausibility") +
  xlab(TeX(r"($\theta^*$)"))+
  annotate("text", x = 0.10, y = 0.3, label = "95% PI is", size = 5)+
  annotate("text", x = 0.10, y = 0.25, label = TeX(r"($(-0.012,\,0.068)$)"), size = 5)+
  geom_point(data = ci.DF,aes(x=x, y=y))+
  annotate("segment", x = lower, xend = 0.075, y = lower.y, yend = 0.235, linetype = 2) +
  annotate("segment", x = upper, xend = 0.075, y = upper.y, yend = 0.235, linetype = 2) +
  geom_point(data = pl.0.DF,aes(x=x, y=y))+
  annotate("text", x = -0.05, y = 0.3, label = "plausibility", size = 5)+
  annotate("text", x = -0.05, y = 0.25, label = TeX(r"($\theta^* \leq 0 \, is \, 0.132$)"), size = 5)+
  annotate("segment", x = -0.025, xend = 0.0, y = 0.23, yend = 0.1318248, linetype = 2) + My_Theme

################# Student t intervals


stdt.pi <- mles$par[1] + sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(n_i^2))+exp(mles$par[3])/n)*qt(c(0.025, 0.975), df = n-1)
stdt.pi

stdt.pi.dfn3 <- mles$par[1] + sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(n_i^2))+exp(mles$par[3])/n)*qt(c(0.025, 0.975), df = n-3)
stdt.pi.dfn3

stdt.pi.new <- mles$par[1] + sqrt(exp(mles$par[2])*(1+(1/(n^2))*sum(n_i^2))+exp(mles$par[3])*(1+1/n))*qt(c(0.025, 0.975), df = n-1)
stdt.pi.new
#[1] -0.09605573  0.15226116

################ nonpar bootstrap

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




################## conformal



#### Functions for conformal prediction


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


############### Parametric bootstrap intervals
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


###################### Bayesian n-trial method and new trial method 


## This is the Bayes approach, default priors
fm1 <- stan_lmer(lrr ~ 1 + (1 | Trial_ID), data = soysff)

## Preiction intervals for Stan
fm1.pdi0 <- posterior_epred(fm1, newdata = data.frame(Trial_ID = "new-trial"), re.form = NULL)
fm1.pdi <- quantile(fm1.pdi0, probs = c(0.025, 0.975))

fm1.pdin <- posterior_predict(fm1, newdata = data.frame(Trial_ID = "new-trial"), re.form = NULL)
fm1.pdi.new <- quantile(fm1.pdin, probs = c(0.025, 0.975))

fm1.pdie <- posterior_epred(fm1, newdata = data.frame(Trial_ID = "ST2015IA385A"), re.form = NULL)
fm1.pdi.exs <- quantile(fm1.pdie, probs = c(0.025, 0.975))

fm1.pdie <- posterior_predict(fm1, newdata = data.frame(Trial_ID = "ST2015IA385A"), re.form = NULL)
fm1.pdi.exs <- quantile(fm1.pdie, probs = c(0.025, 0.975))

ST2015IA315A 


fm1.pdie <- posterior_predict(fm1, newdata = data.frame(Trial_ID = "ST2015IA315A"), re.form = NULL)
fm1.pdi.exs <- quantile(fm1.pdie, probs = c(0.025, 0.975))


## or
fm1.pdi2 <- predintma:::pred_int_stanreg(fm1, pmethod = "ntrial")






##  error bars plot


soysff.tsum <- tsum(soysff, var.names = c("lrr", "Trial_ID"))

fm0 <- lmer(lrr ~ 1 + (1 | Trial_ID), data = soysff)

ggplot(soysff.tsum, aes(x = m, y = id)) + 
  geom_point() + 
  geom_errorbarh(aes(xmin = lb, xmax = ub)) +
  geom_vline(xintercept = 0) + 
  geom_point(aes(x = fixef(fm0), y = -1)) + 
  geom_errorbarh(aes(xmin = stdt.pi[1], xmax = stdt.pi[2], y = -1), lty = 1) +
  geom_point(aes(x = fixef(fm0), y = -2)) + 
  geom_errorbarh(aes(xmin = pboot.pi[1], xmax = pboot.pi[2], y = -2), lty=2) +
  geom_point(aes(x = fixef(fm0), y = -3)) + 
  geom_errorbarh(aes(xmin = fm1.pdi2[2], xmax = fm1.pdi2[3], y = -3), lty=3) +
  geom_point(aes(x = fixef(fm0), y = -4)) + 
  geom_errorbarh(aes(xmin = lower, xmax = upper, y = -4), lty=1) +
  geom_point(aes(x = fixef(fm0), y = -5)) + 
  geom_errorbarh(aes(xmin = pdi.95[1], xmax = pdi.95[2], y = -5), lty=2) +
 geom_point(aes(x = fixef(fm0), y = -6)) + 
  geom_errorbarh(aes(xmin = conf.pdi.theta[1], xmax = conf.pdi.theta[2], y = -6), lty=3) +
#  ggtitle("Response to fungicide in soybeans. Student-t interval (Red), Parametric bootsrap (Orange) \n
#          Bayesian New Trial (Blue), Nonparametric Bootstrap (Purple)")
ggsave("./figs/Soysff-PDI.png", width = 10)










