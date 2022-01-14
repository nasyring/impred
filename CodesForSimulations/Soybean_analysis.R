
# Load packages
library(lme4)
library(predintma)
library(impred)
library(tidyverse)
library(matlib)
library(latex2exp)

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

min.plauses <- apply(soy.grid.plauses$all.plaus.T,2,min)
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
