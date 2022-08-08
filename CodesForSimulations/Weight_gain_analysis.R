###############
#
#	AvgDailyGain data
#
#
###############


library(SASmixed)
library(lme4)
library(impred)
library(rstanarm)

# lmer model fit
model <- lmer(adg ~ InitWt + Treatment + (1 | Block), AvgDailyGain)
model
summ <- summary(model)
summ
reml.sige <- summ$sigma^2
reml.siga <- summ$varcor$Block[1]


# Data regression plots

plot(1,1,type = 'p', col = 'white', xlim = c(300, 500), ylim = c(0, 3), xlab = 'Initial Weight', ylab = 'Average Daily Gain')
x.seq <- seq(from = 300, to = 500, length.out = 100)
y.seq <- x.seq*0.0027797+0.2490348+0.4835076
y.seq2 <- x.seq*0.0027797+0.2490348+0.4639445
y.seq3 <- x.seq*0.0027797+0.2490348+0.5520736 
lines(x.seq, y.seq)
lines(x.seq, y.seq2)
lines(x.seq, y.seq3)
points(AvgDailyGain$InitWt[AvgDailyGain$Block==1], AvgDailyGain$adg[AvgDailyGain$Block==1], col = 'red', pch = 19)
points(AvgDailyGain$InitWt[AvgDailyGain$Block==2], AvgDailyGain$adg[AvgDailyGain$Block==2], col = 'blue', pch = 19)
points(AvgDailyGain$InitWt[AvgDailyGain$Block==3], AvgDailyGain$adg[AvgDailyGain$Block==3], col = 'green', pch = 19)
points(AvgDailyGain$InitWt[AvgDailyGain$Block==4], AvgDailyGain$adg[AvgDailyGain$Block==4], col = 'orange', pch = 19)
points(AvgDailyGain$InitWt[AvgDailyGain$Block==5], AvgDailyGain$adg[AvgDailyGain$Block==5], col = 'purple', pch = 19)
points(AvgDailyGain$InitWt[AvgDailyGain$Block==6], AvgDailyGain$adg[AvgDailyGain$Block==6], col = 'cyan', pch = 19)
points(AvgDailyGain$InitWt[AvgDailyGain$Block==7], AvgDailyGain$adg[AvgDailyGain$Block==7], col = 'darkgreen', pch = 19)
points(AvgDailyGain$InitWt[AvgDailyGain$Block==8], AvgDailyGain$adg[AvgDailyGain$Block==8], col = 'darkgrey', pch = 19)




# lmer predict point prediction
Id = 999
Block = 9
Treatment = '10'
InitWt = 400
Trt = 10
newdata <- data.frame(Id, Block, Treatment, InitWt, Trt)
predict(model, newdata, type = 'response', allow.new.levels = TRUE)


# bootMer

fun <- function(model.fit){
	summ <- summary(model.fit)
	betas <- summ$coefficients[1:3,1]
	fitvar <- summ$varcor[1]$Block[1]
	return(c(betas[1]+400*betas[2]+betas[3] + rnorm(1,0,sqrt(fitvar))))
}

my.boot <- bootMer(model, fun, nsim = 1000, type = 'parametric')

confint(my.boot,type="norm")
confint(my.boot,type="basic")
confint(my.boot,type="perc")
> confint(my.boot,type="norm")
                2.5 %   97.5 %
(Intercept) -0.497681 1.649913
> confint(my.boot,type="basic")
                 2.5 %   97.5 %
(Intercept) -0.5634584 1.645049
> 
> 
> confint(my.boot,type="perc")
                2.5 %   97.5 %
(Intercept) 0.8026439 3.011151



my.boot <- bootMer(model, fun, nsim = 1000, type = 'parametric', use.u = TRUE)

confint(my.boot,type="norm")
confint(my.boot,type="basic")
confint(my.boot,type="perc")

> my.boot <- bootMer(model, fun, nsim = 1000, type = 'parametric', use.u = TRUE)
> 
> confint(my.boot,type="norm")
               2.5 %   97.5 %
(Intercept) 1.882987 3.877679
> confint(my.boot,type="basic")
               2.5 %   97.5 %
(Intercept) 1.906208 3.902161
> confint(my.boot,type="perc")
                2.5 %   97.5 %
(Intercept) 0.8265774 2.822531

my.boot <- bootMer(model, fun, nsim = 1000, type = 'semiparametric', use.u = TRUE)

confint(my.boot,type="norm")
confint(my.boot,type="basic")
confint(my.boot,type="perc")

> my.boot <- bootMer(model, fun, nsim = 1000, type = 'semiparametric', use.u = TRUE)
> 
> confint(my.boot,type="norm")
                2.5 %  97.5 %
(Intercept) 0.9683186 2.99444
> confint(my.boot,type="basic")
                2.5 %   97.5 %
(Intercept) 0.9566139 3.027746
> confint(my.boot,type="perc")
                2.5 %   97.5 %
(Intercept) 0.7951648 2.866297
> 

# IM
X <- cbind(AvgDailyGain$InitWt, AvgDailyGain$Trt/10)
Y <- AvgDailyGain$adg
Z <- make.Z(rep(4,8))
thetaseq <- seq(from = -2, to = 6, length.out = 800)
data <- list(Y=Y,Z=Z,X=X)
x = c(400,1)
z = 1
fun <- function(model){
	summ <- summary(model)
	return(c(summ$varcor$Block[1],summ$sigma^2))
}
boot.var <- bootMer(model, fun, nsim = 100, type = "parametric")
seq.var.a <- quantile(boot.var$t[,1], seq(from = 0.05, to = 0.95, length.out = 19))
seq.var.0 <- quantile(boot.var$t[,2], seq(from = 0.05, to = 0.95, length.out = 19))
var.grid = expand.grid(seq.var.a, seq.var.0)
var.grid <- boot.var$t
plauses <- plaus_two_stage_lme(data, thetaseq, var.grid, x)

seq.var.a <- quantile(boot.var$t[,1], seq(from = 0.2, to = 0.8, length.out = 11))
seq.var.0 <- quantile(boot.var$t[,2], seq(from = 0.2, to = 0.8, length.out = 11))
var.grid = expand.grid(seq.var.a, seq.var.0)
var.grid <- boot.var$t
plauses2 <- plaus_two_stage_lme(data, thetaseq, var.grid, x)




plot(1,1, col = 'white', ylab = expression(pi[n](theta)), xlab = expression(theta), xlim = c(-1, 5), ylim = c(0,1))
lines(thetaseq, plauses$plauses.theta)
min.plaus <- apply(plauses$plauses.all.theta,1,min)
polygon(c(thetaseq, rev(thetaseq)), c(min.plaus, rev(plauses$plauses.theta)),col = 'grey', lty = 0)
lines(thetaseq, plauses$plauses.theta, lwd = 2, lty = 2)
lines(thetaseq, plauses2$plauses.theta, lwd = 2)


plaus.max <- which.max(plauses$plauses.theta)
lower <- thetaseq[which.min(abs(0.05 - plauses$plauses.theta[1:plaus.max]))]
upper <- thetaseq[plaus.max+which.min(abs(0.05 - plauses$plauses.theta[(plaus.max+1):length(thetaseq)]))]

c(lower, upper)
> c(lower, upper)
[1] -0.3078849  3.6770964


plaus.max <- which.max(plauses2$plauses.theta)
lower <- thetaseq[which.min(abs(0.05 - plauses2$plauses.theta[1:plaus.max]))]
upper <- thetaseq[plaus.max+which.min(abs(0.05 - plauses2$plauses.theta[(plaus.max+1):length(thetaseq)]))]

c(lower, upper)
> c(lower, upper)
[1] -0.2678348  3.6370463




# Bayes using stan_lmer

## This is the Bayes approach, default priors
fm1 <- stan_lmer(adg ~ InitWt + Treatment + (1 | Block), AvgDailyGain)

## Preiction intervals for Stan
fm1.pdi0 <- posterior_epred(fm1, newdata, re.form = NULL)
fm1.pdi <- quantile(fm1.pdi0, probs = c(0.025, 0.975))

fm1.pdin <- posterior_predict(fm1, newdata, re.form = NULL)
fm1.pdi.new <- quantile(fm1.pdin, probs = c(0.025, 0.975))

> fm1.pdi
     2.5%     97.5% 
0.7561323 2.9347069 
> fm1.pdi.new
     2.5%     97.5% 
0.6431642 2.9904011 










