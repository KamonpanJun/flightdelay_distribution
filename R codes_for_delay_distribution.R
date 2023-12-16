####created by Miss Kamonpan Junlat####
####This file for analysis delay distribution###
### central region #####
### 2017 ####
### for flight departure domestic at CNX ###

setwd("C:/Users/Lenovo/OneDrive - Chiang Mai University/Master degree/IS/เอาไปงาน conference/for github")

central2017 <- read.csv("central_2560_1to45.csv")
str(central2017)

dcentral2017 <- central2017$delay.time
is.numeric(dcentral2017)
x_bar <- mean(dcentral2017)
min(dcentral2017)
max(dcentral2017)

library(fitdistrplus)
library(MASS)
library(survival)
library(moments)
skewness(dcentral2017)
kurtosis(dcentral2017)

###1.exponential
fitexpocentral2017 <- fitdist(data = dcentral2017, distr = 'exp',method = 'mle')
summary(fitexpocentral2017) ##parameter rate = 0.06981604 
gofstat(fitexpocentral2017)

###find error for expodietribution
xi <- dcentral2017
y_e <- (1/14.32336)*(exp((-xi)/14.32336)) ##beta = 1/0.06981604 

library(Metrics)
mae(xi,y_e)
library(MLmetrics)
MAPE <- mean(abs((xi-y_e)/xi)) * 100
MAPE
rmse(xi,y_e)

###2.gamma distribution
fit_g_central_2017  <- fitdist(data = dcentral2017,distr = "gamma",method = 'mle')
summary(fit_g_central_2017)
gofstat(fit_g_central_2017)

###find error for gamma distribution
y_g <- (1/(gamma(1.6586571)*((8.636466)^1.6586571)))*((xi)^(1.6586571-1))*(exp(-xi/8.636466))
library(Metrics)
mae(xi,y_g)
library(MLmetrics)
MAPE_g <- mean(abs((xi-y_g)/xi)) * 100
MAPE_g
rmse(xi,y_g)

###3.Three-Parameter Gamma Distribution
library(FAdist)
fitthreegam <- fitdist(data = dcentral2017 , distr = "gamma3" , method = 'mle', start = list(shape = 0.1 , scale = 1000 , thres = -0.00000002))
summary(fitthreegam)
gofstat(fitthreegam)

####to find error
shape_3gam <- 2.075847 
scale_3gam <- 7.531365   
thres_3gam <-  -1.310950  

y_3gam <- (1/(gamma(shape_3gam)*((scale_3gam)^shape_3gam)))*((xi-thres_3gam)^(shape_3gam-1))*(exp((-xi-thres_3gam)/scale_3gam))

mae(xi,y_3gam)
MAPE_3gam <- mean(abs((xi-y_3gam)/xi)) * 100 
MAPE_3gam
rmse(xi,y_3gam) 

###4.Generalized gamma Distribution

dggamma <- function(t, theta, kappa, delta, log = FALSE){val <- log(delta) - kappa*log(theta) - lgamma(kappa/delta) + (kappa - 1)*log(t) - (t/theta)^delta 
if(log) return(val) else return(exp(val))}

pggamma <- function(t, theta, kappa, delta, log.p = FALSE){val <- pgamma( t^delta, shape = kappa/delta, scale = theta^delta, log.p = TRUE) 
if(log.p) return(val) else return(exp(val))}

qggamma <- function(p, theta, kappa, delta){out <- qgamma(p, shape = kappa/delta, scale = theta^delta)^(1/delta) 
return(out)}

fitgergam <- fitdist(data = dcentral2017, distr = 'ggamma',method = 'mle', start = list(theta = 1, kappa = 2, delta = 0.5))
summary(fitgergam)
gofstat(fitgergam)

###to find error
theta_gg <- 13.219675
kappa_gg <- 1.455802
delta_gg <- 1.228525

y_GG <- ((delta_gg/(theta_gg^kappa_gg))*(xi^(kappa_gg-1))*exp(-(xi/theta_gg)^delta_gg))/gamma(kappa_gg/delta_gg)

mae(xi,y_GG)
MAPE_GG <- mean(abs((xi-y_GG)/xi)) * 100 
MAPE_GG
rmse(xi,y_GG) 

###5.log normal distribution
fitlnorm <- fitdist(data = dcentral2017, distr = "lnorm",method = 'mle')
summary(fitlnorm)
gofstat(fitlnorm)

###to find error
library(SciViews)
y_log <- (1/(xi*0.8989893)*((2*pi)^(1/2)))*(exp((-(ln(xi)-2.3310009)^2)/(2*(0.8989893)^2)))
mae(xi,y_log)
MAPE_log <- mean(abs((xi-y_log)/xi)) * 100
MAPE_log
rmse(xi,y_log)

###6.weibull Distribution
fitweil <- fitdist(data = dcentral2017,distr = "weibull",method = 'mle')
summary(fitweil)
gofstat(fitweil)

###to find error
gamma_w <- 1.362971
alpha_w <- 15.667395

y_weilbull <- gamma_w/alpha_w*(xi/alpha_w)^(gamma_w-1)*exp(-(xi/alpha_w)^gamma_w)

mae(xi,y_weilbull)
MAPE_weilbull <- mean(abs((xi-y_weilbull)/xi)) * 100 
MAPE_weilbull
rmse(xi,y_weilbull) 

###7.exponential weilbull distribution 
pexpweibull<- function(t,lambda,kappa,alpha,log.p=FALSE){
  log.cdf <- alpha*pweibull(t,scale=lambda,shape=kappa,log.p=TRUE)
  ifelse(log.p, return(log.cdf), return(exp(log.cdf)))
}  

dexpweibull<- function(t,lambda,kappa,alpha,log=FALSE){
  log.pdf <-  log(alpha) + (alpha-1)*pweibull(t,scale=lambda,shape=kappa,log=TRUE) + 
    dweibull(t,scale=lambda,shape=kappa,log=TRUE)
  ifelse(log, return(log.pdf), return(exp(log.pdf)))
}

qexpweibull<- function(p,lambda,kappa,alpha){
  quant <-  qweibull(p^(1/alpha),scale=lambda,shape=kappa)
  return(quant)
}  

fitexpoweil <- fitdist(data = dcentral2017 , distr = 'expweibull', method = 'mle',start = list(lambda = 5, kappa = 0.01, alpha = 0.5 ))
summary(fitexpoweil)
gofstat(fitexpoweil)

###to find error
alpha_ew <- 1.282796
seta_ew <- 1.110705
sigma_ew <- 14.600379

y_ew <- ((alpha_ew*seta_ew)/sigma_ew)*((1-exp(-(xi/sigma_ew)^alpha_ew))^(seta_ew-1)) * exp(-(xi/sigma_ew)^alpha_ew)*(xi/sigma_ew)^(alpha_ew-1)

mae(xi,y_ew)
MAPE_ew <- mean(abs((xi-y_ew)/xi)) * 100 
MAPE_ew
rmse(xi,y_ew)

###8.Generalized Weibull Distribution
dgenweibull<- function(t,sigma,alpha,lambda,log.p=FALSE){
  tsig <- t/sigma
  log.pdf <-   (1/alpha -1)*log(tsig) + (1/lambda - 1)*log( 1 - lambda*tsig^(1/alpha)) -
    log(alpha) - log(sigma)
  ifelse(log.p, return(log.pdf), return(exp(log.pdf)))
}

pgenweibull<- function(t,sigma,alpha,lambda,log.p=FALSE){
  log.cdf <- log( 1 - ( 1 - lambda*( t/sigma )^(1/alpha)  )^(1/lambda) )
  ifelse(log.p, return(log.cdf), return(exp(log.cdf)))
}  

qgenweibull<- function(p,sigma,alpha,lambda){
  if(lambda==0) quant <- sigma*( -log(1-p) )^alpha
  else quant <-  sigma*( (1-(1-p)^lambda)/lambda )^alpha
  return(quant)
}  

fitgenwei1 <- fitdist(data = dcentral2017 , distr = 'genweibull', method = 'mle', start = list(sigma = 12 , alpha = 4, lambda = 0.001 ))
summary(fitgenwei1)
gofstat(fitgenwei1)

###to find error
sigma_gw <-  17.1189093
alpha_gw <-  0.7969206
lambda_gw <-   0.1503987 

y_gw <- (((xi/sigma_gw)^((1/alpha_gw)-1))*(1-(lambda_gw*(xi/sigma_gw)^(1/alpha_gw)))^((1/lambda_gw)-1))/(alpha_gw*sigma_gw)

mae(xi,y_gw)
MAPE_gw <- mean(abs((xi-y_gw)/xi)) * 100 
MAPE_gw
rmse(xi,y_gw) 

###9.Burr Type III distribution
library(CoSMoS)
library(ggplot2)
library(data.table)

fitburiii <- fitdist(data = dcentral2017 , distr = "burrIII", method = 'mle', start = list(scale = 10 , shape1 = 2, shape2 = 3))
summary(fitburiii)
gofstat(fitburiii)

###to find error
y_testbu3 <- dburrIII(x = dcentral2017 , scale = 22.0534363,shape1 = 22.0534363 ,shape2 =  0.3020421)
mae(xi,y_testbu3)
MAPE_testbu3 <- mean(abs((xi-y_testbu3)/xi)) * 100 
MAPE_testbu3
rmse(xi,y_testbu3)

###10.Burr type XII Distribution
library(actuar)
fitburrxii <- fitdist(data = dcentral2017, distr = "burr", method = 'mle', start = list(shape1 = 10, shape2 = 0.009, rate = 0.001))
summary(fitburrxii)
gofstat(fitburrxii)

###to find error
shape1 <- 16.611302248 
shape2 <- 1.405895519
scale_b <- 1/0.008917437

y_burr <- ((shape2*shape1)/scale_b)*((xi/scale_b)^(shape2-1))*((1+((xi/scale_b)^shape2)))^(-(shape1+1))

mae(xi,y_burr)
MAPE_burr <- mean(abs((xi-y_burr)/xi)) * 100
MAPE_burr
rmse(xi,y_burr)

###11.Inverse Burr distribution
library(actuar)
fitinverseburr <- fitdist(data = dcentral2017, distr = 'invburr',method = 'mle',start = list(shape1 = 1, shape2 = 2, rate = 0.005))
summary(fitinverseburr)
gofstat(fitinverseburr)

###to find error
tao_ib <- 0.35514020
gamma_ib <- 3.31325874
seta_ib <- 21.00926

y_inverseburr <- ((tao_ib*gamma_ib*(xi/seta_ib))^(gamma_ib*tao_ib))/(xi*(1+(xi/seta_ib)^gamma_ib)^(tao_ib+1))

mae(xi,y_inverseburr)
MAPE_inverseburr <- mean(abs((xi-y_inverseburr)/xi)) * 100 
MAPE_inverseburr
rmse(xi,y_inverseburr) 

###12.Log pearson type 3 Distribution
dcentral2017log <- log10(dcentral2017)

##tofind parameter
library(EnvStats)
g <- skewness(dcentral2017log)
s <- sd(dcentral2017log)
m <- mean(dcentral2017log)

alpha <- (4/(g^2))
beta <- ((s*g)/2)
tao <- (m - (2*(s/g)))

my.param <- list(shape = alpha, scale = beta, location = tao)

dPIII<-function(x, shape, location, scale) PearsonDS::dpearsonIII(x, shape, location, scale, log=FALSE)
pPIII<-function(q, shape, location, scale) PearsonDS::ppearsonIII(q, shape, location, scale, lower.tail = TRUE, log.p = FALSE)
qPIII<-function(p, shape, location, scale) PearsonDS::qpearsonIII(p, shape, location, scale, lower.tail = TRUE, log.p = FALSE)

fitlogpearsontype3 <- fitdist(data = dcentral2017log , distr = "PIII", method = 'mle',start = my.param )

summary(fitlogpearsontype3)
gofstat(fitlogpearsontype3)

###to find error
xi_log <- dcentral2017log
xi_log_bar <- mean(dcentral2017log)

shape_p <-  3.7228126
scale_p <- -0.2097746 
location_p <- 1.7932910

y_pear <- (1/((abs(scale_p)^shape_p)*gamma(shape_p)))*(abs(xi_log-location_p)^(shape_p-1))*(exp(-(xi_log-location_p)/scale_p))

mae(xi_log,y_pear)
MAPE_pear <- mean(abs((xi_log-y_pear)/xi_log)) * 100
MAPE_pear
rmse(xi_log,y_pear)

###13.log logistic distribution
library(flexsurv)
fitloglogis <- fitdist(data = dcentral2017, distr = "llogis", method = 'mle')
summary(fitloglogis)
gofstat(fitloglogis)

###to find error
a <- 1.953738 
b <- 10.954721

y_loglo <- (a/b)*(((xi/b)^(a-1))/(1+((xi/b)^a))^2)

mae(xi,y_loglo)
MAPE_loglo <- mean(abs((xi-y_loglo)/xi)) * 100
MAPE_loglo
rmse(xi,y_loglo)

###14.generalized extreme value distribution
library(EnvStats)
library(lmomco)
dgev <- function(x,xi,alpha,kappa) {pdfgev(x,list(type="gev",para=c(xi,alpha,kappa),source="pargev"))}

pgev <- function(q,xi,alpha,kappa) {cdfgev(q,list(type="gev",para=c(xi,alpha,kappa),source="pargev"))}

qgev <- function(p,xi,alpha,kappa) {quagev(p,list(type="gev",para=c(xi,alpha,kappa),source="pargev"))}

lmr <- lmoms(dcentral2017 ,nmom = 5)
para.g <- pargev(lmr,checklmom = TRUE)

a=para.g[[2]][['xi']]
b=para.g[[2]][['alpha']]
c=para.g[[2]][['kappa']]

fitgev2 <- fitdist(dcentral2017, "gev",method='mle',start = list(xi=a,alpha=b,kappa=c))
summary(fitgev2)
gofstat(fitgev2)

###to find error
y_gev2 <- dgev(x = dcentral2017,xi = 8.6910061, alpha = 6.9238427,kappa = -0.2205725)

mae(xi,y_gev2)
MAPE_gev2 <- mean(abs((xi-y_gev2)/xi)) * 100 
MAPE_gev2 
rmse(xi,y_gev2) 

###15.Inverse gaussian distribution
library(actuar)
fitigcentral_2017 <- fitdist(data = dcentral2017, distr = "invgauss",method = 'mle',start = list(mean = 10,shape = 0.001))
summary(fitigcentral_2017)
gofstat(fitigcentral_2017)

###to find error
y_inverse <- (11.63639/(2*pi*(xi^3)))^(1/2)*exp((-11.63639 *(xi-14.32827)^2)/(2*((14.32827)^2)*xi))

mae(xi,y_inverse)
MAPE_in <- mean(abs((xi-y_inverse)/xi)) * 100
MAPE_in
rmse(xi,y_inverse)

####16.Generalized Inverse Gaussian Distribution
library(GeneralizedHyperbolic)
fitgig1 <- gigFit(x = dcentral2017 , freq = NULL ,paramStart = c(chi = 3 ,psi = 0.001 ,lambda = 0.09),method = 'Nelder-Mead')
summary(fitgig1)
fitgig1$param
###start parameter from estimate by Nelder-Mead

###use code for gig dis
my_param <- list(chi = 0.1387574, psi = 0.2275896 , lambda = 1.6193816 )
fitgig2 <- fitdist(data = dcentral2017, distr = "gig", method = 'mle',start = my_param)
summary(fitgig2)

###goodness of fit
library(DescTools)
AndersonDarlingTest(dcentral2017, null = "pgig") 

####to find error
y_gig <- dgig(x= dcentral2017, chi = 0.8491121,psi = 0.2130659, lambda = 1.4596320)

mae(xi,y_gig)
MAPE_gig <- mean(abs((xi-y_gig)/xi)) * 100 
MAPE_gig
rmse(xi,y_gig) 

###17.Beta prime distribution
library(extraDistr)
fitbetaprim <- fitdist(data = dcentral2017 , distr = 'betapr', method = 'mle' , start = list(shape1 = 0.0005 , shape2 = 50, scale = 100))
summary(fitbetaprim)
gofstat(fitbetaprim)

###to find error
alpha_bep <- 1.791996
beta_bep <- 16.919855 
sigma_bep <- 127.861665

y_bep <- (((xi/sigma_bep)^(alpha_bep-1))*((1+(xi/sigma_bep))^(-alpha_bep-beta_bep)))/(beta(alpha_bep,beta_bep)*sigma_bep)

mae(xi,y_bep)
MAPE_bep <- mean(abs((xi-y_bep)/xi)) * 100 
MAPE_bep
rmse(xi,y_bep) 

###18.Generalized beta distribution
library(actuar)
fitgenbeta <- fitdist(data = dcentral2017, distr = 'genbeta',method = 'mle', start = list(shape1 = 0.1, shape2 = 2, shape3 = 0.005, rate = 0.001))
summary(fitgenbeta)
gofstat(fitgenbeta)

###to find error
alpha1_GB <- 2.33826384
beta2_GB <- 4.87525627
tao3_GB <- 0.66253770
seta_GB <- 71.35004

y_GB <- ((gamma(alpha1_GB+beta2_GB))/(gamma(alpha1_GB)*gamma(beta2_GB)))*((xi/seta_GB)^(alpha1_GB*tao3_GB))*((1-((xi/seta_GB)^tao3_GB))^(beta2_GB-1))*(tao3_GB/xi)

mae(xi,y_GB)
MAPE_GB <- mean(abs((xi-y_GB)/xi)) * 100 
MAPE_GB
rmse(xi,y_GB) 

###19.The Generalized Beta Distribution of the Second Kind
library(GB2)
fitgenb2 <- fitdist(data = dcentral2017, distr = 'gb2',method = 'mle', start = list(shape1 = 0.1, scale = 10, shape2 = 5, shape3 = 3))
summary(fitgenb2)
gofstat(fitgenb2)

###to find error
r1 <-  0.8936827
r2 <- 8.1747145
r3 <- 1.5285312
b_gb2 <- 64.5940695

y_gb2 <- (r3/(b_gb2*beta(r1,r2)))*((xi/b_gb2)^(r1*r3-1))*(1+(xi/b_gb2)^r3)^(-(r2+r1))

mae(xi,y_gb2)
MAPE_gb2 <- mean(abs((xi-y_gb2)/xi)) * 100 
MAPE_gb2
rmse(xi,y_gb2) 

####20.Generalized pareto distribution
library(extraDistr)

fitgenpare <- fitdist(data = dcentral2017, distr = 'gpd', method = 'mle', start = list(mu = -0.1, sigma = 0.00000008 , xi = 20))
summary(fitgenpare)
gofstat(fitgenpare)

###to find error
k_gp <-  -0.2491820
sigma_gp <- 15.6945328
mu_gp <- 0.8719303

y_gp <- (1/sigma_gp)*(1+(k_gp*((xi-mu_gp)/sigma_gp)))^(-1-(1/k_gp))

mae(xi,y_gp)
MAPE_gp <- mean(abs((xi-y_gp)/xi)) * 100 
MAPE_gp
rmse(xi,y_gp) 

###for year 2018 and year 2019. you can change data ###
### R codes to find each distributions same year 2017 ####
