suppressMessages(library(magrittr))
suppressMessages(library(readr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(Xmisc))
suppressMessages(library(pROC))
suppressMessages(library(rms))
#TODO unit testing on each of the functions.
#The stuff we regress on.
prs <- read_tsv("/work-zfs/abattle4/ashton/prs_dev/scratch/calc_prs_testing/8-28-tests/results_1.tsv_0.0001.tsv")

phenos <- read_tsv("/work-zfs/abattle4/ashton/prs_dev/CHS_race1_toy/pheno_mendrand.txt") %>% rename(ID = IID)

comp_dat <- inner_join(prs, phenos) %>% select(ID, Score, MI)
y <- comp_dat$MI
g <- comp_dat$Score
ncase <- sum(y)
nt <- length(g)
ncont <- nt-ncase
K <- 0.05 #population prevelance

thd <- -qnorm(K,0,1) #K of the population to the left
zv <- dnorm(thd) #the z score, normal desnity
mv <- zv/K  #mean liability
mv2 <- -mv*K/(1-K) #mean liability for healthy/control cases

#Calculate R2 on observed scale using the linear model
#@aram y: the outcome we are hoping to predict
#@param g: the score
#@param ncase: the number of cases that have the condition
#@param nt: the total number in the sample
#@param ncont: number of controls (healthy)

R2ObservedLinear <- function(y,g,n,ncase,ncont)
{
  lmv <- lm(y~g)
  R2_lin <- var(lmv$fitted.values)/(ncase/n*ncont/n)
  
  return(R2_lin)
}

#Get the Cox and Snell pseudo R on the observed scale
#@param n the number of ones in teh sample
R2ObservedCoxSnell <- function(y,g,n)
{ 
  #Cox and Snell R2 <this checks out with the paper>
  logf = logLik(glm(y~g, family = binomial(logit)))
  logn = logLik(glm(y~1, family = binomial(logit)))
  R2_cs = 1-exp((logn-logf)*(2/nt))
  return(R2_cs)
}

#Calculate Nagelkereke Psuedo-R2 on the observed scale
R2ObservedNagelkerke <-function(y,g)
{
  #Nagelkerke R2 <checks out>
  lrmv2 <- lrm(y~g) #logistic regression model
  R2_nagel <- lrmv2$stats[10] #what is this
  return(R2_nagel)
  
}

#Calcualte Pseudo-R2 on the probit liability scale
R2LiabilityProbit(y,g)
{
  #R2 on the probit liability sclae using probit model <checks out I think...>
  pmv <- glm(y~g, family = binomial(probit))
  R2_probit <- var(pmv$linear.predictors)/(var(pmv$linear.predictors) + 1)
  return(R2_probit)
}

#R2 on the logit liability scale
R2LiabilityLogit <- function(y,g)
{
  lls <- glm(y~g, family=binomial(logit))
  R2_logit <- var(lls$linear.predictors)/(var(lls$linear.predictors) + (pi^2)/3)
  return(R2_logit)
}

#R2 on the liability scale using AUC
R2LiabilityAUC <- function(y,g)
{
  aucv = auc(y, pmv$linear.predictors)
  qv = qnorm(aucv[1])
  R2_auc = 2*qv^2/((mv2-mv)^2 + qv^2*mv*(mv-thd) + mv2*(mv2-thd))
  return(R2_auc)
}

#Calculte the pseudo correlation on the liability scale using the transformation
#@param P: proportion of cases in the case-control samples
#@param K: population prevelance of the condition
R2LiabilityAscertained <- function(y,g,n,ncase,ncont,P,K)
{
  #Pre calculations
  thres <- -qnorm(K,0,1) #The value along normal distribution that captures the proportion in the population with condition 
  zval <- dnorm(thresh) #the y coordinate of this threshold, point on density curve
  mv <- zval/K #the mean liability for this case (#don't totally understand this)
  R2_lin <- R2ObservedLinear(y,g,n,ncase,ncont)
  theta <- mv * (P-K)/(1-K)*(mv*(P-K)/(1-K)-thd)
  C <- K*(1-K)/zv^2*K*(1-K)/(P*(1-P)) #C another correction factor
  R2_obs <- R2_lin*C / (1+R2_lin * theta*C)
  return(R2_obs)
}
