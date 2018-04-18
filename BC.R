
setwd("C:/Users/HimChan/Desktop/DAction")
library(nnet)
library(MASS)
train1 <- read.csv("train.csv")
train <- train1[,c(2:4,10:15,17:21)]
rm(train1)
head(train)
test1 <- read.csv("test.csv")
test <- test1[,c(2:4,10:15,17:21)]
rm(test1)
head(test)

set.seed(1)

ltFreq <- train$FreqBC
adjltFreq <- (ltFreq - min(ltFreq)) / (max(ltFreq)-min(ltFreq))
nnet.freq.fit  <- nnet(adjltFreq~.,data=train[-c(1,2,3,6,13,14)], size=10,linout=TRUE)
nnet.freq.pred <- predict(nnet.freq.fit, newdata = test[-c(1,2,3,6,13,14)])*(max(ltFreq)-min(ltFreq))+min(ltFreq)
nnet.freq.pred <- nnet.freq.pred *(nnet.freq.pred >0)
glm.freq.pois  <- glm(FreqBC~.,data=train[-c(1,2,3,6,13)], family="poisson")
glm.freq.zip   <- zeroinfl(FreqBC~.,data=train[-c(1,2,3,6,13)])
glm.freq.zip_pred  <- predict(glm.freq.zip, test[-c(1,2,3,6,13,14)])

glm.freq.nb    <- glm.nb(FreqBC~.,data=train[-c(1,2,3,6,13)])
glm.freq.pois_pred  <- exp(predict(glm.freq.pois, test[-c(1,2,3,6,13,14)]))
glm.freq.nb_pred    <- exp(predict(glm.freq.nb,   test[-c(1,2,3,6,13,14)]))


## MSE

trainp <- subset(train,log(yAvgBC)>0)
head(trainp)
ltClaim <- log(trainp$yAvgBC)
adjltClaim <- (ltClaim - min(ltClaim)) / (max(ltClaim)-min(ltClaim))
wg <- trainp[14]
glm.avgsev_indep <- glm(trainp$yAvgBC~.,data=trainp[-c(1,2,3,6,13,14)], family=Gamma(link="log"))
glm.avgsev_indep.pred <- exp(predict(glm.avgsev_indep, test[-c(1,2,3,6,13,14)]))

glm.avgsev_dep <- glm(trainp$yAvgBC~.,data=trainp[-c(1,2,3,6,13)], family=Gamma(link="log"))
test.avg.pois <- cbind(test[-c(1,2,3,6,13,14)],glm.freq.pois_pred)
colnames(test.avg.pois)[9] <- "FreqBC"
glm.avgsev.pois_pred <- exp(predict(glm.avgsev_dep, test.avg.pois))

test.avg.zip <- cbind(test[-c(1,2,3,6,13,14)],glm.freq.zip_pred)
colnames(test.avg.zip)[9] <- "FreqBC"
glm.avgsev.zip_pred <- exp(predict(glm.avgsev_dep, test.avg.zip))

test.avg.nb <- cbind(test[-c(1,2,3,6,13,14)],glm.freq.nb_pred)
colnames(test.avg.nb)[9] <- "FreqBC"
glm.avgsev.nb_pred <- predict(glm.avgsev_dep, test.avg.nb)
glm.avgsev.nb_pred[22] <- 0

nnet.avgsev_dep.fit <- nnet(adjltClaim~.,data=trainp[-c(1,2,3,6,13)], size=10,linout=TRUE)
test.avg_dep.nnet <- cbind(test[-c(1,2,3,6,13,14)],nnet.freq.pred)
colnames(test.avg_dep.nnet)[9] <- "FreqBC"
nnet.avgsev_dep.pred <- exp(predict(nnet.avgsev_dep.fit, test.avg_dep.nnet)/
                            (max(ltClaim)-min(ltClaim)) +min(ltClaim))

nnet.avgsev_indep.fit <- nnet(adjltClaim~.,data=trainp[-c(1,2,3,6,13,14)], size=10,linout=TRUE)
test.avg_indep.nnet <- cbind(test[-c(1,2,3,6,13,14)])
nnet.avgsev_indep.pred <- exp(predict(nnet.avgsev_indep.fit, test.avg_indep.nnet)/
                                (max(ltClaim)-min(ltClaim)) +min(ltClaim))

tClaim <- log(train$ClaimBC+1)
adjtClaim <- (tClaim - min(tClaim)) / (max(tClaim)-min(tClaim))
nnet.tsev.fit <- nnet(adjtClaim~.,data=train[-c(1,2,3,6,13,14)], size=5,linout=TRUE)
nnet.tsev.pred <- exp(predict(nnet.tsev.fit, newdata = test[-c(1,2,3,6,13,14)])* (max(tClaim)-min(tClaim)))-1
nnet.tsev.pred <- nnet.tsev.pred * (nnet.tsev.pred > 0)

library(tweedie)
require(statmod)

twBC<-glm(ClaimBC~TypeCity+TypeCounty+TypeSchool+TypeTown
          +TypeVillage+CoverageBC+lnDeductBC+NoClaimCreditBC,
          data=train[-c(1,2,13,14)],family=tweedie(var.power=1.5, link.power=0))
glm.tsev.pred <- exp(predict(twBC, newdata = test[-c(1,2,13,14)]))

summary(twBC)
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
par(mfrow=c(1,1),mar=c(1,1,4,1),mgp=c(2.2,1,0),oma=c(0,0,0,0))
plot.nnet(nnet.tsev.fit)

sqrt(mean((glm.avgsev_indep.pred *glm.freq.pois_pred - test$ClaimBC)^2)) # indep two-part glm (pois + gamma)
sqrt(mean((glm.avgsev.pois_pred  *glm.freq.pois_pred - test$ClaimBC)^2)) # dep two-part glm (pois + gamma)
sqrt(mean((glm.avgsev_indep.pred *glm.freq.zip_pred - test$ClaimBC)^2)) # indep two-part glm (zip + gamma)
sqrt(mean((glm.avgsev.zip_pred  *glm.freq.zip_pred - test$ClaimBC)^2)) # dep two-part glm (zip + gamma)
sqrt(mean((nnet.avgsev_indep.pred*nnet.freq.pred - test$ClaimBC)^2))     # indep two-part nnet
sqrt(mean((nnet.avgsev_dep.pred  *nnet.freq.pred - test$ClaimBC)^2))     # dep two-part nnet
sqrt(mean((nnet.tsev.pred - test$ClaimBC)^2))                            # one-part nnet
sqrt(mean((glm.tsev.pred - test$ClaimBC)^2))                             # tweedie glm

source("Gindex.R")

Px <- rep(1,length(test$ClaimBC))
par(mfrow=c(2,3),mar=c(1,1,4,1),mgp=c(2.2,1,0),oma=c(0,0,0,0))
gingraph.poisgam_indep <- gini.index.graphic(glm.avgsev_indep.pred*glm.freq.pois_pred,test$ClaimBC,Px,Title="Indep Pois-Gamma GLM")
gingraph.poisgam_dep <- gini.index.graphic(glm.avgsev.pois_pred*glm.freq.pois_pred,test$ClaimBC,Px,Title="Dep Pois-Gamma GLM")
gingraph.zipgam_indep <- gini.index.graphic(glm.avgsev_indep.pred*glm.freq.pois_pred,test$ClaimBC,Px,Title="Indep ZIP-Gamma GLM")
gingraph.zipgam_dep <- gini.index.graphic(glm.avgsev.pois_pred*glm.freq.pois_pred,test$ClaimBC,Px,Title="Dep ZIP-Gamma GLM")
gingraph.2nnet_indep <- gini.index.graphic(nnet.avgsev_indep.pred*nnet.freq.pred,test$ClaimBC,Px,Title="Indep Two-parts NN")
gingraph.2nnet_dep <- gini.index.graphic(nnet.avgsev_dep.pred*nnet.freq.pred,test$ClaimBC,Px,Title="Dep Two-parts NN")
gingraph.Tweedglm <- gini.index.graphic(glm.tsev.pred,test$ClaimBC,Px,Title="Tweedie GLM")
gingraph.1nnet <- gini.index.graphic(nnet.tsev.pred,test$ClaimBC,Px,Title="One-part NN")

####

MSEs <- rep(0,300)

for (i in 2007:2010) {
  
  trn <- subset(train,Year=!i)
  trnp <- subset(trainp,Year=!i)
  tstp <- subset(trainp,Year==i)
  
  sevglm <- glm(trnp$yAvgBC~.,data=trnp[-c(1,2,3,6,13)], family=Gamma(link="log"),weights = trnp$FreqBC)
  glmbeta <- coefficients(sevglm)
  
  xx <- tstp[-c(1,2,3,6,13)]
  xx <- cbind(rep(1,nrow(xx)),xx)
  c <- tstp$yAvgBC
  
  for (j in 1:300) {
    coef  <- GPEst(c=trnp$yAvgBC,x=trnp[-c(1,2,3,6,13)],id=trnp$PolicyNum,n=trnp$FreqBC,init.beta=glmbeta,k = j)$coef
    s_pred <- exp(as.matrix(xx) %*% coef[1:ncol(xx)]) * tstp$FreqBC
    MSEs[j] <- MSEs[j] + sqrt(mean(s_pred-tstp$FreqBC*c)^2) / 4
  }
}

GPEst <- function(c,x,n,id,init.beta,init.phi=1,k) {
  x <- cbind(rep(1,nrow(x)),x)
  colnames(x)[1] <- "intercept"
  
  e <- ncol(x)
  "negll.GB2" <- function(parm) {
    e <- ncol(x);
    reg_eqn <- as.matrix(x) %*% parm[1:e];
    data <- cbind(id,n*c/exp(reg_eqn),n);
    colnames(data)[2] <- "sv";
    
    temp1 = (sum(lgamma(as.matrix(aggregate(n~id,data,sum))[,2]/parm[e+1]+k/parm[e+1]+1))-sum(log(c))
            +sum(n*(log(n*c)-reg_eqn-log(parm[e+1])))/parm[e+1]+length(unique(id))*(k/parm[e+1]+1)*log(k/parm[e+1])
            -sum(lgamma(n/parm[e+1]))-length(unique(id))*lgamma(k/parm[e+1]+1))
    temp2 = (-sum((as.matrix(aggregate(n~id,data,sum))[,2]/parm[e+1]+k/parm[e+1]+1)*
             log(as.matrix(aggregate(sv~id,data,sum))[,2]/parm[e+1]+k/parm[e+1])))
    result = -temp1-temp2 + (k/parm[e+1]<=1)*10000000000
    return(result) }
  init.est <- as.vector(c(init.beta,init.phi))
  
  fit.GB2 <- optim(init.est, negll.GB2, NULL)
  parm.hat <- fit.GB2$par

  # next estimate the standard errors.
  library(nlme)
  #negll.GB2.Hess <- fdHess(parm.hat, negll.GB2);
  #inv.GB2.Hess <- solve(negll.GB2.Hess$Hessian);
  #parm.se <- sqrt(diag(inv.GB2.Hess));
  # put together the model with the est, se, t, pval, AIC, BIC
  #dfe <- length(c-length(parm.hat));
  #t_ratio<-parm.hat/parm.se;
  ##test if diff. from 1 t_ratio[1:3]<-(parm.hat[1:3]-1)/parm.se[1:3];
  #pval <- pf(t_ratio*t_ratio,df1=1,df2=dfe,lower.tail=F);
  #ttable <- cbind(parm.hat,parm.se,t_ratio,pval) 
  #ttable <- round(ttable,digits=4)
  
  #rownames(ttable)<- c(colnames(x),"k")
  #colnames(ttable)<- c("estimate", "std error", "t-val","Pr>|t|");
  
  AIC<- 2*negll.GB2(parm.hat) + 2*length(parm.hat);
  BIC<- 2*negll.GB2(parm.hat) + log(length(c))*length(parm.hat);
  loglik <- -negll.GB2(parm.hat)
  #return(list(ttable=ttable,AIC=AIC,BIC=BIC,loglik=loglik,coef=parm.hat));
  return(list(AIC=AIC,BIC=BIC,loglik=loglik,coef=parm.hat));
}

x <- trainp[-c(1,2,3,6,13)]
id <- trainp$PolicyNum
x <- cbind(rep(1,nrow(x)),x)
n <- trainp$FreqBC
c <- trainp$yAvgBC
glmbeta <- coefficients(glm.avgsev_dep)

reg_eqns <- as.matrix(x) %*% glmbeta
datas <- cbind(id,n*c/exp(reg_eqns),n);
colnames(datas)[2] <- "sv";

GPm <- GPEst(c=trainp$yAvgBC,x=trainp[-c(1,2,3,6,13)],id=trainp$PolicyNum,n=trainp$FreqBC,init.beta=glmbeta,k=KK)
gp.avgsev.pois_pred <- exp(as.matrix(test.avg.pois) %*% GPm$coef[1:ncol(test.avg.pois)])

reg_eqn <- as.matrix(x) %*% GPm$coef[1:ncol(x)]
data <- cbind(id,n*c/exp(reg_eqn),n);
colnames(data)[2] <- "sv";
post <- aggregate(n~id,data,sum)
post$sv <- aggregate(sv~id,data,sum)[,2]
post$weight <- (post$sv + KK) / (post$n + KK)
colnames(post)[1] <- "PolicyNum"
post$n <- NULL
post$sv <- NULL

postweight <- merge(x = test, y = post, by = "PolicyNum", all.x = TRUE)
postweight <- postweight[,ncol(postweight)]
postweight[is.na(postweight)] <- 1

test.avg.pois <- cbind(rep(1,nrow(test)),test[-c(1,2,3,6,13,14)],glm.freq.pois_pred)
colnames(test.avg.pois)[1] <- "intercept"
colnames(test.avg.pois)[ncol(test.avg.pois)] <- "FreqBC"

gp.avgsev.pois_pred <- exp(as.matrix(test.avg.pois) %*% GPm$coef[1:ncol(test.avg.pois)])
sqrt(mean((gp.avgsev.pois_pred  *glm.freq.pois_pred - test$ClaimBC)^2)) # dep two-part (pois + gp)
sqrt(mean((gp.avgsev.pois_pred*glm.freq.pois_pred*postweight - test$ClaimBC)^2)) # dep two-part (pois + gp)


########3########3########3########3########3########3########3########3########3########3########3

GB2Est <- function(c,x,n,id,init.p=1,init.beta,init.phi=1,k) {
  x <- cbind(rep(1,nrow(x)),x)
  colnames(x)[1] <- "intercept"
  
  e <- ncol(x)
  
  "negll.GB2" <- function(parm) {
      e <- ncol(x);
      reg_eqn <- as.matrix(x) %*% parm[1:e];
      sv <- (c/exp(reg_eqn)*exp(lgamma(n/parm[e+1]+1/parm[e+2]) - lgamma(n/parm[e+1])) )^parm[e+2] ;
      data <- cbind(id,sv,n);
      colnames(data)[2] <- "sv";
      
      temp1 = ( length(unique(id))*( parm[e+2]*(k/parm[e+1]+1)*(lgamma(k/parm[e+1]+1)-lgamma(k/parm[e+1]+1-1/parm[e+2]))  -lgamma(k/parm[e+1]+1))        # (4)
                -sum(log(c/parm[e+2]))+sum(lgamma(as.matrix(aggregate(n~id,data,sum))[,2]/parm[e+1]+k/parm[e+1]+1))-sum(lgamma(n/parm[e+1])) # (1)
                +parm[e+2]*sum( n*( log(c)-reg_eqn+ lgamma(n/parm[e+1]+1/parm[e+2])-lgamma(n/parm[e+1]))) /parm[e+1] )                             # (2)
      
      temp2 = (-sum((as.matrix(aggregate(n~id,data,sum))[,2]/parm[e+1]+k/parm[e+1]+1)*log(as.matrix(aggregate(sv~id,data,sum))[,2]
                                                                                        +exp(lgamma(k/parm[e+1]+1)-lgamma(k/parm[e+1]+1-1/parm[e+2]))^parm[e+2])))                                                                  # (3)
      result = -temp1-temp2 + (k/parm[e+1]<=1)*10000000000 + (parm[e+1]<=0)*10000000000
      return(result) }
  init.est <- as.vector(c(init.beta,init.phi,init.p))
  
  fit.GB2 <- optim(init.est, negll.GB2, NULL)
  parm.hat <- fit.GB2$par
  
  # next estimate the standard errors.
  library(nlme)
  #negll.GB2.Hess <- fdHess(parm.hat, negll.GB2);
  #inv.GB2.Hess <- solve(negll.GB2.Hess$Hessian);
  #parm.se <- sqrt(diag(inv.GB2.Hess));
  # put together the model with the est, se, t, pval, AIC, BIC
  #dfe <- length(c-length(parm.hat));
  #t_ratio<-parm.hat/parm.se;
  ##test if diff. from 1 t_ratio[1:3]<-(parm.hat[1:3]-1)/parm.se[1:3];
  #pval <- pf(t_ratio*t_ratio,df1=1,df2=dfe,lower.tail=F);
  #ttable <- cbind(parm.hat,parm.se,t_ratio,pval) 
  #ttable <- round(ttable,digits=4)
  
  #rownames(ttable)<- c(colnames(x),"Phi","k")
  #colnames(ttable)<- c("estimate", "std error", "t-val","Pr>|t|");
  
  AIC<- 2*negll.GB2(parm.hat) + 2*length(parm.hat);
  BIC<- 2*negll.GB2(parm.hat) + log(length(c))*length(parm.hat);
  loglik <- -negll.GB2(parm.hat)
  #return(list(ttable=ttable,AIC=AIC,BIC=BIC,loglik=loglik,coef=parm.hat));
  return(list(AIC=AIC,BIC=BIC,loglik=loglik,coef=parm.hat)); }

gb2aic <- rep(0,100)
for (i in 1:100) {
  tryCatch({
    temp <- GB2Est(c=trainp$yAvgBC,x=trainp[-c(1,2,3,6,13)],id=trainp$PolicyNum,p=i/10,n=trainp$FreqBC,init.beta=glmbeta)
    gb2aic[i] <- temp$AIC
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

gb2aicx <- rep(0,20)
for (i in 1:20) {
  tryCatch({
    temp <- GB2Est(c=trainp$yAvgBC,x=trainp[-c(1,2,3,6,13)],id=trainp$PolicyNum,p=i/100+0.4,n=trainp$FreqBC,init.beta=glmbeta)
    gb2aicx[i] <- temp$AIC
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

x <- trainp[-c(1,2,3,6,13)]
id <- trainp$PolicyNum
x <- cbind(rep(1,nrow(x)),x)
n <- trainp$FreqBC
c <- trainp$yAvgBC
glmbeta <- coefficients(glm.avgsev_dep)

reg_eqns <- as.matrix(x) %*% glmbeta
datas <- cbind(id,n*c/exp(reg_eqns),n);
colnames(datas)[2] <- "sv";

GB2m <- GB2Est(c=trainp$yAvgBC,x=trainp[-c(1,2,3,6,13)],id=trainp$PolicyNum,n=trainp$FreqBC,init.beta=glmbeta,k=KK)
GB2.avgsev.pois_pred <- exp(as.matrix(test.avg.pois) %*% GB2m$coef[1:ncol(test.avg.pois)])

phi <- GB2m$coef[11]
p <- GB2m$coef[12]

reg2_eqn <- as.matrix(x) %*% GB2m$coef[1:ncol(x)]
data2 <- cbind(id,exp(p*log(c)-p*reg2_eqn+p*lgamma(n/phi+1/p)-p*lgamma(n/phi)) ,n);
colnames(data2)[2] <- "sv";
post2 <- aggregate(n~id,data2,sum)
post2$sv <- aggregate(sv~id,data2,sum)[,2]
post2$weight <- ( exp(lgamma(post2$n/phi+KK/phi+1-1/p)-lgamma(post2$n/phi+KK/phi+1))*
                  (post2$sv+exp(p*lgamma(KK/phi+1)-p*lgamma(KK/phi+1-1/p)) )^(1/p) )

colnames(post2)[1] <- "PolicyNum"
post2$n <- NULL
post2$sv <- NULL

post2weight <- merge(x = test, y = post2, by = "PolicyNum", all.x = TRUE)
post2weight <- post2weight[,ncol(post2weight)]
post2weight[is.na(post2weight)] <- 1

sqrt(mean((GB2.avgsev.pois_pred  *glm.freq.pois_pred - test$ClaimBC)^2)) # dep two-part glm (pois + GB2)
sqrt(mean((GB2.avgsev.pois_pred  *glm.freq.pois_pred*post2weight - test$ClaimBC)^2)) # dep two-part glm (pois + GB2)


residu <- function(x,c,n,id,beta,wfactor=0) {
  x <- cbind(rep(1,nrow(x)),x)
  mu <- exp(as.matrix(x) %*% beta[1:ncol(x)])
  data <- cbind(id,n*c/mu,n);
  post <- aggregate(n~id,data,sum)
  colnames(data)[2] <- "sv";
  colnames(data)[1] <- "PolicyNum";
  post$sv <- aggregate(sv~id,data,sum)[,2]
  weight <- (post$sv + wfactor) / (post$n + wfactor)
  return(list(stable=post,weight=weight)); }


wgt <- residu(x=trainp[-c(1,2,3,6,13)], id=trainp$PolicyNum, n=trainp$FreqBC,
              c=trainp$yAvgBC,beta=coefficients(glm.avgsev_dep))
wtable <- as.matrix(wgt$stable)
M <- which.max(wgt$weight)