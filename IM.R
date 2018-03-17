setwd("C:/Users/HimChan/Desktop/DAction")
library(nnet)
library(MASS)
train1 <- read.csv("train.csv")
train <- train1[,c(2:3,5,10:15,22:26)]
rm(train1)
head(train)

test1 <- read.csv("test.csv")
test <- test1[,c(2:3,5,10:15,22:26)]
rm(test1)
head(test)


set.seed(50)

ltFreq <- train$FreqIM
adjltFreq <- (ltFreq - min(ltFreq)) / (max(ltFreq)-min(ltFreq))
nnet.freq.fit  <- nnet(adjltFreq~.,data=train[-c(1,2,3,6,13,14)], size=10,linout=TRUE)
nnet.freq.pred <- predict(nnet.freq.fit, newdata = test[-c(1,2,3,6,13,14)])*(max(ltFreq)-min(ltFreq))+min(ltFreq)
nnet.freq.pred <- nnet.freq.pred *(nnet.freq.pred >0)
glm.freq.pois  <- glm(FreqIM~.,data=train[-c(1,2,3,6,13)], family="poisson")
glm.freq.nb    <- glm.nb(FreqIM~.,data=train[-c(1,2,3,6,13)])
glm.freq.pois_pred  <- exp(predict(glm.freq.pois, test[-c(1,2,3,6,13,14)]))
glm.freq.nb_pred    <- exp(predict(glm.freq.nb,   test[-c(1,2,3,6,13,14)]))


## MSE

trainp <- subset(train,log(yAvgIM)>0)
head(trainp)
ltClaim <- log(trainp$yAvgIM)
adjltClaim <- (ltClaim - min(ltClaim)) / (max(ltClaim)-min(ltClaim))

glm.avgsev_indep <- glm(trainp$yAvgIM~.,data=trainp[-c(1,2,3,6,13,14)], family=Gamma(link="log"))
glm.avgsev_indep.pred <- exp(predict(glm.avgsev_indep, test[-c(1,2,3,6,13,14)]))

glm.avgsev_dep <- glm(trainp$yAvgIM~.,data=trainp[-c(1,2,3,6,13)], family=Gamma(link="log"))
test.avg.pois <- cbind(test[-c(1,2,3,6,13,14)],glm.freq.pois_pred)
colnames(test.avg.pois)[9] <- "FreqIM"
glm.avgsev.pois_pred <- exp(predict(glm.avgsev_dep, test.avg.pois))

test.avg.nb <- cbind(test[-c(1,2,3,6,13,14)],glm.freq.nb_pred)
colnames(test.avg.nb)[9] <- "FreqIM"
glm.avgsev.nb_pred <- predict(glm.avgsev_dep, test.avg.nb)
glm.avgsev.nb_pred[22] <- 0

nnet.avgsev_dep.fit <- nnet(adjltClaim~.,data=trainp[-c(1,2,3,6,13)], size=10,linout=TRUE)
test.avg_dep.nnet <- cbind(test[-c(1,2,3,6,13,14)],nnet.freq.pred)
colnames(test.avg_dep.nnet)[9] <- "FreqIM"
nnet.avgsev_dep.pred <- exp(predict(nnet.avgsev_dep.fit, test.avg_dep.nnet))

nnet.avgsev_indep.fit <- nnet(adjltClaim~.,data=trainp[-c(1,2,3,6,13,14)], size=10,linout=TRUE)
test.avg_indep.nnet <- cbind(test[-c(1,2,3,6,13,14)])
nnet.avgsev_indep.pred <- exp(predict(nnet.avgsev_indep.fit, test.avg_indep.nnet))

plot(log(test$ClaimIM+1), log(glm.avgsev.pred*glm.freq.pois_pred+1),
     main="GLM predictions vs actual",
     xlab="Actual",ylim=c(1,16),xlim=c(1,16))

plot(log(test$ClaimIM+1), log(nnet.avgsev.pred*nnet.freq.pred+1),
     main="Neural network predictions vs actual",
     xlab="Actual",ylim=c(1,16),xlim=c(1,16))

tClaim <- log(train$ClaimIM+1)
adjtClaim <- (tClaim - min(tClaim)) / (max(tClaim)-min(tClaim))
nnet.tsev.fit <- nnet(adjtClaim~.,data=train[-c(1,2,3,6,13,14)], size=5,linout=TRUE)
nnet.tsev.pred <- exp(predict(nnet.tsev.fit, newdata = test[-c(1,2,3,6,13,14)])* (max(tClaim)-min(tClaim)))-1
nnet.tsev.pred <- nnet.tsev.pred * (nnet.tsev.pred > 0)


library(tweedie)
require(statmod)

twIM<-glm(ClaimIM~TypeCity+TypeCounty+TypeSchool+TypeTown
          +TypeVillage+CoverageIM+lnDeductIM+NoClaimCreditIM,
          data=train[-c(1,2,13,14)],family=tweedie(var.power=1.5, link.power=0))
glm.tsev.pred <- exp(predict(twIM, newdata = test[-c(1,2,13,14)]))

summary(twIM)

sqrt(mean((glm.avgsev_indep.pred *glm.freq.pois_pred - test$ClaimIM)^2)) # indep two-part glm (pois + gamma)
sqrt(mean((glm.avgsev.pois_pred  *glm.freq.pois_pred - test$ClaimIM)^2)) # dep two-part glm (pois + gamma)
#sqrt(mean((glm.avgsev_indep.pred *glm.freq.nb_pred - test$ClaimIM)^2))   # indep two-part glm (nb + gamma)
#sqrt(mean((glm.avgsev.nb_pred    *glm.freq.nb_pred - test$ClaimIM)^2))   # dep two-part glm (nb + gamma)
sqrt(mean((nnet.avgsev_indep.pred*nnet.freq.pred - test$ClaimIM)^2))     # indep two-part nnet
sqrt(mean((nnet.avgsev_dep.pred  *nnet.freq.pred - test$ClaimIM)^2))     # dep two-part nnet
sqrt(mean((nnet.tsev.pred - test$ClaimIM)^2))                            # one-part nnet
sqrt(mean((glm.tsev.pred - test$ClaimIM)^2))                             # tweedie glm

source("Gindex.R")

Px <- rep(1,length(test$ClaimIM))
par(mfrow=c(2,3),mar=c(1,1,4,1),mgp=c(2.2,1,0),oma=c(0,0,0,0))
gingraph.poisgam_indep <- gini.index.graphic(glm.avgsev_indep.pred*glm.freq.pois_pred,test$ClaimIM,Px,Title="Indep Pois-Gamma GLM")
gingraph.poisgam_dep <- gini.index.graphic(glm.avgsev.pois_pred*glm.freq.pois_pred,test$ClaimIM,Px,Title="Dep Pois-Gamma GLM")
#gingraph.nbgam_indep <- gini.index.graphic(glm.avgsev_indep.pred*glm.freq.nb_pred,test$ClaimIM,Px,Title="Indep NB-Gamma GLM")
#gingraph.nbgam_dep <- gini.index.graphic(glm.avgsev.nb_pred*glm.freq.nb_pred,test$ClaimIM,Px,Title="Dep NB-Gamma GLM")
gingraph.2nnet_indep <- gini.index.graphic(nnet.avgsev_indep.pred*nnet.freq.pred,test$ClaimIM,Px,Title="Indep Two-parts NN")
gingraph.2nnet_dep <- gini.index.graphic(nnet.avgsev_dep.pred*nnet.freq.pred,test$ClaimIM,Px,Title="Dep Two-parts NN")
gingraph.Tweedglm <- gini.index.graphic(glm.tsev.pred,test$ClaimIM,Px,Title="Tweedie GLM")
gingraph.1nnet <- gini.index.graphic(nnet.tsev.pred,test$ClaimIM,Px,Title="One-part NN")
