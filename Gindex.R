#
#  COMPUTE SINGLE SCORE GINI, STANDARD ERROR, PROJECTION
# mx = Estimated Loss
# Sx = True Loss
ScoreFunc <- function(mx,Sx){
  n = length(mx)
  origorder = (1:n)
  PSRmat <- data.frame(cbind(mx,Sx,origorder))
  PSRmatOrder <- PSRmat[order(mx),]  #  Sort by relativity
  #  PREMIUM, LOSS DFs
  DFPredLoss = cumsum(PSRmatOrder$mx)/sum(mx)
  DFActualLoss = cumsum(PSRmatOrder$Sx)/sum(Sx)
  Bound = PSRmatOrder$mx*DFPredLoss*(sum(mx)/sum(Sx))
  #  Bound = ((1:n)/n)*DFPrem*(sum(PIx)/sum(Sx))
  #  GINI CALC
  DFdiff = DFPredLoss[2:n]-DFPredLoss[1:(n-1)]
  DFPredavg  = (DFPredLoss[2:n]+DFPredLoss[1:(n-1)])/2
  DFActualavg  = (DFActualLoss[2:n]+DFActualLoss[1:(n-1)])/2
  (Gini = 2*crossprod(DFdiff,DFPreavg-DFActualavg))
  #  PROJECTION CALC
  meany = mean(Sx)
  meanpi = mean(PIx)
  h1 = 0.5* (mean(Sx)*PSRmatOrder$PIx*DFLoss +
               PSRmatOrder$Sx*meanpi*(1-DFPrem) )
  #  STANDARD ERROR CALC
  h1bar = mean(h1)
  sigmah = var(h1)
  sigmahy = cov(h1,PSRmatOrder$Sx)
  sigmahpi = cov(h1,PSRmatOrder$PIx)
  sigmay = var(Sx)
  sigmapi = var(PIx)
  sigmaypi = cov(PSRmatOrder$Sx,PSRmatOrder$PIx)
  temp1 = 4*sigmah + (h1bar/meany)^2*sigmay + (h1bar/meanpi)^2*sigmapi -
    4*(h1bar/meany)*sigmahy -4 *(h1bar/meanpi)*sigmahpi +
    2*(h1bar^2/(meany*meanpi))*sigmaypi
  sigmaGini = 4*temp1/(meany^2*meanpi^2)
  stderrGini = sqrt(sigmaGini/n)
  check = var(PIx-mx)
  Gini = Gini*(check != 0)
  stderrGini = stderrGini*(check != 0)
  h1mat <- data.frame(cbind(h1,PSRmatOrder$origorder))
  colnames(h1mat) = c("h1","origorder")
  h1matOrder <- h1mat[order(h1mat$origorder),]  #  Sort by original order
  Retmat <- data.frame(cbind(DFPrem,DFLoss,Bound))
  RetmatGini<-list(Retmat,Gini,stderrGini,h1matOrder$h1)
  return(RetmatGini)
}
#
#  COMPUTE DIFFERENCES IN GINIs FROM TWO SCORES
#
DiffScores <- function(PIx,SAx,SBx,Sx){
  n = length(PIx)
  tempA=ScoreFunc(PIx=PIx,mx=SAx ,Sx=Sx)
  GiniA = tempA[[2]]
  h1A=tempA[[4]]
  sigmaGiniA = n*(tempA[[3]])^2
  checkA=(sigmaGiniA != 0)
  tempB=ScoreFunc(PIx=PIx,mx=SBx ,Sx=Sx)
  GiniB = tempB[[2]]
  h1B=tempB[[4]]
  sigmaGiniB = n*(tempB[[3]])^2
  checkB=(sigmaGiniB != 0)
  #  ASYMPTOTIC COVARIANCE MATRIX
  h1Abar = mean(h1A);h1Bbar = mean(h1B)
  meany = mean(Sx); meanpi = mean(PIx)
  sigmahAhB = cov(h1A,h1B)
  sigmahAy = cov(h1A,Sx)
  sigmahBy = cov(h1B,Sx)
  sigmahApi = cov(h1A,PIx)
  sigmahBpi = cov(h1B,PIx)
  sigmay = var(Sx)
  sigmapi = var(PIx)
  sigmaypi = cov(Sx,PIx)
  muA = meany*meanpi*(1 - GiniA)/2
  muB = meany*meanpi*(1 - GiniB)/2
  sigmaGiniAB=(4*(4*sigmahAhB - 2*muB*sigmahAy/meany - 2*muB*sigmahApi/meanpi -
                    2*muA*sigmahBy/meany - 2*muA*sigmahBpi/meanpi +
                    muA*muB*sigmay/(meany^2) + 2*muA*muB*sigmaypi/(meany*meanpi) +
                    muA*muB*sigmapi/(meanpi^2) ) / (meany^2 * meanpi^2))*
    checkA*checkB
  stderrdiff = 100*sqrt((sigmaGiniA+sigmaGiniB-2*sigmaGiniAB)/n)
  stderrGiniA = 100*sqrt(sigmaGiniA/n)
  stderrGiniB = 100*sqrt(sigmaGiniB/n)
  tAB = 100*(GiniB - GiniA) / stderrdiff
  corrAB = sigmaGiniAB/sqrt(sigmaGiniA*sigmaGiniB)
  RetDiff <- cbind(100*GiniA,stderrGiniA,100*GiniB,stderrGiniB,
                   sigmaGiniAB,stderrdiff,tAB,corrAB)
  colnames(RetDiff) = c("GiniA","stderrGiniA","GiniB","stderrGiniB",
                        "sigmaGiniAB","stderrdiff","tAB","corrAB")
  return(RetDiff)
}
gini.index.graphic<-function(mx,Sx,Px,Title=NULL,premlab="Predicted pure premium",rellab="Relativity"){
  
  o<-order(mx/Px)
  pct.pred<-cumsum(Px[o])/sum(Px)
  pct.loss<-cumsum(Sx[o])/sum(Sx)
  
  delta.loss<-pct.loss-c(0,pct.loss[-(length(pct.loss))])
  riemann.segment<-(pct.pred-pct.loss)*delta.loss*200
  # the 50 is there to express gini as a percentage
  gini<-sum(riemann.segment)
  par(plt=c(.2,.8,.2,.8),font.axis=1,font.lab=1)
  plot(pct.pred,pct.loss,type="n",ylab="Actual Loss",xlim = c(0,1),ylim=c(0,1),
       xlab=premlab,main=Title,font.main=1)
  mtext(paste("Gini index = ",round(gini,digits=1)),cex=0.75)
  polygon(pct.pred,pct.loss,col="grey")
  lines(pct.pred,pct.loss,lwd=1,col="blue")
  abline(0,1,lwd=1)
 # legend("topleft",c("Actual","Predicted"),col=c("black","blue"), lwd=c(3,3))
  return(gini)
} # end gini function
