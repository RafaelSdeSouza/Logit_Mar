library(arm)
library(caret)
library(pROC)
require(plyr)
require(gam)
library(glmnet)
require(mgcv)
# Read and format data
data     <- read.table("..//data/matched.txt",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]


# Run for Spirals 
data_cut   <- subset(data_cut, zoo=="S")
x1<-data_cut$logM200_L
x2<-data_cut$RprojLW_Rvir
y<-data_cut$bpt
zoo<-as.factor(data_cut$zoo)
g<-gam(y~s(x2)+s(x1),family=binomial)
vis.gam(g,theta=-175,color="heat")
vis.gam(g,se=2,theta=-35,type="response",xlabel=c("oi"))

vis.gam(g, view=c("x1","x2"),plot.type="contour",color="topo")
plot(fit,pages=1,seWithMean=TRUE)




x<-as.matrix(data2[,2:5])
fit<-glmnet(x,y=data2$bpt,alpha=1,family="binomial")
plot(fit,xvar="lambda",label = TRUE)
plot(fit, xvar = "dev", label = TRUE)
cv.glmmod <- cv.glmnet(x,y=data2$bpt,alpha=1,family="binomial",type.measure = "class")
plot(cv.glmmod)
best_lambda <- cv.glmmod$lambda.min

coef.min = coef(cv.glmmod, s = "lambda.min")
active.min = coef.min[which(coef.min != 0)]
index.min = coef.min[active.min]






library(popbio)
logi.hist.plot(data2$logMstar,data2$bpt,boxp=FALSE,type="hist",col="gray")



ROCF<- data.frame(True=data2$bpt,predicted=predict(fit,type = "response"))
F1 <-roc(ROCF$True,ROCF$predicted)
coords(F1,x="best")[1]

ROCF$class<-ROCF$predicted
ROCF$class[which(ROCF$class>=coords(F1,x="best")[1])]<-1
ROCF$class[which(ROCF$class<coords(F1,x="best")[1])]<-0

plot(F1)
confusionMatrix(ROCF$True, ROCF$class)




x <-range(data2$logMstar)
x <- seq(0.5*x[1], 1.5*x[2], length.out=50)
y <- range(data2$vlos_sigma)
y <- seq(0.5*y[1], 1.5*y[2], length.out=50)

z <- outer(x,y,
           function(logMstar,vlos_sigma)
             predict(fit, data.frame(logMstar,vlos_sigma),type = 'response'))
library(rsm)
library(lattice)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#p<-persp(x,y,z, theta=150, phi=20,
#         expand = 0.5,shade = 0.1,
#         xlab="Z", ylab=expression(NII.Ha), zlab=expression(log10.EW.Ha),ticktype='detailed',
#         col = YlOrBr,border=NA,xlog=T,ylog=T)
cor = cm.colors(100)

cairo_pdf("logit3D_2.pdf")
trellis.par.set("axis.line",list(axis.text=list(cex=20),col=NA,lty=1,lwd=2))
par(mar=c(1,1,1,1))
wireframe(z~x+y,data=data.frame(x=x, y=rep(y, each=length(x)), z=z),
          par.settings = list(regions=list(alpha=0.6)),
          col.regions =cor,drape=T,light.source = c(5,5,5),colorkey = FALSE,
          xlab=list(label=expression(log[Mstar]),cex=1.25), ylab=list(label=expression(V[sigma]),cex=1.25),
          zlab=list(rot=90,label=expression(pi),cex=1.25,dist=-1,rot=0),
          scale=list(tck=0.75,arrows=FALSE,distance =c(0.75, 0.75, 0.75)))


dev.off()





market.size <- 800
icecream$opportunity <- market.size - icecream$units
bin.glm <- glm(cbind(units, opportunity) ~ temp, data=icecream,
               family=binomial(link = "logit"))
par(mfrow=c(2,2))
plot(bin.glm)
title(outer=TRUE, line = -1,
      main = list("Binomial (logit) GLM",
                  cex=1.25,col="black", font=2))

meanProb <- predict(bin.glm, type="response")
meanPred <- meanProb*market.size
UpPred <- qbinom(.95, market.size, meanProb)
LwPred <- qbinom(.05, market.size, meanProb)

plotData <- lapply(
  seq(along=icecream$temp),
  function(i){
    y = ylim[1]:ylim[2]
    x = rep(icecream$temp[i], length(y))
    z0 = rep(0, length(y))
    z = dbinom(y, market.size, meanProb[i])
    return(list(x=x, y=y, z0=z0, z=z))
  }
)



icecream <- data.frame(
  temp=c(11.9, 14.2, 15.2, 16.4, 17.2, 18.1,
         18.5, 19.4, 22.1, 22.6, 23.4, 25.1),
  units=c(185L, 215L, 332L, 325L, 408L, 421L,
          406L, 412L, 522L, 445L, 544L, 614L)
)
glmModelPlot <- function(x, y, xlim,ylim, meanPred,  LwPred, UpPred,
                         plotData, main=NULL){
  ## Based on code by Arthur Charpentier:
  ## http://freakonometrics.hypotheses.org/9593
  par(mfrow=c(1,1))
  n <- 2
  N <- length(meanPred)
  zMax <- max(unlist(sapply(plotData, "[[", "z")))*1.5
  mat <- persp(xlim, ylim, matrix(0, n, n), main=main,
               zlim=c(0, zMax), theta=-30,
               ticktype="detailed",box=FALSE)
  C <- trans3d(x, UpPred, rep(0, N),mat)
  lines(C, lty=2)
  C <- trans3d(x, LwPred, rep(0, N), mat)
  lines(C, lty=2)
  C <- trans3d(c(x, rev(x)), c(UpPred, rev(LwPred)),
               rep(0, 2*N), mat)
  polygon(C, border=NA, col=adjustcolor("yellow", alpha.f = 0.5))
  C <- trans3d(x, meanPred, rep(0, N), mat)
  lines(C, lwd=2, col="grey")
  C <- trans3d(x, y, rep(0,N), mat)
  points(C, lwd=2, col="#00526D")
  for(j in N:1){
    xp <- plotData[[j]]$x
    yp <- plotData[[j]]$y
    z0 <- plotData[[j]]$z0
    zp <- plotData[[j]]$z
    C <- trans3d(c(xp, xp), c(yp, rev(yp)), c(zp, z0), mat)
    polygon(C, border=NA, col="light blue", density=40)
    C <- trans3d(xp, yp, z0, mat)
    lines(C, lty=2)
    C <- trans3d(xp, yp, zp, mat)
    lines(C, col=adjustcolor("blue", alpha.f = 0.5))
  }
}

market.size <- 800
icecream$opportunity <- market.size - icecream$units
bin.glm <- glm(cbind(units, opportunity) ~ temp, data=icecream,
               family=binomial(link = "logit"))
par(mfrow=c(2,2))
plot(bin.glm)
title(outer=TRUE, line = -1,
      main = list("Binomial (logit) GLM",
                  cex=1.25,col="black", font=2))

meanProb <- predict(bin.glm, type="response")
meanPred <- meanProb*market.size
UpPred <- qbinom(.95, market.size, meanProb)
LwPred <- qbinom(.05, market.size, meanProb)

plotData <- lapply(
  seq(along=icecream$temp),
  function(i){
    y = ylim[1]:ylim[2]
    x = rep(icecream$temp[i], length(y))
    z0 = rep(0, length(y))
    z = dbinom(y, market.size, meanProb[i])
    return(list(x=x, y=y, z0=z0, z=z))
  }
)

glmModelPlot(x = icecream$temp, y=icecream$units,
             xlim=xlim, ylim=ylim,
             meanPred = meanPred, LwPred = LwPred,
             UpPred = UpPred, plotData = plotData,
             main = "Binomial (logit) GLM")
