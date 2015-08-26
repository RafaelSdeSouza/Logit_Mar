library(mixtools)
library(mclust)
library(arm)
library(caret)
library(pROC)
data<-read.csv("sample.csv",header=TRUE)


sSFR<-data$sSFR
D1<-densityMclust(sSFR,G=2)
plot(D1, what = "density", data = sSFR, breaks = 25)
D1$classification

data$SFR_class<-D1$classification-1


fit<-bayesglm(SFR_class~r_rvir,family=binomial,data=data)

library(popbio)
logi.hist.plot(data$r_rvir,data$SFR_class,boxp=FALSE,type="hist",col="gray")



ROCF<- data.frame(True=data$SFR_class,predicted=predict(fit,type = "response"))
F1 <-roc(ROCF$True,ROCF$predicted)
coords(F1,x="best")[1]

ROCF$class<-ROCF$predicted
ROCF$class[which(ROCF$class>=coords(F1,x="best")[1])]<-1
ROCF$class[which(ROCF$class<coords(F1,x="best")[1])]<-0

plot(F1)
confusionMatrix(ROCF$True, ROCF$class)




x <-range(data$logMstar)
x <- seq(x[1], x[2], length.out=50)    
y <- range(data$logMhalo)
y <- seq(y[1], y[2], length.out=50)

z <- outer(x,y, 
           function(logMstar,logMhalo)
             predict(fit, data.frame(logMstar,logMhalo),type = 'response'))
library(rsm)
library(lattice)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#p<-persp(x,y,z, theta=150, phi=20, 
#         expand = 0.5,shade = 0.1,
#         xlab="Z", ylab=expression(NII.Ha), zlab=expression(log10.EW.Ha),ticktype='detailed',
#         col = YlOrBr,border=NA,xlog=T,ylog=T)
#cor = topo.colors(200)

cairo_pdf("logit3D.pdf")
trellis.par.set("axis.line",list(axis.text=list(cex=20),col=NA,lty=1,lwd=2))
par(mar=c(1,1,1,1))
wireframe(z~x+y,data=data.frame(x=x, y=rep(y, each=length(x)), z=z),
          par.settings = list(regions=list(alpha=0.4)),
          col.regions =cor,drape=T,light.source = c(5,5,5),colorkey = FALSE,
          xlab=list(label=expression(log10.NII.Ha.),cex=1.25), ylab=list(label=expression(log10.EW.Ha..),cex=1.25), 
          zlab=list(rot=90,label=expression(pi),cex=1.25,dist=-1,rot=0),
          scale=list(tck=0.75,arrows=FALSE,distance =c(0.75, 0.75, 0.75)))


dev.off() 
