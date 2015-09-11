library(arm)
library(caret)
library(pROC)
require(plyr)
require(gam)
library(glmnet)
data<-read.csv("sample_agn.csv",header=TRUE,na.strings="")
data2<-na.omit(data)
data2<-data2[data2$logMstar>0,]
data2<-data2[which(data2$vlos_sigma!=Inf),]

#data2$bpt<-as.factor(data2$bpt)
data2$bpt <- revalue(data2$bpt,c("SF"="0","Composite"="0",
                               "LINER"="1","Seyfert/LINER"="1",
                               "Seyfert"="1","BLANK"="0"))




#sSFR<-data$sSFR
#D1<-densityMclust(sSFR,G=2)
#plot(D1, what = "density", data = sSFR, breaks = 25)
#D1$classification

#data$SFR_class<-D1$classification-1

x<-as.matrix(data2[,2:5])
fit<-glmnet(x,y=data2$bpt,alpha=1,family="binomial")
plot(fit,xvar="lambda")
plot(fit, xvar = "dev", label = TRUE)
cv.glmmod <- cv.glmnet(x,y=data2$bpt,alpha=1,family="binomial",type.measure = "class")
plot(cv.glmmod)
best_lambda <- cv.glmmod$lambda.min

coef.min = coef(cv.glmmod, s = "lambda.min")
active.min = which(coef.min != 0)
index.min = coef.min[active.min]

fit<-gam(bpt~s(logMstar,3)+s(vlos_sigma,3),family=binomial,data=data2)


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
