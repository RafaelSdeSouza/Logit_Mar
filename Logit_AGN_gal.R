# Libraries
library(arm)
library(caret)
library(pROC)
require(plyr)
require(gam)
library(glmnet)
require(mgcv)
require(reshape2)
require(ggthemes)
require(Cairo)
data<-read.csv("sample_agn.csv",header=TRUE,na.strings="")
data2<-na.omit(data)
data2<-data2[data2$logMstar>0,]
data2<-data2[which(data2$vlos_sigma!=Inf),]

#data2$bpt<-as.factor(data2$bpt)
data2$bpt <- revalue(data2$bpt,c("SF"="0","Composite"="0",
                               "LINER"="1","Seyfert/LINER"="1",
                               "Seyfert"="1","BLANK"="0"))


# Scale variables
#data2$logMstar<-(data2$logMstar-mean(data2$logMstar))/sd(data2$logMstar)
#data2$logMhalo<-(data2$logMhalo-mean(data2$logMhalo))/sd(data2$logMhalo)
#data2$vlos_sigma<-(data2$vlos_sigma-mean(data2$vlos_sigma))/sd(data2$vlos_sigma)
#data2$r_rvir<-(data2$r_rvir-mean(data2$r_rvir))/sd(data2$r_rvir)

# Variable selection via LASSO

x<-as.matrix(data2[,2:5])
fit<-glmnet(x,y=data2$bpt,alpha=1,family="binomial")
predict(fit, type="coefficients")

plot(fit,xvar="lambda",label = TRUE)
plot(fit, xvar = "dev", label = TRUE)

# Cross-validation

cv.glmmod <- cv.glmnet(x,y=data2$bpt,alpha=1,family="binomial",type.measure = "auc")
plot(cv.glmmod)
best_lambda <- cv.glmmod$lambda.min

coef.min = coef(cv.glmmod, s = "lambda.min")
active.min = coef.min[which(abs(coef.min) > 0.05 )]






# Fit and plot remained variables

fit=glm(bpt~logMstar+r_rvir,data=data2,family=binomial("logit"))
fit2<-gam(bpt~s(logMstar,3)+s(vlos_sigma,3),family=binomial,data=data2)


# Diagnostic 

ROCF<- data.frame(True=data2$bpt,predicted=predict(fit,type = "response"))
F1 <-roc(ROCF$True,ROCF$predicted)
coords(F1,x="best")[1]

ROCF$class<-ROCF$predicted
ROCF$class[which(ROCF$class>=coords(F1,x="best")[1])]<-1
ROCF$class[which(ROCF$class<coords(F1,x="best")[1])]<-0

plot(F1)
confusionMatrix(ROCF$True, ROCF$class)


# Plot 

x <-range(data2$logMstar)
x <- seq(0.95*x[1], 1.25*x[2], length.out=50)
y <- range(data2$r_rvir)
y <- seq(y[1], 1.25*y[2], length.out=50)

z <- outer(x,y,
           function(logMstar,r_rvir)
             predict(fit, data.frame(logMstar,r_rvir),type = 'response'))
library(rsm)
library(lattice)
YlOrBr <- c("#00A3DB")
#p<-persp(x,y,z, theta=150, phi=20,
#         expand = 0.5,shade = 0.1,
#         xlab="Z", ylab=expression(NII.Ha), zlab=expression(log10.EW.Ha),ticktype='detailed',
#         col = YlOrBr,border=NA,xlog=T,ylog=T)
cairo_pdf("logit3D.pdf")
trellis.par.set("axis.line",list(axis.text=list(cex=20),col=NA,lty=1,lwd=2))
par(mar=c(1,1,1,1))
wireframe(z~x+y,data=data.frame(x=x, y=rep(y, each=length(x)), z=z),
          par.settings = list(regions=list(alpha=0.6)),
          col.regions =YlOrBr,drape=T,light.source = c(5,5,5),colorkey = FALSE,
          xlab=list(label=expression(log~M[star]/M['\u0298']),cex=1.25),
          ylab=list(label=expression(R[vir](R['\u0298'])),cex=1.25),
          zlab=list(rot=90,label=expression(P[AGN]),cex=1.25,dist=-1,rot=0),
          scale=list(tck=0.75,arrows=FALSE,distance =c(0.75, 0.75, 0.75)))

dev.off()


g.dat<-melt(data2,id="bpt")
mf_labeller <- function(var, value){
  value <- as.character(value)
  if (var=="sex") { 
    value[value=="Female"] <- "Woman"
    value[value=="Male"]   <- "Man"
  }
  return(value)
}
g.dat$bpt<-as.factor(g.dat$bpt)
g.dat$bpt<-revalue(g.dat$bpt,c("0"="No AGN","1"="AGN"))


facet_wrap_labeller <- function(gg.plot,labels=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  require(gridExtra)
  
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  
  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }
  
  g$grobs <- gg
  class(g) = c("arrange", "ggplot",class(g)) 
  g
}
give.n <- function(x){
  
  return(c(y = 0.5, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}




g1<-ggplot(aes(y=value,x=bpt),data=g.dat)+geom_boxplot(aes(fill=bpt))+
  facet_wrap(~variable,ncol=2,scales="free")+theme_hc()+
  scale_fill_economist()+ylab("")+xlab("")+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=20),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))
g2<-facet_wrap_labeller(g1,labels=c(expression(log~M[star]/M['\u0298']),expression(log~M[h]/M['\u0298']),
                                    expression(sigma~(km/s)),expression(R[vir](R['\u0298']))))

CairoPDF("box.pdf",width = 9,height = 8)
g2
dev.off()