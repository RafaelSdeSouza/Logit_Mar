require(binomTools)
require(ggplot2)
require(ggthemes)
require(arm)
require(blmeco)
dataE     <- read.table("..//data/matched_E.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S.txt",header=TRUE,na.strings="")
data<-rbind(dataE,dataS)
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]




data_cut_E   <- subset(data_cut, zoo == "E")
mod1<- glm(bpt ~ logM200_L+RprojLW_Rvir, family=binomial(link = "logit"),data=data_cut_E)
WAIC<-WAIC(mod1)$WAIC2  
WAIC0<-WAIC(glm(bpt ~ 1, family=binomial(link = "logit"),data=data_cut_E))$WAIC2

# GoF Visual 
binagem<-0.1

t.breaks <-cut(fitted(mod1), seq(0,1, by=binagem))
means <-tapply(data_cut_E$bpt, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
means.se <-tapply(data_cut_E$bpt, t.breaks, semean)
y<-data_cut_E$bpt
#points(seq(0.05, 0.95, by=0.1), means, pch=16, col="orange")
#segments(seq(0.05, 0.95, by=0.1), means-2*means.se,
#         seq(0.05, 0.95,by=0.1), means+2*means.se,lwd=2, col="orange")
#semean <-function(x) sd(x)/sqrt(length(x))
gdata<-data.frame(fit=fitted(mod1),obs=y)
xrange<-range(fitted(mod1))
bined<-data.frame(x=seq(binagem,1, by=binagem),y=means)
ggplot(aes(x=fit,y=obs),data=gdata)+geom_point(colour="#00CED1",size=3,position = position_jitter (h = 0.025))+
  geom_point(aes(x=x,y=y),size=3,data=bined,colour="#de2d26")+
  geom_errorbar(data=bined,guide="none",aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),alpha=0.7,
                colour="#de2d26",width=0.005)+
  geom_abline(size=1.5,intercept = 0, slope = 1,linetype="dashed",colour="gray65")+coord_cartesian(xlim=c(0.75*xrange[1],1.1*xrange[2]))+
                 theme_bw()+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                                  axis.title.y=element_text(vjust=0.75),
                                  axis.title.x=element_text(vjust=-0.25),
                                  text = element_text(size=25))+xlab("Predicted  AGN fraction")+
  ylab("Observed AGN Fraction")+
  annotate("text",color="#de2d26",size=5,label=paste("This model: WAIC = ",round(WAIC,1),sep=""),x=0.65,y=0.2,hjust = 0)+
  annotate("text",color="gray60",size=5,label=paste("Intercept only: WAIC = ",round(WAIC0,1),sep=""),x=0.65,y=0.25,hjust = 0)+
  coord_cartesian(xlim=c(0.25,0.85),ylim=c(-0.1,1.1))


quartz.save(type = 'pdf', file = '..//figures/GoF_E.pdf',width = 9, height = 8)




data_cut_S   <- subset(data_cut, zoo == "S")
mod2<- glm(bpt ~ logM200_L+RprojLW_Rvir, family=binomial(link = "logit"),data=data_cut_S)
WAIC<-WAIC(mod2)$WAIC2  
WAIC0<-WAIC(glm(bpt ~ 1, family=binomial(link = "logit"),data=data_cut_S))$WAIC2 

# GoF Visual 
binagem2<-0.005
t.breaks <-cut(fitted(mod2), seq(0,1, by=binagem2))
means <-tapply(data_cut_S$bpt, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
means.se <-tapply(data_cut_S$bpt, t.breaks, semean)
y<-data_cut_S$bpt
#points(seq(0.05, 0.95, by=0.1), means, pch=16, col="orange")
#segments(seq(0.05, 0.95, by=0.1), means-2*means.se,
#         seq(0.05, 0.95,by=0.1), means+2*means.se,lwd=2, col="orange")
#semean <-function(x) sd(x)/sqrt(length(x))
gdata<-data.frame(fit=fitted(mod2),obs=y)
xrange<-range(fitted(mod2))
bined<-data.frame(x=seq(binagem2,1, by=binagem2),y=means)
ggplot(aes(x=fit,y=obs),data=gdata)+geom_point(colour="#00CED1",size=3,position = position_jitter (h = 0.025))+
  geom_point(aes(x=x,y=y),size=3,data=bined,colour="#de2d26")+
  geom_errorbar(data=bined,guide="none",aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),alpha=0.7,
                colour="#de2d26",width=0.005)+
  geom_abline(size=1.5,intercept = 0, slope = 1,linetype="dashed",colour="gray65")+coord_cartesian(xlim=c(0.75*xrange[1],1.1*xrange[2]))+
  theme_bw()+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                   axis.title.y=element_text(vjust=0.75),
                   axis.title.x=element_text(vjust=-0.25),
                   text = element_text(size=25))+xlab("Predicted  AGN fraction")+
  ylab("Observed AGN Fraction")+
  annotate("text",color="gray60",size=5,label=paste("This model: WAIC = ",round(WAIC,1),sep=""),x=0.505,y=0.2,hjust = 0)+
  annotate("text",color="#de2d26",size=5,label=paste("Intercept only: WAIC = ",round(WAIC0,1),sep=""),x=0.505,y=0.25,hjust = 0)+
  coord_cartesian(xlim=c(0.479,0.519),ylim=c(-0.1,1.1))


quartz.save(type = 'pdf', file = '..//figures/GoF_S.pdf',width = 9, height = 8)


# ROC curve 

source("ROC.R")
ROCtest(mod1)


require(LogisticDx)
GOF<-gof(mod1)

residuals(mod1, "pearson") 

E1 <- resid(mod1,type="pearson")
sum(E1^2)/(mod1$df.residual)

mod0<- bayesglm(bpt ~ 1, family=binomial(link = "logit"),data=data_cut)
mod1 <-glm(bpt ~ logM200_L+RprojLW_Rvir, family=binomial(link = "logit"),data=data_cut)
PseudoR2(mod1)
PseudoR2(mod0)



plot(predict(mod1),residuals(mod1),col=c("blue","red")[1+data_cut$bpt])
 abline(h=0,lty=2,col="grey")
 lines(lowess(predict(mod1),residuals(mod1)),col="black",lwd=2)

 rl=lm(residuals(mod1)~bs(predict(mod1),8))
 #rl=loess(residuals(reg)~predict(reg))
   y=predict(rl,se=TRUE)
 segments(predict(mod1),y$fit+2*y$se.fit,predict(mod1),y$fit-2*y$se.fit,col="green")


#mod2<- gam(bpt ~ s(RprojLW_Rvir,3)+s(logM200_L,3), family=binomial("logit"),data=data_cut)
#plot(mod2,pages=1,residuals=TRUE)
#t.breaks <-cut(fitted(mod2), seq(0,1, by=0.1))
#means <-tapply(data_cut$bpt, t.breaks, mean)
#semean <-function(x) sd(x)/sqrt(length(x))
#means.se <-tapply(data_cut$bpt, t.breaks, semean)
#y<-data_cut$bpt
#gdata2<-data.frame(fit=fitted(mod2),obs=y)
#xrange<-range(fitted(mod2))
#bined<-data.frame(x=seq(0.05, 0.95, by=0.1),y=means)
#ggplot(aes(x=fit,y=obs),data=gdata2)+geom_point(colour="#00CED1",size=3,position = position_jitter (height = 0.05))+
#  geom_point(aes(x=x,y=y),size=3,data=bined,colour="#de2d26")+
#  geom_errorbar(data=bined,guide="none",aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),alpha=0.7,
#                colour="#de2d26",width=0.005)+
#  geom_abline(size=1.5,intercept = 0, slope = 1,linetype="dashed",colour="gray65")+coord_cartesian(xlim=c(0.75*xrange[1],1.1*xrange[2]))+
#  theme_bw()+theme(legend.position="top",plot.title = element_text(hjust=0.5),
#                   axis.title.y=element_text(vjust=0.75),
#                   axis.title.x=element_text(vjust=-0.25),
#                   text = element_text(size=25))+xlab("Predicted  AGN fraction")+
#  ylab("Observed AGN Fraction")


#AIC


#BIC




