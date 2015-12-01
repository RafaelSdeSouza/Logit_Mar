require(rstan)
require(lme4)
require(glmer2stan)
require(sjmisc)
require(binomTools)
require(sjPlot)
require(ggplot2)
require(ggthemes)
require(gam)
data     <- read.table("..//data/matched.txt",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]
data_cut   <- subset(data_cut, zoo == "E")


mod1<- glm(bpt ~ logM200_L+RprojLW_Rvir, family=binomial("logit"),data=data_cut)
mod2<- gam(bpt ~ s(RprojLW_Rvir)+s(logM200_L), family=binomial("logit"),data=data_cut)


plot(fitted(mod2), jitter(data_cut$bpt, amount=0.05), xlab="Fitted values",
                      ylab="Probability of presence", las=1, cex.lab=1.2, cex=0.8)
abline(0,1, lty=3)
t.breaks <-cut(fitted(mod2), seq(0,1, by=0.1))
means <-tapply(data_cut$bpt, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
means.se <-tapply(data_cut$bpt, t.breaks, semean)
points(seq(0.05, 0.95, by=0.1), means, pch=16, col="orange")
segments(seq(0.05, 0.95, by=0.1), means-2*means.se,
         seq(0.05, 0.95,by=0.1), means+2*means.se,lwd=2, col="orange")


gdata<-data.frame(fit=fitted(mod),obs=data_cut$bpt)
xrange<-range(fitted(mod))

bined<-data.frame(x=seq(0.05, 0.95, by=0.1),y=means)


ggplot(aes(x=fit,y=obs),data=gdata)+geom_point(colour="#00CED1",size=3,position = position_jitter (h = 0.05))+
  geom_point(aes(x=x,y=y),size=3,data=bined,colour="#de2d26")+
  geom_errorbar(data=bined,guide="none",aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),alpha=0.7,
                colour="#de2d26",width=0.005)+
  geom_abline(size=1.5,intercept = 0, slope = 1,linetype="dashed",colour="gray65")+coord_cartesian(xlim=c(0.75*xrange[1],1.1*xrange[2]))+
                 theme_bw()+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                                  axis.title.y=element_text(vjust=0.75),
                                  axis.title.x=element_text(vjust=-0.25),
                                  text = element_text(size=25))+xlab("Predicted  AGN fraction")+
  ylab("Observed AGN Fraction")






## summary(budworm.lg)

(Rsq.budworm <- Rsq(fit))

plot(Rsq.budworm, "ROC")



fit <- glm(bpt ~ logM200_L+RprojLW_Rvir ,
             data_cut,
             family = binomial("logit"))
sjp.glm(fit)

sjp.glm(fit, type = "qq")

m1_lme4 <- lmer( Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE )

# construct subject index --- glmer2stan forces you to manage your own cluster indices
sleepstudy$subject_index <- as.integer(as.factor(sleepstudy$Subject))

# fit with glmer2stan
m1_g2s <- lmer2stan( Reaction ~ Days + (Days | subject_index), data=sleepstudy )
m1_g2s@stanmodel