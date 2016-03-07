# Plots for each parameter
#read data
require(ggplot2)
require(ggthemes)
require(Cairo)
require(plyr)

##------------------------------------------------------------------
gMx_S<-read.table("..//data/gMx_S.dat",header=TRUE)
gMx_E<-read.table("..//data/gMx_E.dat",header=TRUE)
#gMx_E$x<-gMx_S$x
gMx_S$x_o<-gMx_S$x*0.289949+13.86258
gMx_E$x_o<-gMx_E$x*0.289949+13.86258

gRx_S<-read.table("..//data/gRx_S.dat",header=TRUE)
gRx_E<-read.table("..//data/gRx_E.dat",header=TRUE)

gRx_S$x_o<-gRx_S$x*2.061402+2.180115
gRx_E$x_o<-gRx_E$x*2.061402+2.180115



# Prepare points ------------------------------------------------------------------------
dataE     <- read.table("..//data/matched_E_original.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S_original.txt",header=TRUE,na.strings="")

MmeanE <-mean(dataE$logM200_L)
MsdE <-sd(dataE$logM200_L)
MmeanS <-mean(dataS$logM200_L)
MsdS <-sd(dataS$logM200_L)

dataE_R<-dataE[which(dataE$logM200_L<=MmeanE+MsdE & dataE$logM200_L >= MmeanE - MsdE),]
dataS_R<-dataS[which(dataS$logM200_L<=MmeanS+MsdS & dataS$logM200_L >= MmeanS - MsdS),]

t.breaks <-cut(dataE_R$RprojLW_Rvir, breaks=c(0,1.5,3,4.5,6,8))
means <-tapply(dataE_R$bpt, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
#semean <-function(x) sqrt(x*(1-x)/length(x))
means.se <-tapply(dataE_R$bpt, t.breaks, semean)
bins<-levels(t.breaks)


gdata<-data.frame(x=bins,y=means)
gdata$gal<-rep("E",nrow(gdata))


t.breaks2 <-cut(dataS_R$RprojLW_Rvir, breaks=c(0,1.5,3,4.5,6,8))
means2 <-tapply(dataS_R$bpt, t.breaks2, mean)
semean2 <-function(x) sd(x)/sqrt(length(x))
#semean2 <-function(x) sqrt(x*(1-x)/length(x))
means.se2 <-tapply(dataS_R$bpt, t.breaks2, semean2)

gdata2<-data.frame(x=bins,y=means2)
gdata2$gal<-rep("S",nrow(gdata2))

gdata$xc<-c(0.75,2.25,3.75,5.25,7)
gdata2$xc<-c(0.8,2.5,4,5.5,7.25)






#gRx_E$x<-gRx_S$x





#c) R_proj
#PRx<-
  ggplot(aes(x=x_o,y=mean),data=gRx_S)+
    geom_ribbon(data=gRx_E,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.9, fill=c("#fcbba1")) +
    geom_ribbon(data=gRx_E,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.8, fill=c("#de2d26")) +
    geom_line(aes(x=x_o,y=mean),data=gRx_E,size=1,linetype="dotted")+
  geom_ribbon(data=gRx_S,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.6,  fill=c("#E0FFFF")) +
  geom_ribbon(data=gRx_S,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.5,  fill=c("#00CED1")) +
  geom_line(size=1,linetype="dashed")+
    geom_point(aes(x=xc,y=y),size=3,data=gdata,colour="red3")+
    geom_errorbar(data=gdata,guide="none",aes(x=xc,y=y,ymin=y-2*means.se,ymax=y+2*means.se),
                  colour="red3",width=0.05)+
    geom_point(aes(x=xc,y=y),size=3,data=gdata2,colour="cyan3")+
    geom_errorbar(data=gdata2,guide="none",aes(x=xc,y=y,ymin=y-1.96*means.se2,ymax=y+1.96*means.se2),
                  colour="cyan3",width=0.05)+
   theme_bw()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(R/R[vir]))+ylab(expression(P[AGN]))+coord_cartesian(ylim=c(0.1,1))
#  +coord_cartesian(xlim=c(0,10))

  

    
    
quartz.save(type = 'pdf', file = '..//figures/P_Rx.pdf',width = 9.5, height = 9)
#CairoPDF("..//figures/P_Rx.pdf",width = 9.25, height = 9)
#PRx
#dev.off()




##------------------------------------------------------------------
RmeanE <-mean(dataE$RprojLW_Rvir)
RsdE <-sd(dataE$RprojLW_Rvir)
RmeanS <-mean(dataS$RprojLW_Rvir)
RsdS <-sd(dataS$RprojLW_Rvir)

dataE_M<-dataE[which(dataE$RprojLW_Rvir<=RmeanE+RsdE & dataE$RprojLW_Rvir >= RmeanE - RsdE),]
dataS_M<-dataS[which(dataS$RprojLW_Rvir<=RmeanS+RsdS & dataS$RprojLW_Rvir >= RmeanS - RsdS),]

t.breaks <-cut(dataE_M$logM200_L, breaks=c(13,13.3,13.6,13.9,14.2,14.6))
means <-tapply(dataE_M$bpt, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
#semean <-function(x) sqrt(x*(1-x)/length(x))
means.se <-tapply(dataE_M$bpt, t.breaks, semean)
bins<-levels(t.breaks)


gdata<-data.frame(x=bins,y=means)
gdata$gal<-rep("E",nrow(gdata))


t.breaks2 <-cut(dataS_M$logM200_L, breaks=c(13,13.3,13.6,13.9,14.2,14.6))
means2 <-tapply(dataS_M$bpt, t.breaks2, mean)
semean2 <-function(x) sd(x)/sqrt(length(x))
#semean2 <-function(x) sqrt(x*(1-x)/length(x))
means.se2 <-tapply(dataS_M$bpt, t.breaks2, semean2)

gdata2<-data.frame(x=bins,y=means2)
gdata2$gal<-rep("S",nrow(gdata2))

gdata$xc<-c(13.2,13.5,13.8,14.1,14.4)
gdata2$xc<-c(13.225,13.525,13.825,14.125,14.425)



##------------------------------------------------------------------
#b) M_halo
#PMx<-
ggplot(aes(x=x_o,y=mean),data=gMx_S)+
  geom_ribbon(data=gMx_E,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2), alpha=0.9, fill=c("#fcbba1")) +
  geom_ribbon(data=gMx_E,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1), alpha=0.8, fill=c("#de2d26")) +
  geom_line(aes(x=x_o,y=mean),data=gMx_E,size=1,linetype="dotted")+
  geom_ribbon(data=gMx_S,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.6,  fill=c("#E0FFFF")) +
  geom_ribbon(data=gMx_S,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.5,  fill=c("#00CED1")) +
  geom_line(size=1,linetype="dashed")+
  geom_point(aes(x=xc,y=y),size=3,data=gdata,colour="red3")+
  geom_errorbar(data=gdata,guide="none",aes(x=xc,y=y,ymin=y-1.96*means.se,ymax=y+1.96*means.se),
                colour="red3",width=0.05)+
  geom_point(aes(x=xc,y=y),size=3,data=gdata2,colour="cyan3")+
  geom_errorbar(data=gdata2,guide="none",aes(x=xc,y=y,ymin=abs(y-1.96*means.se2),ymax=y+1.96*means.se2),
                colour="cyan3",width=0.05)+
  theme_bw()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  #  xlab(expression(log~M[halo]/M['\u0298']))+
  xlab(expression(log~M[halo]))+
  ylab(expression(P[AGN]))+coord_cartesian(xlim=c(13.3,14.7),ylim=c(0.025,1))


quartz.save(type = 'pdf', file = '..//figures/P_Mx.pdf',width = 9.5, height = 9)

#cairo_pdf("..//figures/P_Mx.pdf",width = 9.25, height = 9)
#PMx
#dev.off()




# Plots beta posteriors

#gplot_S<-read.table("gplot_S.dat",header=TRUE)
#gplot_E<-read.table("gplot_E.dat",header=TRUE)
#gplot_S$Parameter<-as.factor(gplot_S$Parameter)
#gplot_E$Parameter<-as.factor(gplot_E$Parameter)


#gplot<-rbind(gplot_S,gplot_E,deparse.level = 2)
gplot<-read.table("gplot.dat",header=TRUE)

gplot$gal<-as.factor(gplot$gal)

gplot$Parameter<-revalue(gplot$Parameter, c("beta[1,1]"= "beta[1]","beta[1,2]"="beta[1]", "beta[2,1]" ="beta[2]", "beta[2,2]"="beta[2]", "beta[3,1]"="beta[3]",
             "beta[3,2]"="beta[3]"))


pL<-ggplot(data=gplot,aes(x=value,group=gal,fill=gal))+
  geom_density(colour="white",size=0.01,alpha=0.8)+facet_grid(Parameter~gal,labeller = label_parsed)+
  theme_bw()+
  theme(legend.position="none",panel.background = element_rect(fill = "white"),plot.background = element_rect(
    fill = "white"),plot.title = element_text(hjust=0.5),
    axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=22),
    axis.text.y=element_text(size=22),
    strip.text.x=element_text(size=25),
    axis.title.x=element_text(vjust=-0.25),
    text = element_text(size=22),axis.title.x=element_text(size=rel(1)),strip.background=element_rect(
      fill = "white"),strip.text=element_text(
        size = 25))+
  geom_vline(xintercept=0,size=1,linetype="dashed",colour=c("grey50")) +
  #  scale_fill_manual(values=c("#E0FFFF","#00CED1","cyan4"))+
  scale_fill_manual(values=c("#de2d26","#00CED1"))+
  ylab("Density")+xlab("Parameter value")


quartz.save(type = 'pdf', file = '..//figures/betas.pdf',width = 8.5, height = 9)

#CairoPDF("..//figures/betas.pdf",width = 8.5, height = 9)
#pL
#dev.off()




