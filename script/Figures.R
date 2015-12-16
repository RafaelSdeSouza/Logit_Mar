# Plots for each parameter
#read data
require(ggplot2)
require(ggthemes)
require(Cairo)
require(plyr)

gMx_S<-read.table("..//data/gMx_S.dat",header=TRUE)
gMx_E<-read.table("..//data/gMx_E.dat",header=TRUE)
#gMx_E$x<-gMx_S$x
gMx_S$x_o<-gMx_S$x*0.289949+13.86258
gMx_E$x_o<-gMx_E$x*0.289949+13.86258

gRx_S<-read.table("..//data/gRx_S.dat",header=TRUE)
gRx_E<-read.table("..//data/gRx_E.dat",header=TRUE)

gRx_S$x_o<-gRx_S$x*2.061402+2.180115
gRx_E$x_o<-gRx_E$x*2.061402+2.180115




#gRx_E$x<-gRx_S$x




#b) M_halo
#PMx<-
  ggplot(aes(x=x_o,y=mean),data=gMx_S)+
    geom_ribbon(data=gMx_E,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2), alpha=0.9, fill=c("#fcbba1")) +
    geom_ribbon(data=gMx_E,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1), alpha=0.8, fill=c("#de2d26")) +
    geom_line(aes(x=x_o,y=mean),data=gMx_E,size=1,linetype="dotted")+
  geom_ribbon(data=gMx_S,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.6,  fill=c("#E0FFFF")) +
  geom_ribbon(data=gMx_S,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.5,  fill=c("#00CED1")) +
  geom_line(size=1,linetype="dashed")+
   
  theme_bw()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
#  xlab(expression(log~M[halo]/M['\u0298']))+
  xlab(expression(log~M[halo]))+
  ylab(expression(P[AGN]))+coord_cartesian(ylim=c(0.15,0.85))
quartz.save(type = 'pdf', file = '..//figures/P_Mx.pdf',width = 9.5, height = 9)

#cairo_pdf("..//figures/P_Mx.pdf",width = 9.25, height = 9)
#PMx
#dev.off()

#c) R_proj
#PRx<-
  ggplot(aes(x=x_o,y=mean),data=gRx_S)+
    geom_ribbon(data=gRx_E,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.9, fill=c("#fcbba1")) +
    geom_ribbon(data=gRx_E,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.8, fill=c("#de2d26")) +
    geom_line(aes(x=x_o,y=mean),data=gRx_E,size=1,linetype="dotted")+
  geom_ribbon(data=gRx_S,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.6,  fill=c("#E0FFFF")) +
  geom_ribbon(data=gRx_S,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.5,  fill=c("#00CED1")) +
  geom_line(size=1,linetype="dashed")+
   theme_bw()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(R/R[vir]))+ylab(expression(P[AGN]))+coord_cartesian(ylim=c(0.15,0.85))
#+coord_cartesian(xlim=c(0,10))

quartz.save(type = 'pdf', file = '..//figures/P_Rx.pdf',width = 9.5, height = 9)
#CairoPDF("..//figures/P_Rx.pdf",width = 9.25, height = 9)
#PRx
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




