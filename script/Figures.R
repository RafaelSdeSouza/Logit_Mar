# Plots for each parameter
#read data
require(ggplot2)
require(ggthemes)
require(Cairo)

gMx_S<-read.table("gMx_S.dat",header=TRUE)
gMx_E<-read.table("gMx_E.dat",header=TRUE)
#gMx_E$x<-gMx_S$x
gMx_S$x_o<-gMx_S$x*0.289949+13.86258
gMx_E$x_o<-gMx_E$x*0.289949+13.86258

gRx_S<-read.table("gRx_S.dat",header=TRUE)
gRx_E<-read.table("gRx_E.dat",header=TRUE)

gRx_S$x_o<-gRx_S$x*2.061402+2.180115
gRx_E$x_o<-gRx_E$x*2.061402+2.180115

#gRx_E$x<-gRx_S$x




#b) M_halo
PMx<-ggplot(aes(x=x_o,y=mean),data=gMx_S)+
  geom_ribbon(data=gMx_S,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.8,  fill=c("#E0FFFF")) +
  geom_ribbon(data=gMx_S,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.7,  fill=c("#00CED1")) +
  geom_line(size=1,linetype="dashed")+
   geom_ribbon(data=gMx_E,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2), alpha=0.8, fill=c("#fcbba1")) +
  geom_ribbon(data=gMx_E,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1), alpha=0.7, fill=c("#de2d26")) +
  geom_line(aes(x=x_o,y=mean),data=gMx_E,size=1,linetype="dotted")+
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
#  xlab(expression(log~M[halo]/M['\u0298']))+
  xlab(expression(log~M[halo]))+
  ylab(expression(P[AGN]))
cairo_pdf("..//figures/P_Mx.pdf",width = 9.25, height = 9)
PMx
dev.off()

#c) R_proj
PRx<-ggplot(aes(x=x_o,y=mean),data=gRx_S)+
  geom_ribbon(data=gRx_S,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.8,  fill=c("#E0FFFF")) +
  geom_ribbon(data=gRx_S,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.7,  fill=c("#00CED1")) +
  geom_line(size=1,linetype="dashed")+
  geom_ribbon(data=gRx_E,aes(x=x_o,y=mean,ymin=lwr2, ymax=upr2),alpha=0.8, fill=c("#fcbba1")) +
    geom_ribbon(data=gRx_E,aes(x=x_o,y=mean,ymin=lwr1, ymax=upr1),alpha=0.7, fill=c("#de2d26")) +
  geom_line(aes(x=x_o,y=mean),data=gRx_E,size=1,linetype="dashed")+
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(R/R[vir]))+ylab(expression(P[AGN]))
#+coord_cartesian(xlim=c(0,10))
CairoPDF("P_Rx.pdf",width = 9.25, height = 9)
PRx
dev.off()



# Plots beta posteriors

gplot_S<-read.table("gplot_S.dat",header=TRUE)
gplot_E<-read.table("gplot_E.dat",header=TRUE)
gplot_S$Parameter<-as.factor(gplot_S$Parameter)
gplot_E$Parameter<-as.factor(gplot_E$Parameter)


gplot<-rbind(gplot_S,gplot_E,deparse.level = 2)

gplot$gal<-as.factor(gplot$gal)
pL<-ggplot(data=gplot,aes(x=value,group=Parameter,fill=gal))+
  geom_density(colour="white",size=0.01,alpha=0.8)+facet_grid(Parameter~gal,labeller = label_parsed)+
  theme_hc()+
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
  scale_fill_manual(values=c("#00CED1","#de2d26"))+
  ylab("Density")+xlab("Parameter value")

CairoPDF("..//figures/beta_E.pdf",width = 8.5, height = 9)
pL
dev.off()





# Plot 3D
pi_AGN<-summary(as.mcmc.list(jags.logit, vars="px"))
pi_AGN<-pi_AGN$quantiles
#gdata1<-data.frame(x=data_cut$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])


z<-as.numeric(pi_AGN[,3])
x<-Mx
y<-Rx

library(rsm)
library(lattice)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#p<-persp(x,y,z, theta=150, phi=20, 
#         expand = 0.5,shade = 0.1,
#         xlab="Z", ylab=expression(NII.Ha), zlab=expression(log10.EW.Ha),ticktype='detailed',
#         col = YlOrBr,border=NA,xlog=T,ylog=T)
cor = topo.colors(200)

cairo_pdf("logit3D.pdf")
trellis.par.set("axis.line",list(axis.text=list(cex=20),col=NA,lty=1,lwd=2))
par(mar=c(1,1,1,1))
wireframe(z~x+y,data=data.frame(x=x, y=rep(y, each=length(x)), z=z),
          par.settings = list(regions=list(alpha=0.4)),
          col.regions =cor,drape=T,light.source = c(5,5,5),colorkey = FALSE,
          xlab=list(label=expression(M[halo]),cex=1.25), ylab=list(label=expression(R/R[vir]),cex=1.25), 
          zlab=list(rot=90,label=expression(P[AGN]),cex=1.25,dist=-1,rot=0),
          scale=list(tck=0.75,arrows=FALSE,distance =c(0.75, 0.75, 0.75)))


dev.off() 




#gdata_pall<-data.frame(x=data_n2$sfr_tot_p50,mean=pi_AGN[,3],lwr1=pi_AGN[,2],lwr2=pi_AGN[,1],upr1=pi_AGN[,4],upr2=pi_AGN[,5])
# Plot fit 

P1<-ggplot(aes(x=x,y=mean),data=gdata)+geom_line()+
  geom_ribbon(data=gdata,aes(x=x,y=mean,ymin=lwr1, ymax=upr1), alpha=0.45, fill=c("#005502")) +
  geom_ribbon(data=gdata,aes(x=x,y=mean,ymin=lwr2, ymax=upr2), alpha=0.35, fill=c("#3A5F0B")) +
  theme_hc()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(SFR))+ylab(expression(P[AGN]))

cairo_pdf("logit_AGN.pdf",width = 8, height = 7)
P1
dev.off()
