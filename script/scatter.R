#scatter 3D

require(plyr)
require(reshape2)
require(ggplot2)
require(ggthemes)

# Original--------------------
# Read and format data
data     <- read.csv("..//data/sample_CRP02_sub.csv",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","lgm_tot_p50","logM200_L","RprojLW_Rvir","sfr_tot_p50","color_gr","zoo")]
data2    <- na.omit(data_cut)
data2    <- data2[data2$lgm_tot_p50>0,]
data2    <- data2[which(data2$logM200_L>0),]
data2    <- data2[which(data2$RprojLW_Rvir>=0),]
data2    <- data2[which(data2$sfr_tot_p50>=-100),]
data3    <- data2
data3    <- subset(data3, bpt!="LINER")    # remove Liners
data3    <- subset(data3, bpt!="BLANK")    # remove BLANK
data3    <- subset(data3, bpt!="Composite")    # remove BLANK
data3$bpt <-droplevels(data3$bpt)
#data3$bpt  <- revalue(data3$bpt,c("Star Forming"="0","Composite"="0",
#                                  "LINER"="1","Seyfert/LINER"="1","Star Fo"="0",
#                                  "Seyfert"="1","BLANK"="0"))
data3$bpt  <- revalue(data3$bpt,c("Star Forming"="0","Seyfert/LINER"="1",
                                  "Seyfert"="1"))
data_n     <- data3
#data_n2    <- subset(data_n, zoo=="E" | zoo == "S")
dataEoriginal<-subset(data_n, zoo=="E")
dataSoriginal<-subset(data_n, zoo == "S")
dataEoriginal$type<-rep("E_full",nrow(dataEoriginal))
dataSoriginal$type<-rep("S_full",nrow(dataSoriginal))

dataEoriginal$PSM<-rep("Before",nrow(dataEoriginal))
dataSoriginal$PSM<-rep("Before",nrow(dataSoriginal))

dataoriginal<-rbind(dataEoriginal,dataSoriginal)


dataraw<-dataoriginal[,c("bpt" ,"lgm_tot_p50","sfr_tot_p50","color_gr","zoo","type","PSM")]
#-------------------------



dataE     <- read.table("..//data/matched_E_original.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S_original.txt",header=TRUE,na.strings="")
dataE$type<-rep("E_cut",nrow(dataE))
dataS$type<-rep("S_cut",nrow(dataS))

dataS$PSM<-rep("After",nrow(dataS))
dataE$PSM<-rep("After",nrow(dataE))

datam<-rbind(dataE,dataS)
datam<-datam[,c("bpt" ,"lgm_tot_p50","sfr_tot_p50","color_gr","zoo","type","PSM")]


Full_data<-rbind(dataraw,datam)




Full_data$zoo<-droplevels(Full_data$zoo)
Full_data$zoo<-revalue(Full_data$zoo, c("E" = "Ellipticals", "S" = "Spirals"))
ggplot(Full_data,aes(y=color_gr,x=lgm_tot_p50,size=bpt,shape=bpt,linetype=bpt,colour=bpt
))+
  geom_point(alpha=0.9,aes(colour=bpt))+
    scale_shape_manual(values=c(16,4))+
  scale_color_manual(name="",values=c("orange2","gray"))+
  stat_density2d( alpha=0.4,size=0.35,colour="black")+
  
  scale_size_manual(name="",values=c(2.5,0.75))+

  scale_linetype_manual(name="",values=c("solid","blank"))+
#  scale_linetype_stata()+
  #  scale_color_gradient()+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  coord_cartesian(ylim=c(0.1,1.25),xlim=c(9.5,12))+xlab(expression(log~M["*"]))+ylab("g-r")+
  facet_wrap(~zoo)

quartz.save(type = 'pdf', file = '..//figures/grMgal.pdf',width = 16, height = 8)


ggplot(Full_data,aes(y=color_gr,x=sfr_tot_p50,size=bpt,shape=bpt,linetype=bpt,colour=bpt
))+
  geom_point(alpha=0.9,aes(colour=bpt))+
  scale_shape_manual(values=c(16,4))+
  scale_color_manual(name="",values=c("orange2","gray"))+
  stat_density2d( alpha=0.4,size=0.35,colour="black")+
  
  scale_size_manual(name="",values=c(2.5,0.75))+
  
  scale_linetype_manual(name="",values=c("solid","blank"))+
  #  scale_linetype_stata()+
  #  scale_color_gradient()+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  coord_cartesian(ylim=c(0.15,1.25),xlim=c(-3,1.5))+
  xlab("log SFR")+ylab("g-r")+
  facet_wrap(~zoo)
quartz.save(type = 'pdf', file = '..//figures/grSFR.pdf',width = 16, height = 8)


ggplot(Full_data,aes(x=color_gr))+
  geom_density(alpha=0.9)+
  xlab("g-r")+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))
  facet_wrap(~zoo)

library(mixtools)
  g_r = Full_data$color_gr
  mixmdl = normalmixEM(g_r)
  plot(mixmdl,which=2,xlim=c(0,2),ylim=c(0,3))
  lines(density(g_r), lty=2, lwd=2)
  
  
Full_data$type<-as.factor(Full_data$type)
ggplot(Full_data,aes(y=lgm_tot_p50,x=sfr_tot_p50,color=type,linetype=type))+
  geom_point(size=0.25,alpha=0.2,shape=3)+
  geom_density2d(alpha=0.7,size=1)+
   scale_color_manual(name="",values=c("red3","coral", "royalblue3","green4"))+
  scale_linetype_manual(name="",values=c("solid", "dashed","solid","dashed"))+
  theme_bw()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  coord_cartesian(ylim=c(9.25,12.25))+xlab(expression(log~SFR~(M['\u0298']*yr^{-1})))+ylab(expression(log~M["*"]~(M['\u0298'])))


quartz.save(type = 'pdf', file = '..//figures/MgalSFR.pdf',width = 9.5, height = 9)




ggplot(Full_data,aes(y=color_gr,x=sfr_tot_p50,color=type,linetype=type,shape=PSM))+
  geom_point(size=0.25,alpha=0.2,shape=3)+
  geom_density2d(alpha=0.7,size=1)+
  scale_color_manual(name="",values=c("red3","coral", "royalblue3","green4"))+
  scale_linetype_manual(name="",values=c("solid", "dashed","solid","dashed"))+
  theme_bw()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  coord_cartesian(ylim=c(0,1.5))+xlab(expression(log~SFR~(M['\u0298']*yr^{-1})))+ylab("g-r")

quartz.save(type = 'pdf', file = '..//figures/grSFR.pdf',width = 9.5, height = 9)







ggplot(datam,aes(y=color_gr,x=lgm_tot_p50,linetype=type,color=zoo))+
  geom_point(aes(alpha=0.25,shape=bpt,size=bpt))+
  #  geom_density2d(alpha=0.7,size=1)+
  scale_color_manual(name="",values=c("red2", "cyan4"))+
  scale_size_manual(name="",values=c(2.5,4))+
  #  scale_linetype_manual(name="",values=c("solid", "dashed","solid","dashed"))+
  theme_bw()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  coord_cartesian(ylim=c(0,1.25),xlim=c(9.5,12))+xlab("Log (stellar mass)")+ylab("(g-r)")


# R proj and Mhalo

datag<-rbind(dataE,dataS)

datag2<-subset(datag,type=="E_cut" | type =="S_cut")
datag2$bpt<-as.factor(datag2$bpt)
ggplot(datag2,aes(y=lgm_tot_p50,x=RprojLW_Rvir,linetype=type,color=type))+
 geom_point(size=2,alpha=0.8,aes(shape=bpt))+
  geom_density2d(size=1.5)+
  scale_color_manual(name="",values=c("red3", "cyan4"))+
  scale_linetype_manual(name="",values=c("solid", "dashed"))+
  scale_shape_manual(values=c(1,19))+
  theme_bw()+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
ylab(expression(log~M["*"]~(M['\u0298'])))+xlab(expression(R/R[vir]))+
  coord_cartesian(ylim=c(10,12))
quartz.save(type = 'pdf', file = '..//figures/RMgal.pdf',width = 9.5, height = 9)




dat_melt0<-datag2[,c("RprojLW_Rvir","zoo","bpt")]
dat_melt1<-melt(dat_melt0,id=c("zoo","bpt"))


sample1<-subset(dat_melt1,zoo=="E")[,"value"]
sample2<-subset(dat_melt1,zoo=="S")[,"value"]
cdf1 <- ecdf(sample1)
cdf2 <- ecdf(sample2) 
# find min and max statistics to draw line between points of greatest distance
minMax <- seq(min(sample1, sample2), max(sample1, sample2), length.out=length(sample1)) 
x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
y0 <- cdf1(x0) 
y1 <- cdf2(x0) 


xE<-uniroot(function(x) cdf1(x)-0.95,lower = 0, upper = 10)$root
xS<-uniroot(function(x) cdf2(x)-0.95,lower = 0, upper = 10)$root

ggplot(dat_melt1,aes(x=value,linetype=zoo,color=zoo,fill=zoo,shape=zoo))+
  stat_ecdf(size=1.5)+
  #  geom_histogram(size=2.25,alpba=0.8,binwidth=0.25)+
  scale_color_manual(name="",values=c("red3", "royalblue3"))+
  scale_linetype_manual(name="",values=c("solid", "dashed"))+
  theme_bw()+
  theme(legend.position= "none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=20),
        strip.text.x=element_text(size=20),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=20),axis.title.x=element_text(size=rel(1)))+xlab(expression(R/R[vir]))+
  ylab("Cumulative Fraction")+
  geom_segment(aes(x = xE, y = 0, xend = xE, yend = 0.95),colour = "gray",linetype="dotted",size=0.75)+
  geom_segment(aes(x = 0, y = 0.95, xend = xE, yend = 0.95),colour = "gray",linetype="dotted",size=0.75)+
  
  geom_segment(aes(x = xS, y = 0, xend = xS, yend = 0.95),colour = "gray",linetype="dotted",size=0.75)+
  geom_segment(aes(x = 0, y =  0.95, xend = xS, yend = 0.95),colour = "gray",linetype="dotted",size=0.75)+
  coord_cartesian(xlim=c(0.05,8),ylim=c(0.025,1))+
  geom_point(x=xE,y=0.95,colour="red3",size=5)+
  geom_point(x=xS,y=0.95,colour="royalblue3",size=5)
quartz.save(type = 'pdf', file = '..//figures/CDF.pdf',width = 9.5, height = 9) 
  




ggplot(gdata,aes(y=value,x=PSM,fill=zoo,linetype=PSM))+geom_boxplot(alpha=0.5)+facet_wrap(zoo~variable,scale="free")

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(values=c("#fcbba1","#a6bddb","#de2d26","#00CED1"))
unlockBinding("ggplot",parent.env(asNamespace("GGally")))
assign("ggplot",ggplot,parent.env(asNamespace("GGally")))
require(GGally)
gp<-ggpairs(
  Full_data,columns=2:4,
  upper = list(continuous = "density", combo = "box"),
  lower = list(continuous = "points", combo = "dot"),
  color = "color ",
   title = ""
)

gp+theme_bw()

ggplot(gdata,aes(y=value,x=PSM))+geom_boxplot()+
  facet_wrap(zoo~variable,scales="free")


x<-Full_data$lgm_tot_p50
y<-Full_data$color_gr
z<-Full_data$sfr_tot_p50
with(Full_data,s3d<-{scatterplot3d(x, y, z,color=Full_data$color, col.axis="gray30",
              col.grid="gray70", pch=17)})


Full_data$color<-as.factor(Full_data$color)

scatter3d(sfr_tot_p50 ~ lgm_tot_p50 + color_gr | color, surface=FALSE, 
          ellipsoid=TRUE, revolutions=3, data=Full_data)
scatter3D(x, y, z,col=Full_data$color,colkey = FALSE,theta=-280)
