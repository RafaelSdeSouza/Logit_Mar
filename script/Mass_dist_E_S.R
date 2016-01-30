# mass dist for early- late-type galaxies

require(plyr)
require(reshape2)
require(ggplot2)
require(ggthemes)


dataE     <- read.table("..//data/matched_E_original.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S_original.txt",header=TRUE,na.strings="")
datam<-rbind(dataE,dataS)


ggplot(datam,aes(x=lgm_tot_p50,color=zoo,linetype=zoo))+
  geom_histogram(size=1.75,binwidth=0.1,fill="white")+ theme_bw()+
  scale_color_manual(name="",values=c("red2", "royalblue2"))+
  theme(legend.position="none",plot.title = element_text(hjust=0.5),
        axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=25),
        strip.text.x=element_text(size=25),
        axis.title.x=element_text(vjust=-0.25),
        text = element_text(size=25),axis.title.x=element_text(size=rel(1)))+
  xlab(expression(log~M["*"]))+ylab("Number of galaxies")

# scatter plot with regression plane
library("plot3D")




x1E<- dataE$lgm_tot_p50
x2E <- dataE$RprojLW_Rvir

color<-as.numeric(datam$zoo)






# Cut x and y variables in bins for counting
x.binE <- seq(10, 11.9, length.out = 10)
y.binE <- seq(min(x2E), max(x2E), length.out = 10)
xyE <- table(cut(x1E, x.binE), cut(x2E, y.binE))
zE <- xyE



hist3D (x.binE[-1],y.binE[-1],z=zE,
        bty = "b2", phi = 30,  theta = -145,
        xlab = "log Mstar", ylab = "R/r200", zlab = "Number of galaxies",
        col = ramp.col(col = c("plum2", "red3"),n = 50, alpha = 0.7), border = "gray", shade = 0.2,
        ticktype = "detailed", space = 0.15, d = 2, add = F,  colkey = list(length=0.45,dist=-0.075,shift=0.1),parse=T,along="y",
        ylim=c(0.2,7.75),xlim=c(10.1,12))
quartz.save(type = 'pdf', file = '..//figures/hist3D_E.pdf',width = 12, height = 10)

x1S<- dataS$lgm_tot_p50
x2S <- dataS$RprojLW_Rvir

# Cut x and y variables in bins for counting
x.binS <- seq(10, 11.9, length.out = 10)
y.binS <- seq(min(x2S), max(x2S), length.out = 10)
xyS <- table(cut(x1S, x.binS), cut(x2S, y.binS))
zS <- xyS



hist3D (x.binS[-1],y.binS[-1],z=zS,
        bty = "b2", phi = 30,  theta = -145,
        xlab = "log Mstar", ylab = "R/r200", zlab = "Number of galaxies",
         border = "gray", shade = 0.2,col=ramp.col(col = c("cyan2", "blue4"), n = 50, alpha = 0.7),
        ticktype = "detailed", space = 0.15, d = 2, add = F,  colkey = list(length=0.45,dist=-0.075,shift=0.1),
        parse=T,along="y",
        ylim=c(0.2,7.75),xlim=c(10.1,12))

quartz.save(type = 'pdf', file = '..//figures/hist3D_S.pdf',width = 12, height = 10)




