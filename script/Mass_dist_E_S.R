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

scatter3D(x2, x1, y, pch = 19, cex = 0.5, cex.lab=1.5,colvar=color,
          theta = 130, phi = 20, ticktype = "detailed",bty = "b2",t="l",
          xlab="log SFR",
          ylab="log Mstar",
          zlab="R/r200")




# Cut x and y variables in bins for counting
x.binE <- seq(10, max(x1E), length.out = 10)
y.binE <- seq(min(x2E), max(x2E), length.out = 10)
xyE <- table(cut(x, x.binE), cut(y, y.binE))
zE <- xyE



hist3D (x.binE[-1],y.binE[-1],z=zE,
        bty = "b2", phi = 30,  theta = -145,
        xlab = "Mstar", ylab = "R/r200", zlab = "Number of galaxies",
        col = c("#de2d26"), border = "gray", shade = 0.2,
        ticktype = "detailed", space = 0.15, d = 2, add = F,  colkey = F,parse=T,along="y",
        ylim=c(0.2,7.5),xlim=c(10,12))


x1S<- dataS$lgm_tot_p50
x2S <- dataS$RprojLW_Rvir

# Cut x and y variables in bins for counting
x.binS <- seq(10, max(x1S), length.out = 10)
y.binS <- seq(min(x2S), max(x2S), length.out = 10)
xyS <- table(cut(x, x.binS), cut(y, y.binS))
zS <- xyS



hist3D (x.binS[-1],y.binS[-1],z=zS,
        bty = "b2", phi = 30,  theta = -145,
        xlab = "Mstar", ylab = "R/r200", zlab = "Number of galaxies",
        col = "#E0FFFF", border = "gray", shade = 0.2,
        ticktype = "detailed", space = 0.15, d = 2, add = F,  colkey = F,parse=T,along="y",
        ylim=c(0.2,7.5),xlim=c(10,11.75))





scatterplot3d(x, y, z, highlight.3d=TRUE,
              col.axis="blue", col.grid="lightblue",
              main="scatterplot3d - 2", pch=20)


+coord_cartesian(xlim=c(10,12))

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

