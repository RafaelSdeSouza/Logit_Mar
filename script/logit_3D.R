# JAGS Code 
#  Required libraries
library(rjags);library(ggmcmc);library(ggplot2);library(ggthemes);
library("plot3D");require(visreg)
uting

source("facet_wrap_labeller.R")
# Read and format data
dataE     <- read.table("..//data/matched_E.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S.txt",header=TRUE,na.strings="")
data<-rbind(dataE,dataS)
data_cutE <- dataE[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]
data_cutS <- dataS[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]


#s3d$plane3d(fit)


# x, y, z variables
x <- data_cutE$RprojLW_Rvir*2.061402+2.180115
Mhalo <- data_cutE$logM200_L*0.289949+13.86258
z <- data_cutE$bpt
# Compute the linear regression (z = ax + by + d)
fit <- glm(z ~ Mhalo+x,family=binomial("logit")) 



visreg2d(fit, "y", "x", plot.type = "persp",scale = "response",col = "deepskyblue2",
         xlab=expression(R/r200),ylab=expression(Mhalo),zlab="AGN fraction",theta=-25,ticktype="detailed",
         col="deepskyblue")

visreg(fit,"x",by="Mhalo",scale = "response",xlab=expression(R/r200),ylab="AGN fraction",breaks=c(13.5,14,14.5))


visreg2d(fit, "y", "x",type="contrast")


# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(1.01*min(x), 0.99*max(x), length.out = grid.lines)
y.pred <- seq(1.01*min(y), 0.99*max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)

z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit,type ="response")

# scatter plot with regression plane


xx <-range(x)
x1 <- seq(0.5*xx[1], 1.5*xx[2], length.out=100)
yy <- range(y)
y1 <- seq(0.5*yy[1], 1.5*yy[2], length.out=100)

f<-function(x,y)
  predict(fit, data.frame(x,y),type = 'response')
z1 <- outer(x1,y1,
           f)
library(rsm)
library(lattice)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#p<-persp(x,y,z, theta=150, phi=20,
#         expand = 0.5,shade = 0.1,
#         xlab="Z", ylab=expression(NII.Ha), zlab=expression(log10.EW.Ha),ticktype='detailed',
#         col = YlOrBr,border=NA,xlog=T,ylog=T)
cor = cm.colors(100)

#cairo_pdf("logit3D_2.pdf")
trellis.par.set("axis.line",list(axis.text=list(cex=20),col=NA,lty=1,lwd=2))
#par(mar=c(1,1,1,1))
wireframe(z~x+y,data=data.frame(x=x1, y=rep(y1, each=length(x1)), z=z1),
          par.settings = list(regions=list(alpha=0.6)),
          col.regions =cor,drape=T,light.source = c(5,5,5),colorkey = FALSE,
#          xlab=list(label=expression(log[Mstar]),cex=1.25), ylab=list(label=expression(V[sigma]),cex=1.25),
#          zlab=list(rot=90,label=expression(pi),cex=1.25,dist=-1,rot=0),
          scale=list(tck=0.75,arrows=FALSE,distance =c(0.75, 0.75, 0.75)))






pdf("..//Figures/Normal3D.pdf",height = 10,width = 9)
scatter3D(x = x, y = x, z=fitpoints, pch = 19, cex = 1, las=3,cex.lab=1.5,
          theta = 120, phi = 20, ticktype = "detailed",col="cyan3",bty = "b2",
          ylab="log SFR",
          xlab="log M",
          zlab="log Z", 
          surf = list(col="gray60",x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints,lwd=1.25,lty=3),colkey = FALSE)
dev.off()


surf3D(x = x.pred, y = y.pred, z = z.pred, pch = 19, cex = 1, las=3,cex.lab=1.5,
          theta = 120, phi = 20, ticktype = "detailed",col="cyan3",bty = "b2",
          ylab="log SFR",
          xlab="log M",
          zlab="log Z", 
          surf = list(col="gray60",x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints,lwd=1.25,lty=3),colkey = FALSE)
