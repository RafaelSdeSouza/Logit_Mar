#scatter 3D
require(scatterplot3d)
require(plyr)
require(plot3D)


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
dataEoriginal$color<-rep("red",nrow(dataEoriginal))
dataSoriginal$color<-rep("cyan",nrow(dataSoriginal))

dataEoriginal$PSM<-rep("Before",nrow(dataEoriginal))
dataSoriginal$PSM<-rep("Before",nrow(dataSoriginal))

dataoriginal<-rbind(dataEoriginal,dataSoriginal)


dataraw<-dataoriginal[,c("bpt" ,"lgm_tot_p50","sfr_tot_p50","color_gr","zoo","color","PSM")]
#-------------------------



dataE     <- read.table("..//data/matched_E_original.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S_original.txt",header=TRUE,na.strings="")
dataE$color<-rep("red4",nrow(dataE))
dataS$color<-rep("cyan4",nrow(dataS))

dataS$PSM<-rep("After",nrow(dataS))
dataE$PSM<-rep("After",nrow(dataE))

datam<-rbind(dataE,dataS)
datam<-datam[,c("bpt" ,"lgm_tot_p50","sfr_tot_p50","color_gr","zoo","color","PSM")]


Full_data<-rbind(dataraw,datam)
gdata<-melt(Full_data,id=c("bpt","zoo","color","PSM"))


ggplot <- function(...) ggplot2::ggplot(...) + scale_color_manual(values=c("#00CED1","#00CED1","#00CED1","#00CED1"))
unlockBinding("ggplot",parent.env(asNamespace("GGally")))
assign("ggplot",ggplot,parent.env(asNamespace("GGally")))
require(GGally)
gp<-ggpairs(
  Full_data,columns=2:4,
  upper = list(continuous = "density", combo = "box"),
  lower = list(continuous = "points", combo = "dot"),
  color = "color ",
  alpha=0.7,
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
