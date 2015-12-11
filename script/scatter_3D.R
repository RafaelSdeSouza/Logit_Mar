#scatter 3D
require(scatterplot3d)

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

unmatched<-data_n2[,c("bpt" ,"lgm_tot_p50","sfr_tot_p50","color_gr","zoo")]
#-------------------------
unmatched


dataE     <- read.table("..//data/matched_E_original.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S_original.txt",header=TRUE,na.strings="")
dataE$color<-rep("red2",nrow(dataE))
dataS$color<-rep("cyan3",nrow(dataS))
datam<-rbind(dataE,dataS)

x<-datam$color_gr
y<-datam$lgm_tot_p50
z<-datam$sfr_tot_p50
with(datam,s3d<-{scatterplot3d(x, y, z,color=datam$color, col.axis="blue",
              col.grid="lightblue", pch=20)})






