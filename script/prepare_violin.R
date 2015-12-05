#prepare data for violin



# Read and format data
data     <- read.csv("..//data/sample_CRP02_sub.csv",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","lgm_tot_p50","logM200_L","RprojLW_Rvir","sfr_tot_p50","color_gr","zoo")]
data2    <- na.omit(data_cut)
data2    <- data2[data2$lgm_tot_p50>0,]
data2    <- data2[which(data2$logM200_L>0),]
data2    <- data2[which(data2$RprojLW_Rvir>=0),]
data2    <- data2[which(data2$sfr_tot_p50>=-100),]

# Standardized variables

data2<-data.frame(bpt=data2$bpt,as.data.frame(scale(data2[,2:6])),zoo=data2$zoo)
#trainIndex <- sample(1:nrow(data2),nrow(data2))
data3      <- data2
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
data_n2    <- subset(data_n, zoo=="E" | zoo == "S")

#before match
unmatched<-data_n2[,c("bpt" ,"lgm_tot_p50","sfr_tot_p50","color_gr","zoo")]
unmatched$type<-rep("Before PSM", nrow(unmatched))

dataE     <- read.table("..//data/matched_E.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S.txt",header=TRUE,na.strings="")
datam<-rbind(dataE,dataS)

matched<-datam[,c("bpt","lgm_tot_p50","sfr_tot_p50","color_gr","zoo")]

matched$type<-rep("After PSM", nrow(matched))

full<-rbind(unmatched,matched,deparse.level=2)
colnames(full)<-c("AGN","lgm_tot_p50","sfr_tot_p50","color_gr","Galaxy","type")
full$AGN<-revalue(full$AGN,c("1"="Yes","0"="No"))

full_E<-subset(full, Galaxy=="E")
full_S<-subset(full, Galaxy=="S")



write.csv(full_E,"..//data/full_E.csv")

write.csv(full_S,"..//data/full_S.csv")

#gmelt0<-melt(unmatched, id.vars = c('bpt'))

#CairoPDF("..//figures/before.pdf",height=12,width = 10)
#ggplot(aes(x=value,group=bpt,fill=bpt),data=gmelt0)+geom_density()+
#  facet_grid(variable~bpt)+scale_fill_stata()+theme_hc()+coord_cartesian(xlim=c(-2,4))
#dev.off()
#after match
#gmelt <- melt(matched[,1:4], id.vars = c('bpt'))







