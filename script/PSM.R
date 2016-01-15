##  R script for propensity score matching 
#  Copyright (C) 2015  COIN
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as published by
#the Free Software Foundation.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

# This function performs  propensity score matching using the package MatchIt

# Required libraries
library(MatchIt);library(plyr);library(reshape2)



# Read and format data
data     <- read.csv("..//data/sample_CRP02_sub.csv",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","lgm_tot_p50","logM200_L","RprojLW_Rvir","sfr_tot_p50","color_gr","zoo","groupID","igal")]
data2    <- na.omit(data_cut)
data2    <- data2[data2$lgm_tot_p50>0,]
data2    <- data2[which(data2$logM200_L>0),]
data2    <- data2[which(data2$RprojLW_Rvir>=0),]
data2    <- data2[which(data2$sfr_tot_p50>=-100),]

# Standardized variables

#data2<-data.frame(bpt=data2$bpt,as.data.frame(scale(data2[,c("lgm_tot_p50","sfr_tot_p50","color_gr")])),zoo=data2$zoo,
#                  RprojLW_Rvir=data2$RprojLW_Rvir,logM200_L=data2$logM200_L,groupID=data2$groupID)

data2<-data.frame(bpt=data2$bpt,data2[,2:6],zoo=data2$zoo,groupID=data2$groupID,igal=data2$igal)

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
#data_n2    <- subset(data_n, zoo=="E" | zoo == "S")
#data_n2$zoo<-droplevels(data_n2$zoo)


data_n2    <- subset(data_n, zoo == "S")
data_n2$zoo<-droplevels(data_n2$zoo)

#m.out <- matchit(formula = bpt ~ lgm_tot_p50 + sfr_tot_p50 + color_gr, data=data_n2[,c("bpt" ,"lgm_tot_p50","sfr_tot_p50","color_gr")], method = 	"nearest", distance = "logit")

m.out <- matchit(formula = bpt ~ lgm_tot_p50 + sfr_tot_p50 + color_gr, data=data_n2, method = 	"nearest", distance = "logit")

matched <- match.data(m.out)

write.matrix(matched,"..//data/matched_S_original.txt",sep=" ")





# Script ends here




ggplot(data=matched,aes(x=sfr_tot_p50,group=bpt))+geom_density()










# plot distribution of Mass before and after match 


#before match
#unmatched<-data_n2[,c("bpt" ,"lgm_tot_p50","sfr_tot_p50","color_gr")]
#write.matrix(unmatched,"..//data/unmatched.txt",sep=" ")

#gmelt0<-melt(unmatched, id.vars = c('bpt'))

#CairoPDF("..//figures/before.pdf",height=12,width = 10)
#ggplot(aes(x=value,group=bpt,fill=bpt),data=gmelt0)+geom_density()+
#  facet_grid(variable~bpt)+scale_fill_stata()+theme_hc()+coord_cartesian(xlim=c(-2,4))
#dev.off()
#after match
#gmelt <- melt(matched[,1:4], id.vars = c('bpt'))

#cairo_pdf("..//figures/after.pdf",height=12,width = 10)
#ggplot(aes(x=value,group=bpt,fill=bpt),data=gmelt)+geom_density()+
#  facet_grid(variable~bpt)+scale_fill_stata()+theme_hc()+coord_cartesian(xlim=c(-2,4))
#dev.off()




#matched$mytype<-rep("matched",nrow(matched))
#unmatched$mytype<-rep("unmatched",nrow(unmatched))
#full<-rbind(matched[,c(1,2,3,4,7)],unmatched)


#full$bpt  <- revalue(full$bpt ,c("0"="AGN 0","1"="AGN 1"))
#full$lgm_tot_p50<-round(full$lgm_tot_p50,2)
#full$sfr_tot_p50<-round(full$sfr_tot_p50,2)
#full$color_gr<-round(full$color_gr,2)

#write.csv(full,"..//data/full.csv",row.names = FALSE)

#require(beanplot)

#gmelt$group<-paste(gmelt$variable,gmelt$bpt,sep=" ")

#beanplot(value~group,border = NA,data=gmelt,side = "both",ll = 0.0001,
#         col = list("blue", c("red", "white")),ylim=c(-4,4))



#gmelt0$group<-paste(gmelt0$variable,gmelt0$bpt,sep=" ")
#gmelt0$value<-as.numeric(gmelt0$value)

#beanplot(value~group,border = NA,data=gmelt0,side = "both",ll = 0.0001,
#         col = list("blue", c("red", "white")),ylim=c(-4,4))
