dataE     <- read.table("..//data/matched_E_original.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S_original.txt",header=TRUE,na.strings="")
datam<-rbind(dataE,dataS)



AGN_R<-datam[,c("bpt","RprojLW_Rvir")]
cut(dataE$RprojLW_Rvir,breaks=c(0,1,2,3,4,5,6,7,8))


t.breaks <-cut(dataE$RprojLW_Rvir, breaks=c(0,1.5,3,4.5,6,8))
means <-tapply(dataE$bpt, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
means.se <-tapply(dataE$bpt, t.breaks, semean)
bins<-levels(t.breaks)



gdata<-data.frame(x=bins,y=means)
gdata$gal<-rep("E",nrow(gdata))


t.breaks2 <-cut(dataS$RprojLW_Rvir, breaks=c(0,1.5,3,4.5,6,8))
means2 <-tapply(dataS$bpt, t.breaks2, mean)
semean2 <-function(x) sd(x)/sqrt(length(x))
means.se2 <-tapply(dataS$bpt, t.breaks2, semean2)

gdata2<-data.frame(x=bins,y=means2)
gdata2$gal<-rep("S",nrow(gdata2))


ggplot(aes(x=x,y=y),data=gdata)+
  geom_point(aes(x=x,y=y),size=3,data=gdata,colour="#de2d26")+
  geom_point(aes(x=x,y=y),size=3,data=gdata2,colour="cyan")+
  geom_errorbar(data=gdata,guide="none",aes(x=x,y=y,ymin=y-2*means.se,ymax=y+2*means.se),alpha=0.7,
                colour="#de2d26",width=0.05)+
  geom_errorbar(data=gdata2,guide="none",aes(x=x,y=y,ymin=y-2*means.se2,ymax=y+2*means.se2),alpha=0.7,
                colour="cyan3",width=0.05)+
  theme_bw()+theme(legend.position="top",plot.title = element_text(hjust=0.5),
                   axis.title.y=element_text(vjust=0.75),
                   axis.title.x=element_text(vjust=-0.25),
                   text = element_text(size=25))+xlab(expression(R/R[vir]))+
  ylab("Observed AGN Fraction")+facet_wrap(~gal)
 

