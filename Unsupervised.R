require(plyr)
library(caret)
require(kernlab)
library(e1071)
require(MASS)
require(mclust)
AGN_data<-read.table("../data/outputdata_diagnostic.txt",header=TRUE,sep="")

# Format data for WHAN test
WHAN<-AGN_data[,c("log10..NII..Ha.","log10.EW.Ha..","WHAN_Class")]
WHAN$WHAN_Class<-as.factor(WHAN$WHAN_Class)
#write.matrix(WHAN,"../data/WHAN.txt")




data<-WHAN[,-3]
mod4 = Mclust(data)
summary(mod4)
plot(mod4, what = "classification")
plot(mod4, what = "boundaries", ngrid = 200)



ggplot(WHAN,aes(x=log10..NII..Ha.,y=log10.EW.Ha..,colour=WHAN_Class))+
  geom_point()+theme_stata()


library(fpc)
# eps is radius of neighborhood, MinPts is no of neighbors
# within eps
 cluster <- dbscan(data, eps=0.6, MinPts=4)
 plot(cluster, data)

# Notice points in cluster 0 are unassigned outliers
 table(cluster$cluster, WHAN$WHAN_Class)