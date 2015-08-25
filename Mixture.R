library(mixtools)
library(mclust)
data<-read.table("data_index_AMADA.txt",header=TRUE)


g.r<-data$g.r
FUV<-data$FUV.NUV


mixmdl = normalmixEM(g.r)
pdf("gr.pdf",height = 8,width = 8)
plot(mixmdl,which=2,xlab2="g-r",ylab2="Densidade")
lines(density(g.r), lty=2, lwd=2)
dev.off()

mixmdl3 = normalmixEM(FUV)
pdf("FUV.NUV.pdf",height = 8,width = 8)
plot(mixmdl3,which=2,xlab2="FUV-NUV",ylab2="Densidade")
lines(density(FUV), lty=2, lwd=2)
dev.off()

data2<-data.frame(data$g.r,data$FUV.NUV)
mixmdl2 = mvnormalmixEM(data2)
pdf("bivariate.pdf",height = 8,width = 8)
plot(mixmdl2,which=2,ylab2="FUV-NUV",xlab2="g-r")
dev.off()


colnames(data2)<-c("g-r","FUV-NUV")
M1<-Mclust(data2)
pdf("Mixture_models.pdf",height = 8,width = 8)
plot(M1, what = "classification")
dev.off()

pdf("BIC.pdf",height = 8,width = 8)
plot(M1, what = "BIC")
dev.off()

plot(M1, what = "density", type = "persp")

M2<-densityMclust(FUV)
plot(M2, what = "density", data = FUV, breaks = 50)

M2$classification

gdata<-data.frame(FUV=M2$data,class=M2$classification)
gdata$class<-as.factor(gdata$class)
ggplot(gdata,aes(x=FUV,group=class))+
  geom_histogram(binwidth=0.1,aes(fill=class,alpha=0.2))

ggplot(gdata,aes(x=FUV))+
  geom_density()
