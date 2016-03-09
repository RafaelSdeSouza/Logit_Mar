#Plot posteriors

gplot<-read.table("..//data/gplot.dat",header=TRUE)

gplot$gal<-as.factor(gplot$gal)

gplot$Parameter<-revalue(gplot$Parameter, c("beta[1,1]"= "beta[1]","beta[1,2]"="beta[1]", "beta[2,1]" ="beta[2]", "beta[2,2]"="beta[2]", "beta[3,1]"="beta[3]",
                                            "beta[3,2]"="beta[3]"))


ggplot(data=gplot,aes(x=value,group=gal,fill=gal))+
  geom_density(colour="white",size=0.01,alpha=0.8)+facet_grid(Parameter~gal,labeller = label_parsed)+
  theme_bw()+
  theme(legend.position="none",panel.background = element_rect(fill = "white"),plot.background = element_rect(
    fill = "white"),plot.title = element_text(hjust=0.5),
    axis.title.y=element_text(vjust=0.75),axis.text.x=element_text(size=22),
    axis.text.y=element_text(size=22),
    strip.text.x=element_text(size=25),
    axis.title.x=element_text(vjust=-0.25),
    text = element_text(size=22),axis.title.x=element_text(size=rel(1)),strip.background=element_rect(
      fill = "white"),strip.text=element_text(
        size = 25))+
  geom_vline(xintercept=0,size=1,linetype="dashed",colour=c("grey50")) +
  #  scale_fill_manual(values=c("#E0FFFF","#00CED1","cyan4"))+
  scale_fill_manual(values=c("#de2d26","#00CED1"))+
  ylab("Posterior")+xlab("Parameter value")


quartz.save(type = 'pdf', file = '..//figures/betas.pdf',width = 8.5, height = 9)