require(rstan)
require(lme4)
require(glmer2stan)
require(sjmisc)
require(binomTools)
require(sjPlot)
data     <- read.table("..//data/matched.txt",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]
data_cut   <- subset(data_cut, zoo == "S")


mod<- glm(bpt ~ RprojLW_Rvir+logM200_L, family=binomial("logit"),data=data_cut)
plot(fitted(mod), jitter(data_cut$bpt, amount=0.05), xlab="Fitted values",
                      ylab="Probability of presence", las=1, cex.lab=1.2, cex=0.8)

abline(0,1, lty=3)
t.breaks <-cut(fitted(mod), seq(0,1, by=0.05))
means <-tapply(data_cut$bpt, t.breaks, mean)
semean <-function(x) sd(x)/sqrt(length(x))
means.se <-tapply(data_cut$bpt, t.breaks, semean)
points(seq(0.05, 0.95, by=0.1), means, pch=16, col="orange")
segments(seq(0.05, 0.95, by=0.1), means-2*means.se,
         seq(0.05, 0.95,by=0.1), means+2*means.se,lwd=2, col="orange")




## summary(budworm.lg)

(Rsq.budworm <- Rsq(fit))

plot(Rsq.budworm, "ROC")



fit <- glm(bpt ~ logM200_L+RprojLW_Rvir ,
             data_cut,
             family = binomial("logit"))
sjp.glm(fit)

sjp.glm(fit, type = "qq")

m1_lme4 <- lmer( Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE )

# construct subject index --- glmer2stan forces you to manage your own cluster indices
sleepstudy$subject_index <- as.integer(as.factor(sleepstudy$Subject))

# fit with glmer2stan
m1_g2s <- lmer2stan( Reaction ~ Days + (Days | subject_index), data=sleepstudy )
m1_g2s@stanmodel