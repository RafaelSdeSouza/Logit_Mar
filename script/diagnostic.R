require(rstan)
require(lme4)
require(glmer2stan)
require(sjmisc)
require(binomTools)
require(sjPlot)
data     <- read.table("..//data/matched.txt",header=TRUE,na.strings="")
data_cut <- data[,c("bpt","logM200_L","RprojLW_Rvir","zoo")]
data_cut   <- subset(data_cut, zoo=="E")


fit <- glm(bpt ~ logM200_L+RprojLW_Rvir+I(logM200_L^2)+I(RprojLW_Rvir^2), family=binomial("logit"),data=data_cut)
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