icecream <- data.frame(
  temp=c(11.9, 14.2, 15.2, 16.4, 17.2, 18.1, 
         18.5, 19.4, 22.1, 22.6, 23.4, 25.1),
  units=c(185L, 215L, 332L, 325L, 408L, 421L, 
          406L, 412L, 522L, 445L, 544L, 614L)
)

glmModelPlot <- function(x, y, xlim,ylim, meanPred,  LwPred, UpPred, 
                         plotData, main=NULL){
  ## Based on code by Arthur Charpentier:
  ## http://freakonometrics.hypotheses.org/9593
  par(mfrow=c(1,1))
  n <- 2
  N <- length(meanPred)
  zMax <- max(unlist(sapply(plotData, "[[", "z")))*1.5
  mat <- persp(xlim, ylim, matrix(0, n, n), main=main,
               zlim=c(0, zMax), theta=-30, 
               ticktype="detailed",box=FALSE)
  C <- trans3d(x, UpPred, rep(0, N),mat)
  lines(C, lty=2)
  C <- trans3d(x, LwPred, rep(0, N), mat)
  lines(C, lty=2)
  C <- trans3d(c(x, rev(x)), c(UpPred, rev(LwPred)),
               rep(0, 2*N), mat)
  polygon(C, border=NA, col=adjustcolor("gray90", alpha.f = 0.5))
  C <- trans3d(x, meanPred, rep(0, N), mat)
  lines(C, lwd=2, col="grey")
  C <- trans3d(x, y, rep(0,N), mat)
  points(C, lwd=2, col="#00526D")
  for(j in N:1){
    xp <- plotData[[j]]$x
    yp <- plotData[[j]]$y
    z0 <- plotData[[j]]$z0
    zp <- plotData[[j]]$z
    C <- trans3d(c(xp, xp), c(yp, rev(yp)), c(zp, z0), mat)
    polygon(C, border=NA, col="gray95", density=40)
    C <- trans3d(xp, yp, z0, mat)
    lines(C, lty=2)
    C <- trans3d(xp, yp, zp, mat)
    lines(C, col=adjustcolor("black", alpha.f = 0.5))
  }
}

market.size <- 800
icecream$opportunity <- market.size - icecream$units

xlim <- c(min(icecream$temp)*0.95, max(icecream$temp)*1.05)
ylim <- c(floor(min(icecream$units)*0.95),
          ceiling(max(icecream$units)*1.05))


bin.glm <- glm(cbind(units, opportunity) ~ temp, data=icecream, 
               family=binomial(link = "logit"))
par(mfrow=c(2,2))
plot(bin.glm)
title(outer=TRUE, line = -1,
      main = list("Binomial (logit) GLM", 
                  cex=1.25,col="black", font=2)) 

meanProb <- predict(bin.glm, type="response")
meanPred <- meanProb*market.size
UpPred <- qbinom(.95, market.size, meanProb)
LwPred <- qbinom(.05, market.size, meanProb)

plotData <- lapply(
  seq(along=icecream$temp),
  function(i){
    y = ylim[1]:ylim[2]
    x = rep(icecream$temp[i], length(y))
    z0 = rep(0, length(y))
    z = dbinom(y, market.size, meanProb[i])
    return(list(x=x, y=y, z0=z0, z=z))
  }
)

glmModelPlot(x = icecream$temp, y=icecream$units,
             xlim=xlim, ylim=ylim,
             meanPred = meanPred, LwPred = LwPred,
             UpPred = UpPred, plotData = plotData,
             main = "Binomial (logit) GLM")