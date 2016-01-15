##################################################################################################
# PURPOSE:
#    Plot BPT diagram of ellipticals and spiral galaxies
# CALLING SEQUENCE:
#    
# INPUTS:
#    ../data/sample_CRP02_sub.csv
# PARAMETERS:
#    
# OUTPUT:
#    ../figures/BPT.pdf
# REQUIRED SCRIPTS:
#   
##################################################################################################

read_data = T

if(read_data == T){
  print('Reading data......')
  remove(list = ls())
  data = read.csv('../data/sample_CRP02_sub.csv')
  data = data.frame(data, logO3O2 = log10(data$oiii_5007_flux / (data$oii_3726_flux + data$oii_3729_flux)))
  print('Done!')
}

#-----------------------
# OUTPUT FILE
#-----------------------
output_file = '../figures/BPT2.pdf'

#-----------------------
# Pos PSM
#-----------------------

dataE     <- read.table("..//data/matched_E_original.txt",header=TRUE,na.strings="")
dataS     <- read.table("..//data/matched_S_original.txt",header=TRUE,na.strings="")


#-----------------------
# SUBSAMPLES
#-----------------------
xx.xxS = data$zoo == 'S' & is.na(data$zoo) == F
xx.xxE = data$zoo == 'E' & is.na(data$zoo) == F
xx.xxS_sey = xx.xxS & data$bpt == "Seyfert" | data$bpt == "Seyfert/LINER" & is.na(data$bpt) == F
xx.xxE_sey = xx.xxE & data$bpt == "Seyfert" | data$bpt == "Seyfert/LINER" & is.na(data$bpt) == F
xx.xxS_SF = xx.xxS & data$bpt == "Star Forming"  & is.na(data$bpt) == F & !is.na(match(data$igal,dataS$igal))
xx.xxE_SF = xx.xxE & data$bpt == "Star Forming"  & is.na(data$bpt) == F & !is.na(match(data$igal,dataE$igal))
xx.xxS_bpt = xx.xxS & is.na(data$logO3Hb) == F & is.na(data$logN2Ha) == F 
xx.xxE_bpt = xx.xxE & is.na(data$logO3Hb) == F & is.na(data$logN2Ha) == F 
#xx.xxS_oii = xx.xxS & is.na(data$logO3O2) == F & is.na(data$logN2Ha) == F 
#xx.xxE_oii = xx.xxE & is.na(data$logO3O2) == F & is.na(data$logN2Ha) == F 





#-----------------------
# PLOT DEFINITIONS
#-----------------------
cc = 0.0
pch1 = 18; cex1 = 0.3; col1 = rgb(cc, cc, cc, 0.8)
pch2 = 18; cex2 = 0.3; col2 = rgb(cc, cc, cc, 0.8)
pch1sey = 18; cex1sey = 1; col1sey = rgb(0.9, 0, 0, 0.8)
pch2sey = 18; cex2sey = 1; col2sey = rgb(0.9, 0, 0, 0.8)

#-----------------------
# OPEN FIGURE FILE
#-----------------------
pdf(output_file, width = 13, height = 8)
par(mfrow = c(1, 2), mar = c(5, 5, 2, 2), cex.lab = 1.5, cex.axis = 1.5)

#-----------------------
# BPT SPIRAL GALAXIES
#-----------------------
plot(data$logN2Ha[xx.xxS], data$logO3Hb[xx.xxS], xlim = c(-1, 0.5), ylim = c(-1.5, 1.5),
     xlab = expression(paste('log ([NII]/H', alpha, ')')), 
     ylab = expression(paste('log ([OIII]/H', beta, ')')),
     xaxt = 'n', yaxt = 'n', pch = pch1, cex = cex1, col = col1, lwd = 1.5)
points(data$logN2Ha[xx.xxS_sey], data$logO3Hb[xx.xxS_sey], pch = pch1sey, cex = cex1sey, col = col1sey)
points(data$logN2Ha[xx.xxS_SF], data$logO3Hb[xx.xxS_SF], pch = 20, cex = 1, col = "cyan3")

axis(1, at = seq(-2, 2, 0.5), labels = T, tcl = -0.6)
axis(1, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)
axis(3, at = seq(-2, 2, 0.5), labels = F, tcl = -0.6)
axis(3, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)

axis(2, at = seq(-2, 2, 0.5), labels = T, tcl = -0.6)
axis(2, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)
axis(4, at = seq(-2, 2, 0.5), labels = F, tcl = -0.6)
axis(4, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)

leg1 = expression(paste('SPIRALS'))
par(font = 2)
legend(-1, 1.5, c(leg1), bty = 'n', pt.cex = c(cex1), pch = c(pch1), 
       col = c(col1), cex = 1.2)
par(font = 1)

xx = seq(-4, 0.0, 0.01)
#Ka = 0.61 / (xx - 0.05) + 1.30
Ka = 0.61 / (xx - 0.05) + 1.30
#Ka = 0.33 / (xx +0.11) + 1.01
lines(xx, Ka, col = 'red', lwd = 1.5)

xx = seq(-4, 0.4, 0.01)
Ke = 0.61 / (xx - 0.47) + 1.19
lines(xx, Ke, col = 'darkgreen', lwd = 1.5)



#xx = seq(-4, 0.4, 0.01)
#Gra=(-30.787 +1.1358*xx + 0.27297*xx^2)*tanh(5.7409*xx)- 31.093
#lines(xx, Gra, col = 'violet', lwd = 1.5)



xx = seq(-0.43, 5, 0.01)
Sey = 1.01 * xx + 0.48
lines(xx, Sey, col = 'black', lwd = 2, lty = 2)

legend(-0.7, -1.2, 'Kauffmann+03', bty = 'n', text.col = 'red')
legend(0.15, -1.2, 'Kewley+01', bty = 'n', text.col = 'darkgreen')
legend(0.15, 0.7, 'Kewley+06', bty = 'n', text.col = 'black')

par(font = 2)
legend(-0.28, -1.2, 'Composite', bty = 'n', cex = 1.2)
par(font = 1)

#-----------------------
# BPT ELLIPTICAL GALAXIES
#-----------------------
plot(data$logN2Ha[xx.xxE], data$logO3Hb[xx.xxE], xlim = c(-1, 0.5), ylim = c(-1.5, 1.5),
     xlab = expression(paste('log ([NII]/H', alpha, ')')), 
     ylab = expression(paste('log ([OIII]/H', beta, ')')),
     xaxt = 'n', yaxt = 'n', pch = pch2, cex = cex2, col = col2, lwd = 1)
points(data$logN2Ha[xx.xxE_sey], data$logO3Hb[xx.xxE_sey], pch = pch2sey, cex = cex2sey, col = col2sey)
points(data$logN2Ha[xx.xxE_SF], data$logO3Hb[xx.xxE_SF], pch = 20, cex = 1, col = "cyan3")

axis(1, at = seq(-2, 2, 0.5), labels = T, tcl = -0.6)
axis(1, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)
axis(3, at = seq(-2, 2, 0.5), labels = F, tcl = -0.6)
axis(3, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)

axis(2, at = seq(-2, 2, 0.5), labels = T, tcl = -0.6)
axis(2, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)
axis(4, at = seq(-2, 2, 0.5), labels = F, tcl = -0.6)
axis(4, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)

leg2 = expression(paste('ELLIPTICALS'))
par(font = 2)
legend(-1, 1.5, c(leg2), bty = 'n', pt.cex = c(cex2), pch = c(pch2), 
       col = c(col2), cex = 1.2)
par(font = 1)

xx = seq(-4, 0.0, 0.01)
Ka = 0.61 / (xx - 0.05) + 1.30
lines(xx, Ka, col = 'red', lwd = 1.5)

xx = seq(-4, 0.4, 0.01)
Ke = 0.61 / (xx - 0.47) + 1.19
lines(xx, Ke, col = 'darkgreen', lwd = 1.5)

xx = seq(-0.43, 5, 0.01)
Sey = 1.01 * xx + 0.48
lines(xx, Sey, col = 'black', lwd = 2, lty = 2)

legend(-0.7, -1.2, 'Kauffmann+03', bty = 'n', text.col = 'red')
legend(0.15, -1.2, 'Kewley+01', bty = 'n', text.col = 'darkgreen')
legend(0.15, 0.7, 'Kewley+06', bty = 'n', text.col = 'black')

par(font = 2)
legend(-0.28, -1.2, 'Composite', bty = 'n', cex = 1.2)
par(font = 1)





dev.off()


# -----------------------------------------
composite = data$logO3Hb[xx.xxE] >  0.61 / (data$logN2Ha[xx.xxE] - 0.05) + 1.3

xx = data$logN2Ha[xx.xxE | xx.xxS]
yy = data$logO3Hb[xx.xxE | xx.xxS]
#print(length(xx[is.na(xx) == F & is.na(yy) == F]))
#print(which(is.na(xx) == T & is.na(yy) == T))
