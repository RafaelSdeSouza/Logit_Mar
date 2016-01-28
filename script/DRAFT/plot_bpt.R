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

read_data = F

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
output_file = '../figures/BPT.pdf'


#-----------------------
# SUBSAMPLES
#-----------------------
xx.xxS = data$zoo == 'S' & is.na(data$zoo) == F
xx.xxE = data$zoo == 'E' & is.na(data$zoo) == F
xx.xxS_sey = xx.xxS & data$bpt == "Seyfert" & is.na(data$bpt) == F
xx.xxE_sey = xx.xxE & data$bpt == "Seyfert" & is.na(data$bpt) == F
xx.xxS_bpt = xx.xxS & is.na(data$logO3Hb) == F & is.na(data$logN2Ha) == F 
xx.xxE_bpt = xx.xxE & is.na(data$logO3Hb) == F & is.na(data$logN2Ha) == F 
xx.xxS_oii = xx.xxS & is.na(data$logO3O2) == F & is.na(data$logN2Ha) == F 
xx.xxE_oii = xx.xxE & is.na(data$logO3O2) == F & is.na(data$logN2Ha) == F 

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
pdf(output_file, width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(5, 5, 2, 2), cex.lab = 1.5, cex.axis = 1.5)

#-----------------------
# BPT SPIRAL GALAXIES
#-----------------------
plot(data$logN2Ha[xx.xxS], data$logO3Hb[xx.xxS], xlim = c(-1, 0.5), ylim = c(-1.5, 1.5),
     xlab = expression(paste('log ([NII]/H', alpha, ')')), 
     ylab = expression(paste('log ([OIII]/H', beta, ')')),
     xaxt = 'n', yaxt = 'n', pch = pch1, cex = cex1, col = col1, lwd = 1.5)
points(data$logN2Ha[xx.xxS_sey], data$logO3Hb[xx.xxS_sey], pch = pch1sey, cex = cex1sey, col = col1sey)

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
Ka = 0.61 / (xx - 0.05) + 1.30
lines(xx, Ka, col = 'red', lwd = 1.5)

xx = seq(-4, 0.4, 0.01)
Ke = 0.61 / (xx - 0.47) + 1.19
lines(xx, Ke, col = 'darkgreen', lwd = 1.5)

xx = seq(-0.43, 5, 0.01)
Sey = 1.01 * xx + 0.48
lines(xx, Sey, col = 'black', lwd = 2, lty = 2)

#-----------------------
# BPT ELLIPTICAL GALAXIES
#-----------------------
plot(data$logN2Ha[xx.xxE], data$logO3Hb[xx.xxE], xlim = c(-1, 0.5), ylim = c(-1.5, 1.5),
     xlab = expression(paste('log ([NII]/H', alpha, ')')), 
     ylab = expression(paste('log ([OIII]/H', beta, ')')),
     xaxt = 'n', yaxt = 'n', pch = pch2, cex = cex2, col = col2, lwd = 1)
points(data$logN2Ha[xx.xxE_sey], data$logO3Hb[xx.xxE_sey], pch = pch2sey, cex = cex2sey, col = col2sey)

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

#############################################
# ALTERNATIVE DIGNOSTIC DIAGRAM
#   LOG([NII]/Ha) VS. LOG([OIII]/[OII])
#   Cid Fernandes+2010
#############################################
#-----------------------
# SPIRAL GALAXIES
#-----------------------
plot(data$logN2Ha[xx.xxS], data$logO3O2[xx.xxS], xlim = c(-1, 0.5), ylim = c(-1.5, 1.5),
     xlab = expression(paste('log ([NII]/H', alpha, ')')), 
     ylab = expression(paste('log ([OIII]/[OII])')),
     xaxt = 'n', yaxt = 'n', pch = pch1, cex = cex1, col = col1, lwd = 1.5)

points(data$logN2Ha[xx.xxS_sey], data$logO3O2[xx.xxS_sey], pch = pch1sey, cex = cex1sey, col = col1sey)

xx = data$logN2Ha[xx.xxS_sey]
yy = data$logO3O2[xx.xxS_sey]
Ke = 0.48 / (xx - 0.21) + 1.25
Sey = 0.64 * xx -0.06
sey_bpt_oxy = (yy > Ke | xx > 0.1) & yy > Sey

xx = data$logN2Ha[xx.xxS]
yy = data$logO3O2[xx.xxS]
Ke = 0.48 / (xx - 0.21) + 1.25
Sey = 0.64 * xx -0.06
sey_oxy = (yy > Ke | xx > 0.1) & yy > Sey

NS = length(data$logN2Ha[xx.xxS])
NS_bpt = length(data$logN2Ha[xx.xxS_bpt]) 
NS_oii = length(data$logN2Ha[xx.xxS_oii])
NS_bpt_oii = length(data$logN2Ha[xx.xxS_bpt & xx.xxS_oii])

print('----------------- SPIRALS -----------------')
print(sprintf('%-50s%8s%11s%11s', '', 'N', 'N/NS', 'N/NS_d'))
print(sprintf('%-50s%8i', '#N SPIRALS: ', NS))
print(sprintf('%-50s%8i%10.2f%1s', '#N spirals in the BPT diagram: ', NS_bpt, NS_bpt / NS * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s', '#N spirals in the OIII/OII diagram: ', NS_oii, NS_oii / NS * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s', '#N spirals in the BPT AND OIII/OII diagram: ', NS_bpt_oii, NS_bpt_oii / NS * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s%10.2f%1s', '#N Seyfert BPT diagram: ', length(data$logN2Ha[xx.xxS_sey]), 
            length(data$logN2Ha[xx.xxS_sey]) / NS * 100, '%',
            length(data$logN2Ha[xx.xxS_sey]) / NS_bpt * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s%10.2f%1s', '#N Seyfert OIII/OII diagram: ', length(data$logN2Ha[xx.xxS_sey][sey_oxy]), 
            length(data$logN2Ha[xx.xxS_sey][sey_oxy]) / NS * 100, '%',
            length(data$logN2Ha[xx.xxS_sey][sey_oxy]) / NS_oii * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s%10.2f%1s', '#N Seyfert BPT & OIII/OII diagram: ', length(data$logN2Ha[xx.xxS_sey][sey_bpt_oxy]),
              length(data$logN2Ha[xx.xxS_sey][sey_bpt_oxy]) / NS * 100, '%',
              length(data$logN2Ha[xx.xxS_sey][sey_bpt_oxy]) / NS_bpt_oii * 100, '%'))

#points(xx[sey_oxy], yy[sey_oxy], pch = pch1sey, cex = cex1sey, col = col1sey)

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

xx = seq(-4, -0.2, 0.01)
Ka = 0.33 / (xx +0.11) + 1.01
lines(xx, Ka, col = 'red', lwd = 1.5)

xx = seq(-4, 0.1, 0.01)
Ke = 0.48 / (xx - 0.21) + 1.25
lines(xx, Ke, col = 'darkgreen', lwd = 1.5)

xx = seq(-0.35, 5, 0.01)
Sey = 0.64 * xx -0.06
lines(xx, Sey, col = 'black', lwd = 2, lty = 2)

#-----------------------
# ELLIPTICAL GALAXIES
#-----------------------
plot(data$logN2Ha[xx.xxE], data$logO3O2[xx.xxE], xlim = c(-1, 0.5), ylim = c(-1.5, 1.5),
     xlab = expression(paste('log ([NII]/H', alpha, ')')), 
     ylab = expression(paste('log ([OIII]/[OII])')),
     xaxt = 'n', yaxt = 'n', pch = pch2, cex = cex2, col = col2, lwd = 1)

points(data$logN2Ha[xx.xxE_sey], data$logO3O2[xx.xxE_sey], pch = pch2sey, cex = cex2sey, col = col2sey)

xx = data$logN2Ha[xx.xxE_sey]
yy = data$logO3O2[xx.xxE_sey]
Ke = 0.48 / (xx - 0.21) + 1.25
Sey = 0.64 * xx -0.06
sey_bpt_oxy = (yy > Ke | xx > 0.1) & yy > Sey

xx = data$logN2Ha[xx.xxE]
yy = data$logO3O2[xx.xxE]
Ke = 0.48 / (xx - 0.21) + 1.25
Sey = 0.64 * xx -0.06
sey_oxy = (yy > Ke | xx > 0.1) & yy > Sey

NE = length(data$logN2Ha[xx.xxE])
NE_bpt = length(data$logN2Ha[xx.xxE_bpt]) 
NE_oii = length(data$logN2Ha[xx.xxE_oii])
NE_bpt_oii = length(data$logN2Ha[xx.xxE_bpt & xx.xxE_oii])

print('----------------- ELLIPTICALS -----------------')
print(sprintf('%-50s%8s%11s%11s', '', 'N', 'N/NE', 'N/NE_d'))
print(sprintf('%-50s%8i', '#N ELLIPTICALS: ', NE))
print(sprintf('%-50s%8i%10.2f%1s', '#N ellipticals in the BPT diagram: ', NE_bpt, NE_bpt / NE * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s', '#N ellipticals in the OIII/OII diagram: ', NE_oii, NE_oii / NE * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s', '#N ellipticals in the BPT AND OIII/OII diagram: ', NE_bpt_oii, NE_bpt_oii / NE * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s%10.2f%1s', '#N Seyfert BPT diagram: ', length(data$logN2Ha[xx.xxE_sey]), 
              length(data$logN2Ha[xx.xxE_sey]) / NE * 100, '%',
              length(data$logN2Ha[xx.xxE_sey]) / NE_bpt * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s%10.2f%1s', '#N Seyfert OIII/OII diagram: ', length(data$logN2Ha[xx.xxE_sey][sey_oxy]), 
              length(data$logN2Ha[xx.xxE_sey][sey_oxy]) / NE * 100, '%',
              length(data$logN2Ha[xx.xxE_sey][sey_oxy]) / NE_oii * 100, '%'))
print(sprintf('%-50s%8i%10.2f%1s%10.2f%1s', '#N Seyfert BPT & OIII/OII diagram: ', length(data$logN2Ha[xx.xxE_sey][sey_bpt_oxy]),
              length(data$logN2Ha[xx.xxE_sey][sey_bpt_oxy]) / NE * 100, '%',
              length(data$logN2Ha[xx.xxE_sey][sey_bpt_oxy]) / NE_bpt_oii * 100, '%'))

#points(xx[sey_oxy], yy[sey_oxy], pch = pch1sey, cex = cex1sey, col = col1sey)

axis(1, at = seq(-2, 2, 0.5), labels = T, tcl = -0.6)
axis(1, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)
axis(3, at = seq(-2, 2, 0.5), labels = F, tcl = -0.6)
axis(3, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)
pch2sey = 18; cex2sey = 0.3; col2sey = rgb(0.9, 0, 0, 0.8)
axis(2, at = seq(-2, 2, 0.5), labels = T, tcl = -0.6)
axis(2, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)
axis(4, at = seq(-2, 2, 0.5), labels = F, tcl = -0.6)
axis(4, at = seq(-2, 2, 0.1), labels = F, tcl = -0.3)

leg2 = expression(paste('ELLIPTICALS'))
par(font = 2)
legend(-1, 1.5, c(leg2), bty = 'n', pt.cex = c(cex2), pch = c(pch2), 
       col = c(col2), cex = 1.2)
par(font = 1)

xx = seq(-4, -0.2, 0.01)
Ka = 0.33 / (xx +0.11) + 1.01
lines(xx, Ka, col = 'red', lwd = 1.5)

xx = seq(-4, 0.1, 0.01)
Ke = 0.48 / (xx - 0.21) + 1.25
lines(xx, Ke, col = 'darkgreen', lwd = 1.5)

xx = seq(-0.35, 5, 0.01)
Sey = 0.64 * xx -0.06
lines(xx, Sey, col = 'black', lwd = 2, lty = 2)

#legend(-0.75, -1.3, 'Kauffmann+03', bty = 'n', text.col = 'red')
#legend(-0.06, -1.3, 'Kewley+01', bty = 'n', text.col = 'darkgreen')

#par(font = 2)
#legend(-0.35, -1.15, 'Composite', bty = 'n', cex = 1.2)
#par(font = 1)

dev.off()


# -----------------------------------------
composite = data$logO3Hb[xx.xxE] >  0.61 / (data$logN2Ha[xx.xxE] - 0.05) + 1.3

xx = data$logN2Ha[xx.xxE | xx.xxS]
yy = data$logO3Hb[xx.xxE | xx.xxS]
#print(length(xx[is.na(xx) == F & is.na(yy) == F]))
#print(which(is.na(xx) == T & is.na(yy) == T))
