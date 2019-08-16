library(readr)

## Colors for tools
pdcol <- '#DF0000'
lumpycol <- '#002094'
dellycol <- '#40B000'

## Scaling
lwd <- 1.5
cexlab <- 0.66
cexaxis <- 0.66

### Without 0-0-0 sites in trios.
popdel <-read_delim("Polaris_kids/popdel_no3Re_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
popdel_deDup <-read_delim("Polaris_kids/popdel_no3Re_Mendel_Q.deDup.txt", " ", escape_double = FALSE, trim_ws = TRUE)
delly <- read_delim("Polaris_kids/delly_no3Re_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
delly_deDup <- read_delim("Polaris_kids/delly_no3Re_deDup_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
lumpy <- read_delim("Polaris_kids/lumpy_no3Re_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
lumpy_deDup <- read_delim("Polaris_kids/lumpy_no3Re_Mendel_Q.deDup.txt", " ", escape_double = FALSE, trim_ws = TRUE)

#pdf("mendel.pdf", width=3, height=3)
par(mar = c(3.5,4,0.5,0.5), mfrow=c(1,1))
plot(x=popdel$consistent/49, y=1-(popdel$relative/100), ylim=c(0,0.08), xlim=c(0, 4200), col=NULL, lwd=lwd, type="l", bty="n", xaxt = "n", yaxt="n", ylab="", xlab="")
lines(x=lumpy$consistent/49, y=1-(lumpy$relative/100),  col=lumpycol, lwd=2, type="l")
lines(x=lumpy_deDup$consistent/49, y=1-(lumpy_deDup$relative/100), col=lumpycol, lwd=lwd, type="l", lty=2)
lines(x=popdel$consistent/49, y=1-(popdel$relative/100), col=pdcol, lwd=lwd, type="l")
lines(x=popdel_deDup$consistent/49, y=1-(popdel_deDup$relative/100), col=pdcol, lty=2, lwd=lwd, type="l")
lines(x=delly$consistent/49, y=1-(delly$relative/100), col=dellycol, lwd=2, type="l")
lines(x=delly_deDup$consistent/49, y=1-(delly_deDup$relative/100), col=dellycol, lty=2, lwd=lwd, type="l")
abline(h=0.3/100, col="grey", lty=3)
xticks <- seq(0, 4000, 1000)
axis(side = 1, at = xticks, cex.axis = cexaxis)
yticks = seq(0, 0.08, 0.01)
labelsY=c("0", "1", "2", "3", "4", "5", "6", "7", "8")
axis(side=2, at=yticks, cex.axis = cexaxis, labels=labelsY)
mtext(side=1, line=1.5, "Consistent deletion sites per trio", cex=cexlab)
mtext(side=2, line=1.5, "Mendelian error rate [%]", cex=cexlab)
legend("left" , inset=0.01, legend=c("PopDel", "Delly", "LUMPY"),
       col = c(pdcol, dellycol, lumpycol),
       cex=cexaxis, lwd=lwd, bty="n")
#dev.off()


## Alternative with consideration of 0-0-0 sites in trios.
#popdel <-read_delim("Polaris_kids/popdel_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
#popdel_deDup <-read_delim("Polaris_kids/popdel_deDup_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
#delly <- read_delim("Polaris_kids/delly_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
#delly_deDup <- read_delim("Polaris_kids/delly_deDup_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
#lumpy <- read_delim("Polaris_kids/lumpy_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
#lumpy_deDup <- read_delim("Polaris_kids/lumpy_deDup_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)

##pdf("mendel.pdf", width=3, height=3)
#par(mar = c(3.5,4,0.5,0.5), mfrow=c(1,1))
#plot(x=popdel$consistent/49, y=1-(popdel$consistent/popdel$total), ylim=c(0,0.02), xlim=c(0, 16000), col=pdcol, lwd=lwd, type="l", bty="n", xaxt = "n", yaxt="n", ylab="", xlab="")
#lines(x=popdel_deDup$consistent/49, y=1-(popdel_deDup$consistent/popdel_deDup$total), col=pdcol, lty=2, lwd=lwd, type="l")
#lines(x=delly$consistent/49, y=1-(delly$consistent/delly$total), col=dellycol, lwd=2, type="l")
#lines(x=delly_deDup$consistent/49, y=1-(delly_deDup$consistent/delly_deDup$total), col=dellycol, lty=2, lwd=lwd, type="l")
#lines(x=lumpy$consistent/49, y=1-(lumpy$consistent/lumpy$total),  col=lumpycol, lwd=2, type="l")
#lines(x=lumpy_deDup$consistent/49, y=1-(lumpy_deDup$consistent/lumpy_deDup$total), col=lumpycol, lwd=lwd, type="l", lty=2)
#abline(h=0.03/100, col="grey", lty=3)
#xticks <- seq(0, 16000, 4000)
#labelsX <-seq(0,16000,4000)
#axis(side = 1, at = xticks, labels=labelsX, cex.axis = cexaxis)
#yticks = seq(0,0.02, 0.005)
#labelsY=c("0.0", "0.5", "1.0", "1.5", "2.0")
#axis(side=2, at=yticks, cex.axis = cexaxis, labels = labelsY)
#mtext(side=1, line=1.5, "Consistent deletion sites per trio", cex=cexlab)
#mtext(side=2, line=1.5, "Mendelian error rate [%]", cex=cexlab)
#legend("left" , inset=0.01, legend=c("PopDel", "Delly", "LUMPY"),
#       col = c(pdcol, dellycol, lumpycol),
#       cex=cexlab, lwd=lwd, bty="n")
##dev.off()
