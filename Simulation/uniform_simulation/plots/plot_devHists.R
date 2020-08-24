pdcol <- '#DF0000'
mantacol <- '#00A0FF'
smoovecol <- '#002094'
dellycol <- '#40B000'
gridsscol <- '#FFCF00'

file_popDel_s0010 <- "eval/popdel/startDev10.hist.tsv"
file_delly_s0010 <- "eval/delly_n/startDev10.hist.tsv"
file_smoove_s0010 <- "eval/smoove/startDev10.hist.tsv"
file_manta_s0010 <- "eval/manta/startDev10.hist.tsv"

file_popDel_s0010 <- "eval/popdel/truth_norm.hist.tsv"
file_delly_s0010 <- "eval/delly_n/truth_norm.hist.tsv"
file_smoove_s0010 <- "eval/smoove/truth_norm.hist.tsv"
file_manta_s0010 <- "eval/manta/truth_norm.hist.tsv"
  
popDel_s0010 <- read.table(file_popDel_s0010, header=TRUE, row.names=NULL, quote="", comment.char="")
delly_s0010 <- read.table(file_delly_s0010, header=TRUE, row.names=NULL, quote="", comment.char="")
smoove_s0010 <- read.table(file_smoove_s0010, header=TRUE, row.names=NULL, quote="", comment.char="")
manta_s0010 <- read.table(file_manta_s0010, header=TRUE, row.names=NULL, quote="", comment.char="")

### Start deviation ###
pdhist <- hist(popDel_s0010$startDev, breaks=seq(-150,150), plot=FALSE)
dellyhist <- hist(delly_s0010$startDev, breaks=seq(-150,150), plot=FALSE)
smoovehist <- hist(smoove_s0010$startDev, breaks=seq(-150,150), plot=FALSE)
mantahist <- hist(manta_s0010$startDev, breaks=seq(-150,150), plot=FALSE)

pdf("figures/startDevHist.pdf", width=8, height=5)
par(mfrow=c(2,2), mar=c(3.5,4.2,0.75,1))
xlim=c(-5,4)
xlab="Start position deviation [bp]"
ylab="Frequency"
cexaxis=1.25
plot(pdhist, bty="n", xaxt="n", cex.axis=cexaxis, col=pdcol, xlim=xlim, xlab="", ylab="", main="")
axis(side=1, pos=0, at=seq(-5,5)-0.5, labels=seq(-5,5))
mtext(side=1, xlab, line=2)
mtext(side=2, ylab, line=2.5)
legend("left" , inset=0.01, legend=c("PopDel", "Delly", "Lumpy", "Manta"),
       col = c(pdcol, dellycol, smoovecol, mantacol),
       cex=cexLegend, lwd=lwd, bty="n")
plot(dellyhist, bty="n", xaxt="n", cex.axis=cexaxis, col=dellycol, xlim=xlim, xlab="", ylab="", main="")
axis(side=1, pos=0, at=seq(-5,5)-0.5, labels=seq(-5,5))
mtext(side=1, xlab, line=2)
mtext(side=2, ylab, line=2.5)
plot(smoovehist, bty="n", xaxt="n", cex.axis=cexaxis, col=smoovecol, xlim=xlim, xlab="", ylab="", main="")
axis(side=1, pos=0, at=seq(-5,5)-0.5, labels=seq(-5,5))
mtext(side=1, xlab, line=2)
mtext(side=2, ylab, line=2.5)
plot(mantahist,bty="n", xaxt="n", cex.axis=cexaxis, col=mantacol, xlim=xlim, xlab="", ylab="", main="")
axis(side=1, pos=0, at=seq(-5,5)-0.5, labels=seq(-5,5))
mtext(side=1, xlab, line=2)
mtext(side=2, ylab, line=2.5)
dev.off()