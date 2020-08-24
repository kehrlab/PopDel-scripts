library(readr)

## Colors for tools
pdcol <- '#DF0000'
smoovecol <- '#002094'
dellycol <- '#40B000'
mantacol <- '#00A0FF'

## Scaling
lwd <- .7
lty <- 1
pch <- 4
dot <- 16
cex <- 0.4
cex2 <- 0.7

### Without 0-0-0 sites in trios.
popdel <-read_delim("eval/mendel/popdel_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
delly <- read_delim("eval/mendel/delly_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
smoove <- read_delim("eval/mendel/smoove_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
manta <- read_delim("eval/mendel/manta_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)

popdel_hc <-read_delim("eval/mendel/popdel_highConf_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
delly_hc <- read_delim("eval/mendel/delly_highConf_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
smoove_hc <- read_delim("eval/mendel/smoove_highConf_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)
manta_hc <- read_delim("eval/mendel/manta_highConf_Mendel_Q.txt", " ", escape_double = FALSE, trim_ws = TRUE)


pdf("mendel_Ashkenazim.pdf", width=3, height=3)
par(mar = c(3.5,4,0.5,0.5), mfrow=c(1,1))
plot(0,0, ylim=c(0,0.03), xlim=c(0, 2000), col=NULL, bty="n", xaxt = "n", yaxt="n", ylab="", xlab="")
xticks <- seq(0, 2000, 500)
axis(side = 1, at = xticks, cex.axis = cexaxis)
yticks = seq(0, 0.03, 0.01)
labelsY=c("0", "1", "2", "3")
axis(side=2, at=yticks, cex.axis = cexaxis, labels=labelsY)
mtext(side=1, line=2, "Consistent deletions per trio", cex=cexlab)
mtext(side=2, line=2, "Inheritance errors [%]", cex=cexlab)
abline(h=0, col="grey", lty=1)
lines(x=popdel$consistent, y=1-(popdel$relative/100), col=pdcol, lwd=lwd, lty=lty)
lines(x=delly$consistent, y=1-(delly$relative/100), col=dellycol, lwd=lwd, lty=lty)
lines(x=smoove$consistent, y=1-(smoove$relative/100),  col=smoovecol, lwd=lwd, lty=lty)
lines(x=manta$consistent, y=1-(manta$relative/100), col=mantacol, lwd=lwd, lty=lty)
points(x=popdel$consistent, y=1-(popdel$relative/100), col=pdcol, pch=dot, cex=cex)
points(x=delly$consistent, y=1-(delly$relative/100), col=dellycol, pch=dot, cex=cex)
points(x=smoove$consistent, y=1-(smoove$relative/100), col=smoovecol, pch=dot, cex=cex)
points(x=manta$consistent, y=1-(manta$relative/100), col=mantacol, pch=dot, cex=cex)
legend("left" , inset=0.01, legend=c("PopDel", "Delly", "Lumpy", "Manta"),
       col = c(pdcol, dellycol, smoovecol, mantacol),
       cex=cexaxis, lwd=lwd, bty="n")
dev.off()

pdf("mendel_highConf_Ashkenazim.pdf", width=3, height=3)
par(mar = c(3.5,4,0.5,0.5), mfrow=c(1,1))
plot(0,0, ylim=c(0,0.02), xlim=c(0, 2000), col=NULL, bty="n", xaxt = "n", yaxt="n", ylab="", xlab="")
xticks <- seq(0, 2000, 500)
axis(side = 1, at = xticks, cex.axis = cexaxis)
yticks = seq(0, 0.03, 0.01)
labelsY=c("0", "1", "2", "3")
axis(side=2, at=yticks, cex.axis = cexaxis, labels=labelsY)
mtext(side=1, line=2, "Consistent deletions per trio", cex=cexlab)
mtext(side=2, line=2, "Inheritance errors [%]", cex=cexlab)
abline(h=0, col="grey", lty=1)
lines(x=popdel_hc$consistent, y=1-(popdel_hc$relative/100), col=pdcol, lwd=lwd, lty=lty)
lines(x=delly_hc$consistent, y=1-(delly_hc$relative/100), col=dellycol, lwd=lwd, lty=lty)
lines(x=smoove_hc$consistent, y=1-(smoove_hc$relative/100),  col=smoovecol, lwd=lwd, lty=lty)
lines(x=manta_hc$consistent, y=1-(manta_hc$relative/100), col=mantacol, lwd=lwd, lty=lty)
points(x=popdel_hc$consistent, y=1-(popdel_hc$relative/100), col=pdcol, pch=dot, cex=cex)
points(x=delly_hc$consistent, y=1-(delly_hc$relative/100), col=dellycol, pch=dot, cex=cex)
points(x=smoove_hc$consistent, y=1-(smoove_hc$relative/100), col=smoovecol, pch=dot, cex=cex)
points(x=manta_hc$consistent, y=1-(manta_hc$relative/100), col=mantacol, pch=dot, cex=cex)
legend("left" , inset=0.01, legend=c("PopDel", "Delly", "Lumpy", "Manta"),
       col = c(pdcol, dellycol, smoovecol, mantacol),
       cex=cexaxis, lwd=lwd, bty="n")
dev.off()
