p <- read.table("popdel.ntrio.hwe.txt", h=T, na.strings = "NA")
d <- read.table("delly.ntrio.hwe.txt", h=T, na.strings = "NA")
l <- read.table("lumpy.ntrio.hwe.txt", h=T, na.strings = "NA")

trMendel <- function(t) {
    m <- t[which(is.na(t$NTrio)),]

    tt <- m
    errors <- m[,4]+m[,5]+m[,8]+m[,9]+m[,11]+m[,15]+m[,18]+m[,19]
    tt$correct <- m[,6]+m[,7]+m[,10]+m[,12]+m[,13]+m[,14]+m[,16]+m[,17]+m[,20]
    tt$mendel <- 100*errors/(errors+tt$correct)
    
    m <- t[which(t$NTrio==1),]
    tt$trc <- m[,6]+m[,7]
    tt$tr <- (m[,7]/tt$trc)*100
    
    tt <- tt[which(tt$trc > 100),]
    return(tt)
}

dd <- trMendel(d)
ll <- trMendel(l)
pp <- trMendel(p)

d_gt <- 28 + 1
l_gt <- 78 + 1
p_gt <- 25 + 1
#por_gt <- 25 + 1

red <- '#DF0000' #"red2"
green <- '#40B000' #"green3"
blue <- '#002094' #"blue2"

lwd <- .7
lty <- 1
pch <- 4
dot <- 16
cex <- 0.4
cex2 <- 0.7

transmissionRate <- "Deletion transmitted [%]"
numTransmissions <- "Number of transmissions"
mendelRate <- "Inheritance errors [%]"
numConsistent <- "Consistent deletions per trio"

pdf("tr-mendel.pdf", 7, 2)

par(mar=c(4,4,1,1), mfrow=c(1,3))

plot(0,0,xlim=c(0,1800), ylim=c(0,7), type="n", xlab=numConsistent, ylab=mendelRate, bty="n")
abline(h=0, col="grey", lty=lty, lwd=lwd)
#lines(dd$correct/49, dd$mendel, col=green, lwd=lwd)
points(dd$correct/49, dd$mendel, col=green, pch=dot, cex=cex)
#lines(ll$correct/49, ll$mendel ,col=blue, lwd=lwd)
points(ll$correct/49, ll$mendel, col=blue, pch=dot, cex=cex)
#lines(pp$correct/49, pp$mendel, col=red, lwd=lwd)
points(pp$correct/49, pp$mendel, col=red, pch=dot, cex=cex)
#points(dd$correct[d_gt]/49, dd$mendel[d_gt], col=green, pch=pch, cex=.5)
#points(ll$correct[l_gt]/49, ll$mendel[l_gt], col=blue, pch=pch, cex=.5)
#points(pp$correct[p_gt]/49, pp$mendel[p_gt], col=red, pch=pch, cex=.5)
legend("topleft", legend=c("PopDel", "Delly", "Lumpy"), col=c(red, green, blue), lty=1, bty="n")

plot(0,0,xlim=c(0,4000), ylim=c(0,100), type="n", xlab=numTransmissions, ylab=transmissionRate, bty="n")
abline(h=50, col="grey", lty=lty, lwd=lwd)
#lines(dd$trc, dd$tr, col=green, lwd=lwd)
points(dd$trc, dd$tr, col=green, pch=dot, cex=cex)
#lines(ll$trc, ll$tr ,col=blue, lwd=lwd)
points(ll$trc, ll$tr, col=blue, pch=dot, cex=cex)
#lines(pp$trc, pp$tr, col=red, lwd=lwd)
points(pp$trc, pp$tr, col=red, pch=dot, cex=cex)
#points(dd$trc[d_gt], dd$tr[d_gt], col=green, pch=pch, cex=.5)
#points(ll$trc[l_gt], ll$tr[l_gt], col=blue, pch=pch, cex=.5)
#points(pp$trc[p_gt], pp$tr[p_gt], col=red, pch=pch, cex=.5)

plot(0,0,xlim=c(0,7), ylim=c(0,100), type="n", xlab=mendelRate, ylab=transmissionRate, bty="n")
abline(h=50, col="grey", lty=lty, lwd=lwd)
abline(v=0, col="grey", lty=lty, lwd=lwd)
points(0, 0.5, col="grey", pch=pch, cex=.5)
#lines(dd$mendel, dd$tr, col=green, lwd=lwd)
points(dd$mendel, dd$tr, col=green, pch=dot, cex=cex)
#lines(ll$mendel, ll$tr ,col=blue, lwd=lwd)
points(ll$mendel, ll$tr, col=blue, pch=dot, cex=cex)
#lines(pp$mendel, pp$tr, col=red, lwd=lwd)
points(pp$mendel, pp$tr, col=red, pch=dot, cex=cex)
#points(dd$mendel[d_gt], dd$tr[d_gt], col=green, pch=pch, cex=.5)
#points(ll$mendel[l_gt], ll$tr[l_gt], col=blue, pch=pch, cex=.5)
#points(pp$mendel[p_gt], pp$tr[p_gt], col=red, pch=pch, cex=.5)

plot(0,0,xlim=c(800,1250), ylim=c(0,.5), type="n", xlab=numConsistent, ylab=mendelRate, bty="n")
abline(h=0, col="grey", lty=lty, lwd=lwd)
lines(dd$correct/49, dd$mendel, col=green, lwd=lwd)
points(dd$correct/49, dd$mendel, col=green, pch=dot, cex=cex)
lines(ll$correct/49, ll$mendel ,col=blue, lwd=lwd)
points(ll$correct/49, ll$mendel, col=blue, pch=dot, cex=cex)
lines(pp$correct/49, pp$mendel, col=red, lwd=lwd)
points(pp$correct/49, pp$mendel, col=red, pch=dot, cex=cex)
#points(dd$correct[d_gt]/49, dd$mendel[d_gt], col=green, pch=pch, cex=cex2)
#points(ll$correct[l_gt]/49, ll$mendel[l_gt], col=blue, pch=pch, cex=cex2)
#points(pp$correct[p_gt]/49, pp$mendel[p_gt], col=red, pch=pch, cex=cex2)
legend("topleft", legend=c("PopDel", "Delly", "Lumpy"), col=c(red, green, blue), lty=1, bty="n")

plot(0,0,xlim=c(2400,3600), ylim=c(44,56), type="n", xlab=numTransmissions, ylab=transmissionRate, bty="n")
abline(h=50, col="grey", lty=lty, lwd=lwd)
lines(dd$trc, dd$tr, col=green, lwd=lwd)
points(dd$trc, dd$tr, col=green, pch=dot, cex=cex)
lines(ll$trc, ll$tr ,col=blue, lwd=lwd)
points(ll$trc, ll$tr, col=blue, pch=dot, cex=cex)
lines(pp$trc, pp$tr, col=red, lwd=lwd)
points(pp$trc, pp$tr, col=red, pch=dot, cex=cex)
#points(dd$trc[d_gt], dd$tr[d_gt], col=green, pch=pch, cex=cex2)
#points(ll$trc[l_gt], ll$tr[l_gt], col=blue, pch=pch, cex=cex2)
#points(pp$trc[p_gt], pp$tr[p_gt], col=red, pch=pch, cex=cex2)

plot(0,0,xlim=c(0,.5), ylim=c(44,56), type="n", xlab=mendelRate, ylab=transmissionRate, bty="n")
abline(h=50, col="grey", lty=lty, lwd=lwd)
abline(v=0, col="grey", lty=lty, lwd=lwd)
lines(dd$mendel, dd$tr, col=green, lwd=lwd)
points(dd$mendel, dd$tr, col=green, pch=dot, cex=cex)
lines(ll$mendel, ll$tr ,col=blue, lwd=lwd)
points(ll$mendel, ll$tr, col=blue, pch=dot, cex=cex)
lines(pp$mendel, pp$tr, col=red, lwd=lwd)
points(pp$mendel, pp$tr, col=red, pch=dot, cex=cex)
#points(dd$mendel[d_gt], dd$tr[d_gt], col=green, pch=pch, cex=cex2)
#points(ll$mendel[l_gt], ll$tr[l_gt], col=blue, pch=pch, cex=cex2)
#points(pp$mendel[p_gt], pp$tr[p_gt], col=red, pch=pch, cex=cex2)
#points(0, 50, col="grey", pch=pch, cex=cex2)

dev.off()