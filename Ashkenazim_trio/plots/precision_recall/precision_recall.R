library(readr)

## Neccessary tables are created by precision_recall.sh ###
truth <- 678
manta <- read_delim("manta_PR.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
delly <- read_delim("delly_PR.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
lumpy <- read_delim("smoove_PR.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
popdel<- read_delim("popdel_PR.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

manta_single <- read_delim("manta_PR.HG002_single.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
delly_single <- read_delim("delly_PR.HG002_single.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
lumpy_single <- read_delim("smoove_PR.HG002_single.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
popdel_single <- read_delim("popdel_PR.HG002_single.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

## Scaling
lwd <-2

pdcol <- '#DF0000'
mantacol <- '#00A0FF'
lumpycol <- '#002094'
dellycol <- '#40B000'
  
pdf("precision-recall.pdf", width=6, height=2.5)

par(mfrow=c(1,1), mar=c(4.0,4.0,0.5,0.5))
plot(y=1-(popdel$FP/(popdel$FP+popdel$TP)), x=popdel$TP/truth, type="l", ylab="Precision", xlab="Recall", col = pdcol,
     ylim = c(0.8,1), xlim=c(0.4,1), bty="n", lwd=lwd)
lines(y=1-(manta$FP/(manta$FP+manta$TP)), x=manta$TP/truth, col = mantacol, lwd=lwd)
lines(y=1-(delly$FP/(delly$FP+delly$TP)), x=delly$TP/truth, col = dellycol, lwd=lwd)
lines(y=1-(lumpy$FP/(lumpy$FP+lumpy$TP)), x=lumpy$TP/truth, col = lumpycol, lwd=lwd)

lines(y=1-(popdel_single$FP/(popdel_single$FP+popdel_single$TP)), x=popdel_single$TP/truth, col = pdcol, lty = 2, lwd=lwd)
lines(y=1-(manta_single$FP/(manta_single$FP+manta_single$TP)), x=manta_single$TP/truth, col = mantacol, lty = 2, lwd=lwd)
lines(y=1-(delly_single$FP/(delly_single$FP+delly_single$TP)), x=delly_single$TP/truth, col = dellycol, lty = 2, lwd=lwd)
lines(y=1-(lumpy_single$FP/(lumpy_single$FP+lumpy_single$TP)), x=lumpy_single$TP/truth, col = lumpycol, lty = 2, lwd=lwd)

legend(x = 0.38, y = 0.98, inset=0.01, legend=c("PopDel", "Delly", "Manta", "Lumpy"),
       col = c(pdcol, dellycol, mantacol, lumpycol),
       lwd=lwd, bty="n")

legend(x = 0.54, y = 0.92 , inset=0.1, legend=c("Single", "Trio"),
       col = "black", lty = c(2,1),
       lwd=lwd, bty="n")
  
dev.off()
  