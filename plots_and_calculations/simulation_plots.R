library(readr)
library(plotrix)

## Colors for tools
pdcol <- '#DF0000'
mantacol <- '#00A0FF'
lumpycol <- '#002094'
dellycol <- '#40B000'
gridsscol <- '#FFCF00'

## Symbols for points
pdchp <- 1
lumpychp <- 2
dellychp <- 3
mantachp <- 4
gridsschp <- 5

## Scaling
lwd <-2.5
cexlab <- 1.25
cexaxis <- 1.25
cexLegend <- 1.25

##Paths to the different files
popdelProfileTimes="time/profile.time"
popdelCallTimes="time/popdel.time"
popdelBufferTimes="time/popdel_buffer0100.time"
dellyTimes="time/delly.time"
dellyNTimes="time/delly_n.time"
lumpyTimes="time/lumpy.time"
lumpyExpressTimes="time/lumpy_express.time"
splitterTimes="time/splitter.time"
splitsortTimes="time/splitsort.time"
discoTimes="time/disco.time"
discosortTimes="time/discosort.time"
svtyperTimes="time/svtyper.time"
gridssTimes="time/gridss.time"
mantaTimes="time/manta.time"

popdelCountsPath <- "TP-FP-FN/popdel/popdel.bed.counts" ##Create by eval_bed.sh
dellyCountsPath  <- "TP-FP-FN/delly/delly.bed.counts"
lumpyCountsPath  <- "TP-FP-FN/lumpy/lumpy.bed.counts"
mantaCountsPath  <- "TP-FP-FN/manta/manta.bed.counts"
gridssCountsPath <- "TP-FP-FN/gridss/gridss.bed.counts"

##### Loading, pre-processing and formating #######

##Runtime and memory for 10-100 simulated samlples##

#Get profiling running times and mem and add them up for all sample numbers
popdelProfile <- read_delim(popdelProfileTimes,
                            "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
colnames(popdelProfile) <- c("real", "usr", "sys", "mem")
popdelProfileSum <- data.frame(c(1, seq(10,100,10), seq(200,1000,100)))
tmp <- c(popdelProfile$real[1], sum(popdelProfile$real[1:10]), sum(popdelProfile$real[1:20]),
         sum(popdelProfile$real[1:30]), sum(popdelProfile$real[1:40]), sum(popdelProfile$real[1:50]),
         sum(popdelProfile$real[1:60]), sum(popdelProfile$real[1:70]), sum(popdelProfile$real[1:80]),
         sum(popdelProfile$real[1:90]), sum(popdelProfile$real[1:100]), sum(popdelProfile$real[1:200]),
         sum(popdelProfile$real[1:300]), sum(popdelProfile$real[1:400]), sum(popdelProfile$real[1:500]),
         sum(popdelProfile$real[1:600]), sum(popdelProfile$real[1:700]), sum(popdelProfile$real[1:800]),
         sum(popdelProfile$real[1:900]), sum(popdelProfile$real[1:1000]))
popdelProfileSum <- cbind(popdelProfileSum, tmp)
tmp <- c(popdelProfile$usr[1], sum(popdelProfile$usr[1:10]), sum(popdelProfile$usr[1:20]),
         sum(popdelProfile$usr[1:30]), sum(popdelProfile$usr[1:40]), sum(popdelProfile$usr[1:50]),
         sum(popdelProfile$usr[1:60]), sum(popdelProfile$usr[1:70]), sum(popdelProfile$usr[1:80]),
         sum(popdelProfile$usr[1:90]), sum(popdelProfile$usr[1:100]), sum(popdelProfile$usr[1:200]),
         sum(popdelProfile$usr[1:300]), sum(popdelProfile$usr[1:400]), sum(popdelProfile$usr[1:500]),
         sum(popdelProfile$usr[1:600]), sum(popdelProfile$usr[1:700]), sum(popdelProfile$usr[1:800]),
         sum(popdelProfile$usr[1:900]), sum(popdelProfile$usr[1:1000]))
popdelProfileSum <- cbind(popdelProfileSum, tmp)
tmp <- c(popdelProfile$sys[1], sum(popdelProfile$sys[1:10]), sum(popdelProfile$sys[1:20]),
         sum(popdelProfile$sys[1:30]), sum(popdelProfile$sys[1:40]), sum(popdelProfile$sys[1:50]),
         sum(popdelProfile$sys[1:60]), sum(popdelProfile$sys[1:70]), sum(popdelProfile$sys[1:80]),
         sum(popdelProfile$sys[1:90]), sum(popdelProfile$sys[1:100]), sum(popdelProfile$sys[1:200]),
         sum(popdelProfile$sys[1:300]), sum(popdelProfile$sys[1:400]), sum(popdelProfile$sys[1:500]),
         sum(popdelProfile$sys[1:600]), sum(popdelProfile$sys[1:700]), sum(popdelProfile$sys[1:800]),
         sum(popdelProfile$sys[1:900]), sum(popdelProfile$sys[1:1000]))
popdelProfileSum <- cbind(popdelProfileSum, tmp)
tmp <- c(popdelProfile$mem[1], max(popdelProfile$mem[1:10]), max(popdelProfile$mem[1:20]),
         max(popdelProfile$mem[1:30]), max(popdelProfile$mem[1:40]), max(popdelProfile$mem[1:50]),
         max(popdelProfile$mem[1:60]), max(popdelProfile$mem[1:70]), max(popdelProfile$mem[1:80]),
         max(popdelProfile$mem[1:90]), max(popdelProfile$mem[1:100]), max(popdelProfile$mem[1:000]),
         max(popdelProfile$mem[1:300]), max(popdelProfile$mem[1:400]), max(popdelProfile$mem[1:500]),
         max(popdelProfile$mem[1:600]), max(popdelProfile$mem[1:700]), max(popdelProfile$mem[1:800]),
         max(popdelProfile$mem[1:900]), max(popdelProfile$mem[1:1000]))
popdelProfileSum <- cbind(popdelProfileSum, tmp)
rm(tmp)
colnames(popdelProfileSum) <- c("samples", "real", "usr", "sys", "mem")

## Get running times and memory values of remaining tools.
cnames <- c("samples", "real", "usr", "sys", "mem")
popdelCall <- read_delim(popdelCallTimes,
                         "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
delly <- read_delim(dellyTimes,
                    "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
dellyN <- read_delim(dellyNTimes,
                    "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
lumpy <- read_delim(lumpyTimes,
                    "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
lumpyExpress <- read_delim(lumpyExpressTimes,
                          "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
lumpy_split <- read_delim(splitterTimes,
                    "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
lumpy_splitsort <- read_delim(splitsortTimes,
                          "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
lumpy_disco <- read_delim(discoTimes,
                          "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
lumpy_discosort <- read_delim(discosortTimes,
                          "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
svtyper <- read_delim(svtyperTimes,
                      "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
gridss <- read_delim(gridssTimes,
                     "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
manta <- read_delim(mantaTimes,
                    "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)

## Add up running times of lumpy's steps ##
lumpyExpressSum <- data.frame(c(1, seq(10,100,10)))
lumpyExpressSum <- cbind(lumpyExpressSum, lumpyExpress$real + svtyper$real[1:11])
lumpyExpressSum <- cbind(lumpyExpressSum, lumpyExpress$usr + svtyper$usr[1:11])
lumpyExpressSum <- cbind(lumpyExpressSum, lumpyExpress$sys + svtyper$sys[1:11])
lumpyExpressSum <- cbind(lumpyExpressSum, pmax(lumpyExpress$mem, svtyper$mem[1:11]))
colnames(lumpyExpressSum) <- c("samples", "real", "usr", "sys", "mem")

lumpySum <- data.frame(c(1, seq(10,100,10), seq(200,500,100)))
lumpyPreprocess <- data.frame(lumpy_split[1], lumpy_split[2:4] + lumpy_splitsort[2:4] + lumpy_disco[2:4] + lumpy_discosort[2:4],
                              pmax(as.matrix(lumpy_split[5]),
                                   as.matrix(lumpy_splitsort[5]),
                                   as.matrix(lumpy_disco[5]),
                                   as.matrix(lumpy_discosort[5])))
tmp <- c(lumpyPreprocess$real[1], sum(lumpyPreprocess$real[1:10]), sum(lumpyPreprocess$real[1:20]),
         sum(lumpyPreprocess$real[1:30]), sum(lumpyPreprocess$real[1:40]), sum(lumpyPreprocess$real[1:50]),
         sum(lumpyPreprocess$real[1:60]), sum(lumpyPreprocess$real[1:70]), sum(lumpyPreprocess$real[1:80]),
         sum(lumpyPreprocess$real[1:90]), sum(lumpyPreprocess$real[1:100]), sum(lumpyPreprocess$real[1:200]),
         sum(lumpyPreprocess$real[1:300]), sum(lumpyPreprocess$real[1:400]), sum(lumpyPreprocess$real[1:500]))
lumpySum <- cbind(lumpySum, tmp + svtyper$real)

tmp <- c(lumpyPreprocess$usr[1], sum(lumpyPreprocess$usr[1:10]), sum(lumpyPreprocess$usr[1:20]),
         sum(lumpyPreprocess$usr[1:30]), sum(lumpyPreprocess$usr[1:40]), sum(lumpyPreprocess$usr[1:50]),
         sum(lumpyPreprocess$usr[1:60]), sum(lumpyPreprocess$usr[1:70]), sum(lumpyPreprocess$usr[1:80]),
         sum(lumpyPreprocess$usr[1:90]), sum(lumpyPreprocess$usr[1:100]), sum(lumpyPreprocess$usr[1:200]),
         sum(lumpyPreprocess$usr[1:300]), sum(lumpyPreprocess$usr[1:400]), sum(lumpyPreprocess$usr[1:500]))
lumpySum <- cbind(lumpySum, tmp + svtyper$usr)
tmp <- c(lumpyPreprocess$sys[1], sum(lumpyPreprocess$sys[1:10]), sum(lumpyPreprocess$sys[1:20]),
         sum(lumpyPreprocess$sys[1:30]), sum(lumpyPreprocess$sys[1:40]), sum(lumpyPreprocess$sys[1:50]),
         sum(lumpyPreprocess$sys[1:60]), sum(lumpyPreprocess$sys[1:70]), sum(lumpyPreprocess$sys[1:80]),
         sum(lumpyPreprocess$sys[1:90]), sum(lumpyPreprocess$sys[1:100]), sum(lumpyPreprocess$sys[1:200]),
         sum(lumpyPreprocess$sys[1:300]), sum(lumpyPreprocess$sys[1:400]), sum(lumpyPreprocess$sys[1:500]))
lumpySum <- cbind(lumpySum, tmp + svtyper$sys)
tmp <- c(max(lumpyPreprocess$mem[1]), max(lumpyPreprocess$mem[1:10]), max(lumpyPreprocess$mem[1:20]),
         max(lumpyPreprocess$mem[1:30]), max(lumpyPreprocess$mem[1:40]), max(lumpyPreprocess$mem[1:50]),
         max(lumpyPreprocess$mem[1:60]), max(lumpyPreprocess$mem[1:70]), max(lumpyPreprocess$mem[1:80]),
         max(lumpyPreprocess$mem[1:90]), max(lumpyPreprocess$mem[1:100]), max(lumpyPreprocess$mem[1:200]),
         max(lumpyPreprocess$mem[1:300]), max(lumpyPreprocess$mem[1:400]), max(lumpyPreprocess$mem[1:500]))
lumpySum <- cbind(lumpySum, pmax(tmp, lumpy$mem, svtyper$mem))
colnames(lumpySum) <- c("samples", "real", "usr", "sys", "mem")
lumpySum <- cbind(lumpySum, pmax(tmp, lumpy$mem, svtyper$mem))
############ Plots #############


par(mfrow=c(2,2), mar=c(3.5,4.2,0.5,1))

##CPU-Time
f = 3600 ## factor for conversion of seconds into hours
plot(x=popdelCall$samples, y=(popdelCall$usr + popdelCall$sys + popdelProfileSum$usr + popdelProfileSum$sys) / f,
     ylim=c(0,25), type="l", col=NULL,
     ylab="",
     xlab="",
     xlim = c(1,1000),
     xaxt="n",
    #main="CPU-Times for growing number of simulated samples (chr21)",
     lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty='n')
tickpos = c (1, seq(100,1000,100))
axis(side=1, at=tickpos, cex.axis = cexaxis)
mtext(side=1, line=2.5, "Number of samples", cex=cexlab)
mtext(side=2, line=2.5, "CPU-time [h]", cex=cexlab)
lines(x=popdelCall$samples, y=(popdelCall$usr + popdelCall$sys + popdelProfileSum$usr + popdelProfileSum$sys) / f,
      col = pdcol, type="l", lwd=lwd)
lines(x=delly$samples, y=(delly$usr + delly$sys) / f, col = dellycol, type="l", lwd=lwd)
lines(x=lumpySum$samples, y=(lumpySum$usr + lumpySum$sys) / f, col = lumpycol, type="l", lwd=lwd)
lines(x=lumpyExpressSum$samples, y=(lumpyExpressSum$usr + lumpyExpressSum$sys) / f, col = lumpycol, type="l", lty="dashed", lwd=lwd)
lines(x=gridss$samples, y=(gridss$usr + gridss$sys) / f, col = gridsscol, type="l", lwd=lwd)
lines(x=manta$samples, y=(manta$usr + manta$sys) / f, col = mantacol, type="l", lwd=lwd)
lines(x=dellyN$samples, y=(dellyN$usr + dellyN$sys) / f, col = dellycol, type="l", lty = "dashed", lwd=lwd)
points(x=c(manta$samples[nrow(manta)],
           lumpySum$samples[nrow(lumpySum)],
           lumpyExpressSum$samples[nrow(lumpyExpressSum)],
           delly$samples[nrow(delly)],
           dellyN$samples[nrow(dellyN)],
           gridss$samples[nrow(gridss)]),
       y=c((manta$sys[nrow(manta)] + manta$usr[nrow(manta)]) / f,
           (lumpySum$sys[nrow(lumpySum)] + lumpySum$usr[nrow(lumpySum)]) / f,
           (lumpyExpressSum$sys[nrow(lumpyExpressSum)] + lumpyExpressSum$usr[nrow(lumpyExpressSum)]) / f,
           (delly$sys[nrow(delly)] + delly$usr[nrow(delly)]) / f,
           (dellyN$sys[nrow(dellyN)] + dellyN$usr[nrow(dellyN)]) / f,
           (gridss$sys[nrow(gridss)] + gridss$usr[nrow(gridss)]) / f),
       col=c(mantacol, lumpycol, lumpycol, dellycol, dellycol, gridsscol),
       pch = 4, lwd = lwd)

##Memory
plot(x=gridss$samples, y=(gridss$mem) / 1000,
     ylim=c(10,100000), type="l", log="xy",
     xlim=c(1,1000),
     yaxt="n",
     xaxt="n",
     col=NULL,
     xlab="",
     ylab="",
     #main="Memory usage for growing number of simulated samples (chr21)",
     lwd=lwd, cex.lab = cexlab,  bty="n")
xticks <- c(seq(1, 10, 1), seq(20, 100, 10), seq(200, 1000, 100))
axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
axis(side = 1, at = c(1, 10, 100, 1000), cex.axis = cexaxis)
yticks = c(10,100, 1000, 10000, 100000)
labelsY=c("10", expression(10^2), expression(10^3), expression(10^4), expression(10^5))
axis(side=2, at=yticks, cex.axis = cexaxis, labels = labelsY)
mtext(side=1, line=2.5, "Number of samples", cex=cexlab)
mtext(side=2, line=2.5, "Max. memory [MB]", cex=cexlab)
lines(x=popdelCall$samples, y=(pmax(popdelCall$mem, popdelProfileSum$mem)) / 1000, col = pdcol, type = "l", lwd=lwd)
lines(x=delly$samples, y=(delly$mem)/ 1000, col = dellycol, type="l", lwd=lwd)
lines(x=lumpySum$samples, y=(lumpySum$mem) / 1000, col = lumpycol, type="l", lwd=lwd)
lines(x=manta$samples, y=(manta$mem) / 1000, col = mantacol, type="l", lwd=lwd)
lines(x=gridss$samples, y=(gridss$mem) / 1000, col = gridsscol, type="l", lwd=lwd)
lines(x=lumpyExpressSum$samples, y=(lumpyExpressSum$mem) / 1000, col = lumpycol, type="l", lty="dashed" ,lwd=lwd)
lines(x=dellyN$samples, y=(dellyN$mem)/ 1000, col = dellycol, type="l", lty="dashed", lwd=lwd)
points(x=c(manta$samples[nrow(manta)],
           lumpySum$samples[nrow(lumpySum)],
           lumpyExpressSum$samples[nrow(lumpyExpressSum)],
           delly$samples[nrow(delly)],
           dellyN$samples[nrow(dellyN)],
           gridss$samples[nrow(gridss)]),
       y=c(manta$mem[nrow(manta)] / 1000,
           lumpySum$mem[nrow(lumpySum)] / 1000,
           lumpyExpressSum$mem[nrow(lumpyExpressSum)] / 1000,
           delly$mem[nrow(delly)] / 1000,
           dellyN$mem[nrow(dellyN)] / 1000,
           gridss$mem[nrow(gridss)] / 1000),
       col=c(mantacol, lumpycol, lumpycol, dellycol, dellycol, gridsscol),
       pch = 4, lwd = lwd)

#### Calculation of TP, FP, FN ####

############## TP-FP-FN ###############
##Function for extracting precision and recall
getPrecRec <- function(infile)
{
    df <- as.data.frame(t(read.table(infile)))
    colnames(df) <- c("SampleNum", "TP","FP","FN")
    df$TotalVar <- df$TP + df$FN
    df$Prec <- df$TP /  (df$TP + df$FP)
    df$Rec <- df$TP / (df$TP + df$FN)
    df$F1 <- 2/(1 / df$Prec + 1 / df$Rec)
    return(df)
}

## Loading the data
popdelCounts <- getPrecRec(popdelCountsPath)
dellyCounts  <- getPrecRec(dellyCountsPath)
lumpyCounts  <- getPrecRec(lumpyCountsPath)
mantaCounts  <- getPrecRec(mantaCountsPath)
gridssCounts <- getPrecRec(gridssCountsPath)

## Precision
#formatC(x, digits = 1, format = "f")
plot(x=popdelCounts$SampleNum, y=popdelCounts$Prec, type="l", col=NULL, ylim= c(0.7,1), log="x",
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n",
     #main = "Precision for growing number of simulated samples (chr21)",
     lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty = "n", pch=pdchp)
axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
axis(side = 1, at = c(1, 10, 100, 1000), cex.axis = cexaxis)
yticks = seq(0.7, 1, 0.1)
labelsY=c(0.7, 0.8, 0.9, 1)
axis(side=2, at=yticks, cex.axis = cexaxis)
mtext(side=1, line=2.5, "Number of samples", cex=cexlab)
mtext(side=2, line=2.5, "Precision", cex=cexlab)
lines(x=popdelCounts$SampleNum, y=popdelCounts$Prec, type="l", col=pdcol, lwd=lwd, pch=pdchp)
lines(x=dellyCounts$SampleNum, y=dellyCounts$Prec, type="l", col=dellycol, lwd=lwd, pch=dellychp)
lines(x=lumpyCounts$SampleNum, y=lumpyCounts$Prec, type="l", col=lumpycol, lwd=lwd, pch=lumpychp)
lines(x=gridssCounts$SampleNum, y=gridssCounts$Prec, type="l", col=gridsscol, lwd=lwd, pch=gridsschp)
lines(x=mantaCounts$SampleNum, y=mantaCounts$Prec, type="l", col=mantacol, lwd=lwd, pch=mantachp)
points(x=c(mantaCounts$SampleNum[nrow(mantaCounts)],
           lumpyCounts$SampleNum[nrow(lumpyCounts)],
           dellyCounts$SampleNum[nrow(dellyCounts)],
           gridssCounts$SampleNum[nrow(gridssCounts)]),
       y=c(mantaCounts$Prec[nrow(manta)],
           lumpyCounts$Prec[nrow(lumpyCounts)],
           dellyCounts$Prec[nrow(dellyCounts)],
           gridssCounts$Prec[nrow(gridssCounts)]),
       col=c(mantacol, lumpycol, dellycol, gridsscol),
       pch = 4, lwd = lwd)
legend("bottomleft" , inset=0.01, legend=c("PopDel", "Delly", "LUMPY", "Manta", "GRIDSS"),
       col = c(pdcol, dellycol, lumpycol, mantacol, gridsscol),
       cex=cexLegend, lwd=lwd, bty="n")

## Recall
plot(x=popdelCounts$SampleNum, y=popdelCounts$Rec, type="l", col=NULL, ylim= c(0,1), log="x",
     xlab ="",
     ylab = "",
     xaxt = "n",
     #main = "Recall for growing number of simulated samples (chr21)",
     lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty="n")
axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
axis(side = 1, at = c(1, 10, 100, 1000), cex.axis = cexaxis)
mtext(side=1, line=2.5, "Number of samples", cex=cexlab)
mtext(side=2, line=2.5, "Recall", cex=cexlab)
lines(x=popdelCounts$SampleNum, y=popdelCounts$Rec, type="l", col=pdcol, lwd=lwd, pch=pdchp)
lines(x=dellyCounts$SampleNum, y=dellyCounts$Rec, type="l", col=dellycol, lwd=lwd, pch=dellychp)
lines(x=lumpyCounts$SampleNum, y=lumpyCounts$Rec, type="l", col=lumpycol, lwd=lwd, pch=lumpychp)
lines(x=gridssCounts$SampleNum, y=gridssCounts$Rec, type="l", col=gridsscol, lwd=lwd, pch=gridsschp)
lines(x=mantaCounts$SampleNum, y=mantaCounts$Rec, type="l", col=mantacol, lwd=lwd, pch=mantachp)
points(x=c(mantaCounts$SampleNum[nrow(mantaCounts)],
           lumpyCounts$SampleNum[nrow(lumpyCounts)],
           dellyCounts$SampleNum[nrow(dellyCounts)],
           gridssCounts$SampleNum[nrow(gridssCounts)]),
       y=c(mantaCounts$Rec[nrow(manta)],
           lumpyCounts$Rec[nrow(lumpyCounts)],
           dellyCounts$Rec[nrow(dellyCounts)],
           gridssCounts$Rec[nrow(gridssCounts)]),
       col=c(mantacol, lumpycol, dellycol, gridsscol),
       pch = 4, lwd = lwd)


## PopDel plot for varying buffer size

popdelBuffer <- read_delim(popdelBufferTimes,
                          "\t", escape_double = FALSE, trim_ws = TRUE, col_names = c("b", "real", "usr", "sys", "mem"))
par(mfrow=c(1,1), mar=c(4,4.5,1,2))
f = 60 # scaling factor for time
plot(x = 1, y = 0,
     #main = "Running time VS memory consumption (parameter -b)",
     col = NULL,
     log = "x",
     ylim = c(1, 3000 / f),
     ylab = "CPU time [min]",
     xlim = c(10,11000),
     xlab = "Max. memory [MB]",
     bty = "n",
     cex.lab = cexlab,
     xaxt = "n")
grid(0, NULL)
tickpos <- c(seq(10, 100, 10),seq(200, 1000, 100),seq(2000, 10000, 1000))
axis(side = 1, at = tickpos, labels=FALSE)
axis(side = 1, at = c(10, 100, 1000, 10000))
abline (v = tickpos, lty = "dotted", col = "lightgray")
lines(x = popdelBuffer$mem / 1000, y = (popdelBuffer$usr + popdelBuffer$sys) / f, type ="b", pch = pdchp,col = pdcol, lwd = lwd)
defaultIdx = which(popdelBuffer$b == 200000)
points(x = popdelBuffer$mem[defaultIdx] / 1000 , y = (popdelBuffer$usr[defaultIdx] + popdelBuffer$sys[defaultIdx]) / f,
       col = pdcol, pch = 16)
