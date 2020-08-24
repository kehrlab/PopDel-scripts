library(readr)
#library(plotrix)

################################Functions#####################################
## Create the sum of the given batches for the desired column
sumBatches <- function(source, batches, col)
{
    out <- vector(mode="numeric")
    for (i in batches)
    {
        out <- c(out, sum(source[1:i, col]))
    }
    return(out)
}
## Take the max of every batch
maxBatches <- function(source, batches, col)
{
    out <- vector(mode="numeric")
    for (i in batches)
    {
        out <- c(out, max(source[1:i, col]))
    }
    return(out)
}
#Sum the times for all single processes into batches of 1,2,...,10,20,...,100,200...,1000 samples.
sumFromSingles <- function(infile, cnames)
{
    times <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
    
    if (ncol(times) == 5) ## Smoove's times have the sample number as first column. Ignore them.
        times <- times[,2:5]
    
    batches <- c(seq(1,10), seq(20,100,10), seq(200,1000,100))
    timeSum <- data.frame(batches)
    for (i in seq(1,3))
        timeSum <- cbind(timeSum, sumBatches(times, batches, i))
    
    timeSum <- cbind(timeSum, maxBatches(times, batches, 4))
    
    colnames(timeSum) <- cnames
    return(timeSum)
}
readSmooveGenotypeBatches <- function(infile, batches)
{
    times <- read_delim(infile, "\t", escape_double = FALSE, trim_ws = TRUE,
                        col_names = c("batch", "sample", "real", "usr", "sys", "mem"))

    c <- vector(mode="numeric") 
    m <- vector(mode="numeric") 
    for (b in batches) 
    {
        if (b > 1)
        {
            c <- rbind(c, colSums(times[strtoi(times$batch, 10)==b & strtoi(times$sample, 10)<=b,
                                                   c(3,4,5)]))
            m <- rbind(m, max(times[strtoi(times$batch, 10)==b & strtoi(times$sample, 10)<=b, 6]))
        }
        else
        {   ## This process was not performed for single sample calling
            c <- rbind(c, c(0,0,0))
            m <- 0
        }
    }
    out <- cbind(c,m)
    out <- as.data.frame(out)
    colnames(out) <- c("real", "usr", "sys", "mem")
    rownames(out) <- batches
    return(out)
}
firstColToNames <- function(df)
{
    tmp <- as.data.frame(df)
    n <- strtoi(tmp[,1], 10)
    out <- tmp[,-1]
    row.names(out) <- n
    return(out)
}
getSmooveAll <- function(smooveCallTimes, smooveMergeTimes, smooveGenotypeTimes, smoovePasteTimes, batches, cnames)
{
    smooveCall <- firstColToNames(sumFromSingles(smooveCallTimes, cnames))
    smooveCall[1,] <- read.table(smooveSingleTimes)[1,-1] ## Replace the time for one sample with the correct time for single sample calling
    smooveMerge <- rbind(c(0,0,0,0), firstColToNames(read_delim(smooveMergeTimes,
                                                                "\t",
                                                                escape_double = FALSE,
                                                                trim_ws = TRUE,
                                                                col_names = cnames)))
    smooveGenotype <- readSmooveGenotypeBatches(smooveGenotypeTimes, batches)
    smoovePaste <- rbind(c(0,0,0,0), firstColToNames(read_delim(smoovePasteTimes,
                                                                "\t",
                                                                escape_double = FALSE,
                                                                trim_ws = TRUE,
                                                                col_names = cnames)))
    
    ## Add up running times of smoove's steps ##
    out <- smooveCall[,-4] + smooveMerge[,-4] + smooveGenotype[,-4] + smoovePaste[,-4]
    out <- cbind(batches, out, pmax(smooveCall$mem, smooveMerge$mem, smooveGenotype$mem, smoovePaste$mem))
    colnames(out) <- cnames
    return(out)
}
appendGtInfo <- function(df, infile)
{
    t <- read.table(infile)
    df$FP1 <- 0
    df$FP2 <- 0
    df$FP3 <- 0
    df$FN1 <- 0
    df$FN2 <- 0
    df$FN3 <- 0
    df$TP1 <- 0
    df$TP2 <- 0
    df$TotalFPGT <- 0
    df$TotalFNGT <- 0
    df$TotalTGT <- 0
    df$TotalFGT <- 0
    df$precGT <- 0
    df$recGT <- 0
    df$F1GT <- 0
    for (j in seq(1, dim(t)[1]))
    {
        n <- as.numeric(t[j,1])
        i <- which(df$SampleNum==n)
        df[i,]$FP1 <- as.numeric(t[j,2])
        df[i,]$FP2 <- as.numeric(t[j,5])
        df[i,]$FP3 <- as.numeric(t[j,9])
        df[i,]$FN1 <- as.numeric(t[j,7])
        df[i,]$FN2 <- as.numeric(t[j,8])
        df[i,]$FN3 <- as.numeric(t[j,6])
        df[i,]$TP1 <- as.numeric(t[j,4])
        df[i,]$TP2 <- as.numeric(t[j,3])
        df[i,]$TotalFPGT <-  df[i,]$FP1 + df[i,]$FP2 +  df[i,]$FP3
        df[i,]$TotalFNGT <- df[i,]$FN1 + df[i,]$FN2 + df[i,]$FN3
        df[i,]$TotalTGT <- df[i,]$TP1 + df[i,]$TP2
        df[i,]$TotalFGT <- df[i,]$TotalFPGT + df[i,]$TotalFNGT
        df[i,]$precGT <- df[i,]$TotalTGT / (df[i,]$TotalTGT + df[i,]$TotalFPGT)
        df[i,]$recGT <- df[i,]$TotalTGT / (df[i,]$TotalTGT + df[i,]$TotalFNGT)
        df[i,]$F1GT <- 2/(1 / df[i,]$precGT + 1 / df[i,]$recGT)
    }
    return(df)
}

###################################################################
############### Options, settings and paths #######################

## Colors for tools
pdcol <- '#DF0000'
mantacol <- '#00A0FF'
smoovecol <- '#002094'
dellycol <- '#40B000'
gridsscol <- '#FFCF00'

## Symbols for points
pdchp <- 1
smoovechp <- 2
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
#dellyTimes="time/delly.time"
dellyNTimes="time/delly_n.time"
smooveCallTimes="time/smoove_call.time"
smooveSingleTimes="time/smoove_single.time"
smooveMergeTimes="time/smoove_merge.time"
smooveGenotypeTimes="time/smoove_genotype.time"
smoovePasteTimes="time/smoove_paste.time"
gridssTimes="time/gridss.time"
mantaTimes="time/manta.time"

popdelCountsPath <- "eval/popdel/paper/popdel.bed.counts" ##Create by eval_bed.sh
dellyCountsPath  <- "eval/delly_n/delly_n.bed.counts"
smooveCountsPath  <- "eval/smoove/smoove.bed.counts"
mantaCountsPath  <- "eval/manta/manta.bed.counts"
gridssCountsPath <- "eval/gridss/gridss.bed.counts"

popdelGTCountsPath <- "eval/popdel/paper/popdel.GT.counts" ##Alternative with genotypes. created by eval.sh
dellyGTCountsPath  <- "eval/delly_n/delly_n.GT.counts"
smooveGTCountsPath  <- "eval/smoove/smoove.GT.counts"
mantaGTCountsPath  <- "eval/manta/manta.GT.counts"
# 
# popdelCountsPath <- "TP-FP-FN/popdel/popdel.py.counts" ##Alternative w/o genotypes. Created by eval.sh
# dellyCountsPath  <- "TP-FP-FN/delly_n/delly_n.py.counts"
# smooveCountsPath  <- "TP-FP-FN/smoove/smoove.py.counts"
# mantaCountsPath  <- "TP-FP-FN/manta/manta.py.counts"
# gridssCountsPath <- "TP-FP-FN/gridss/gridss.py.counts"


############ Loading, pre-processing and formating ############

batches<-c(seq(1,9), seq(10,100,10), seq(200,1000,100))
cnames <- c("samples", "real", "usr", "sys", "mem")
popdelProfile <- sumFromSingles(popdelProfileTimes, cnames)
popdelCall <- read_delim(popdelCallTimes, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
dellyN <- read_delim(dellyNTimes, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
gridss <- read_delim(gridssTimes, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
manta <- read_delim(mantaTimes, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = cnames)
smoove <- getSmooveAll(smooveCallTimes, smooveMergeTimes, smooveGenotypeTimes, smoovePasteTimes, batches, cnames)

############ Plots #############

pdf("figures/figure1.pdf", width=8, height=5)
par(mfrow=c(2,2), mar=c(3.5,4.2,0.5,1))

##CPU-Time
f = 3600 ## factor for conversion of seconds into hours
plot(x=popdelCall$samples, y=(popdelCall$usr + popdelCall$sys + popdelProfile$usr + popdelProfile$sys) / f,
     ylim=c(0,80), type="l", col=NULL,
     ylab="",
     xlab="",
     xlim = c(1,1000),
     xaxt="n",
    #main="CPU-Times for growing number of simulated samples (chr21)",
     lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty='n')
tickpos = c (1, seq(100,1000,100))
axis(side=1, at=tickpos, cex.axis = cexaxis)
mtext(side=1, line=2.5, "Number of genomes", cex=cexlab)
mtext(side=2, line=2.5, "CPU-time [h]", cex=cexlab)
lines(x=popdelCall$samples, y=(popdelCall$usr + popdelCall$sys + popdelProfile$usr + popdelProfile$sys) / f,
      col = pdcol, type="l", lwd=lwd)
#lines(x=delly$samples, y=(delly$usr + delly$sys) / f, col = dellycol, type="l", lwd=lwd)
lines(x=smoove$samples, y=(smoove$usr + smoove$sys) / f, col = smoovecol, type="l", lwd=lwd)
lines(x=gridss$samples, y=(gridss$usr + gridss$sys) / f, col = gridsscol, type="l", lwd=lwd)
lines(x=manta$samples, y=(manta$usr + manta$sys) / f, col = mantacol, type="l", lwd=lwd)
lines(x=dellyN$samples, y=(dellyN$usr + dellyN$sys) / f, col = dellycol, type="l", lwd=lwd)


##Memory
plot(x=gridss$samples, y=(gridss$mem) / 1000,  ##check smoove memory at 400 samples. Check smoove running times
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
mtext(side=1, line=2.5, "Number of genomes", cex=cexlab)
mtext(side=2, line=2.5, "Max. memory [MB]", cex=cexlab)
lines(x=popdelCall$samples, y=(pmax(popdelCall$mem, popdelProfile$mem)) / 1000, col = pdcol, type = "l", lwd=lwd)
#lines(x=delly$samples, y=(delly$mem)/ 1000, col = dellycol, type="l", lwd=lwd)
lines(x=smoove$samples, y=(smoove$mem) / 1000, col = smoovecol, type="l", lwd=lwd)
lines(x=manta$samples, y=(manta$mem) / 1000, col = mantacol, type="l", lwd=lwd)
lines(x=gridss$samples, y=(gridss$mem) / 1000, col = gridsscol, type="l", lwd=lwd)
lines(x=dellyN$samples, y=(dellyN$mem)/ 1000, col = dellycol, type="l", lwd=lwd)


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
smooveCounts  <- getPrecRec(smooveCountsPath)
mantaCounts  <- getPrecRec(mantaCountsPath)
gridssCounts <- getPrecRec(gridssCountsPath)

## Precision
#formatC(x, digits = 1, format = "f")
plot(x=popdelCounts$SampleNum, y=popdelCounts$Prec, type="l", col=NULL, ylim= c(0.8,1), log="x",
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n",
     #main = "Precision for growing number of simulated samples (chr21)",
     lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty = "n", pch=pdchp)
axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
axis(side = 1, at = c(1, 10, 100, 1000), cex.axis = cexaxis)
yticks = seq(0.8, 1, 0.1)
labelsY=c(0.8, 0.9, 1)
axis(side=2, at=yticks, cex.axis = cexaxis)
mtext(side=1, line=2.5, "Number of genomes", cex=cexlab)
mtext(side=2, line=2.5, "Precision", cex=cexlab)
lines(x=dellyCounts$SampleNum, y=dellyCounts$Prec, type="l", col=dellycol, lwd=lwd, pch=dellychp)
lines(x=smooveCounts$SampleNum, y=smooveCounts$Prec, type="l", col=smoovecol, lwd=lwd, pch=smoovechp)
lines(x=gridssCounts$SampleNum, y=gridssCounts$Prec, type="l", col=gridsscol, lwd=lwd, pch=gridsschp)
lines(x=mantaCounts$SampleNum, y=mantaCounts$Prec, type="l", col=mantacol, lwd=lwd, pch=mantachp)
lines(x=popdelCounts$SampleNum, y=popdelCounts$Prec, type="l", col=pdcol, lwd=lwd, pch=pdchp)
points(x=gridssCounts$SampleNum[nrow(gridssCounts)],
       y=gridssCounts$Prec[nrow(gridssCounts)],
       col=gridsscol,
       pch = 4, lwd = lwd)
legend("bottomleft" , inset=0.01, legend=c("PopDel", "Delly", "Lumpy", "Manta", "GRIDSS"),
       col = c(pdcol, dellycol, smoovecol, mantacol, gridsscol),
       cex=cexLegend, lwd=lwd, bty="n")

## Recall
plot(x=popdelCounts$SampleNum, y=popdelCounts$Rec, type="l", col=NULL, ylim= c(0.6,1), log="x",
     xlab ="",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     #main = "Recall for growing number of simulated samples (chr21)",
     lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty="n")
axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
axis(side = 1, at = c(1, 10, 100, 1000), cex.axis = cexaxis)
yticks = seq(0.6, 1, 0.1)
labelsY=c(0.6, 0.8, 1)
axis(side=2, at=yticks, labels=FALSE, cex.axis = cexaxis)
axis(side = 2, at = labelsY, cex.axis = cexaxis)
mtext(side=1, line=2.5, "Number of genomes", cex=cexlab)
mtext(side=2, line=2.5, "Recall", cex=cexlab)
lines(x=popdelCounts$SampleNum, y=popdelCounts$Rec, type="l", col=pdcol, lwd=lwd, pch=pdchp)
lines(x=dellyCounts$SampleNum, y=dellyCounts$Rec, type="l", col=dellycol, lwd=lwd, pch=dellychp)
lines(x=smooveCounts$SampleNum, y=smooveCounts$Rec, type="l", col=smoovecol, lwd=lwd, pch=smoovechp)
lines(x=gridssCounts$SampleNum, y=gridssCounts$Rec, type="l", col=gridsscol, lwd=lwd, pch=gridsschp)
lines(x=mantaCounts$SampleNum, y=mantaCounts$Rec, type="l", col=mantacol, lwd=lwd, pch=mantachp)
points(x=gridssCounts$SampleNum[nrow(gridssCounts)],
       y=gridssCounts$Rec[nrow(gridssCounts)],
       col=gridsscol,
       pch = 4, lwd = lwd)
dev.off()

## With GT consideration

popdelCounts <- getPrecRec(popdelCountsPath)
popdelCounts <- appendGtInfo(popdelCounts, popdelGTCountsPath)
dellyCounts  <- getPrecRec(dellyCountsPath)
dellyCounts  <- appendGtInfo(dellyCounts, dellyGTCountsPath)
smooveCounts  <- getPrecRec(smooveCountsPath)
smooveCounts  <- appendGtInfo(smooveCounts, smooveGTCountsPath)
mantaCounts  <- getPrecRec(mantaCountsPath)
mantaCounts  <- appendGtInfo(mantaCounts, mantaGTCountsPath)
gridssCounts <- getPrecRec(gridssCountsPath)

popdelUndercall <- popdelCounts$FN1 + 2*popdelCounts$FN2 + popdelCounts$FN3
popdelOvercall <- popdelCounts$FP1 + 2*popdelCounts$FP2 + popdelCounts$FP3 
popdelMisscall <- popdelUndercall + popdelOvercall
popdelTotal <- 2 * (popdelCounts$TP1 + popdelCounts$TP2 + popdelCounts$FN1 + popdelCounts$FN2 + popdelCounts$FN3 + popdelCounts$FP1 + popdelCounts$FP2)

dellyUndercall <- dellyCounts$FN1 + 2*dellyCounts$FN2 + dellyCounts$FN3
dellyOvercall <- dellyCounts$FP1 + 2*dellyCounts$FP2 + dellyCounts$FP3 
dellyMisscall <- dellyUndercall + dellyOvercall
dellyTotal <- 2 * (dellyCounts$TP1 + dellyCounts$TP2 + dellyCounts$FN1 + dellyCounts$FN2 + dellyCounts$FN3 + dellyCounts$FP1 + dellyCounts$FP2)

smooveUndercall <- smooveCounts$FN1 + 2*smooveCounts$FN2 + smooveCounts$FN3
smooveOvercall <- smooveCounts$FP1 + 2*smooveCounts$FP2 + smooveCounts$FP3 
smooveMisscall <- smooveUndercall + smooveOvercall
smooveTotal <- 2 * (smooveCounts$TP1 + smooveCounts$TP2 + smooveCounts$FN1 + smooveCounts$FN2 + smooveCounts$FN3 + smooveCounts$FP1 + smooveCounts$FP2)

mantaUndercall <- mantaCounts$FN1 + 2*mantaCounts$FN2 + mantaCounts$FN3
mantaOvercall <- mantaCounts$FP1 + 2*mantaCounts$FP2 + mantaCounts$FP3 
mantaMisscall <- mantaUndercall + mantaOvercall
mantaTotal <- 2 * (mantaCounts$TP1 + mantaCounts$TP2 + mantaCounts$FN1 + mantaCounts$FN2 + mantaCounts$FN3 + mantaCounts$FP1 + mantaCounts$FP2)

lwd<-1.25
cex<-1
cexLegend<-1
cexlab<-1
pdf("uniform_GT.pdf", width=6, height=4)
par(mfrow=c(1,1), mar=c(3.5,4,0.5,0.5))
plot(x=1, y=1, type="l", col=NULL, ylim= c(0.7, 1), xlim=c(1,1000), log="x",
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n",
     lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty = "n", pch=pdchp)
axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
axis(side = 1, at = c(1, 10, 100, 1000), cex.axis = cexaxis)
yticks = seq(0.7, 1, 0.1)
axis(side=2, at=yticks, cex.axis = cexaxis)
mtext(side=1, line=2.5, "Number of genomes", cex=cexlab)
mtext(side=2, line=2.5, "Ratio of correct alleles", cex=cexlab)
lines(x=popdelCounts$SampleNum,y=1 - (popdelMisscall / popdelTotal), type="l", col=pdcol, lwd=lwd, pch=pdchp)
lines(x=dellyCounts$SampleNum, y=1 - (dellyMisscall / dellyTotal), type="l", col=dellycol, lwd=lwd, pch=dellychp)
lines(x=smooveCounts$SampleNum, y=1 - (smooveMisscall / smooveTotal), type="l", col=smoovecol, lwd=lwd, pch=smoovechp)
lines(x=mantaCounts$SampleNum, y=1 - (mantaMisscall / mantaTotal), type="l", col=mantacol, lwd=lwd, pch=mantachp)
legend("bottomright" , inset=0.01, legend=c("PopDel", "Delly", "Lumpy", "Manta"),
       col = c(pdcol, dellycol, smoovecol, mantacol),
       cex=cexLegend, lwd=lwd, bty="n")
dev.off()

##Write data as tables
getTable <- function(counts)
{
    table <- as.data.frame(cbind(counts$SampleNum, counts$TP, counts$FP, counts$FN ,counts$Prec, counts$Rec, counts$F1, counts$precGT, counts$recGT, counts$F1GT))
    colnames(table) <- c("Samples", "TP", "FP", "FN", "Prec.", "Rec.", "F1", "GT-Prec.", "GT-Rec.", "GT-F1")
    return(table)
}

popdelTable <- getTable(popdelCounts)
dellyTable <- getTable(dellyCounts)
mantaTable <- getTable(mantaCounts)
smooveTable <- getTable(smooveCounts)
gridssTable <- as.data.frame(cbind(gridssCounts$SampleNum, gridssCounts$TP, gridssCounts$FP, gridssCounts$FN ,gridssCounts$Prec, gridssCounts$Rec, gridssCounts$F1))
colnames(gridssTable) <- c("Samples", "TP", "FP", "FN", "Prec.", "Rec.", "F1")

#write.table(popdelTable, "supplement/popdel.table.txt", row.names=FALSE)
#write.table(dellyTable, "supplement/delly.table.txt", row.names=FALSE)
#write.table(mantaTable, "supplement/manta.table.txt", row.names=FALSE)
#write.table(smooveTable, "supplement/smoove.table.txt", row.names=FALSE)
#write.table(gridssTable, "supplement/gridss.table.txt", row.names=FALSE)

#Fill gaps of GRIDSS table
for (i in seq(500, 1000, 100))
{
    gridssTable <- rbind(gridssTable, c(i, NA, NA, NA, NA, NA, NA))    
}

TPTable <- cbind(popdelTable$Samples, popdelTable$TP, dellyTable$TP, mantaTable$TP, smooveTable$TP, gridssTable$TP)
FPTable <- cbind(popdelTable$Samples, popdelTable$FP, dellyTable$FP, mantaTable$FP, smooveTable$FP, gridssTable$FP)
FNTable <- cbind(popdelTable$Samples, popdelTable$FN, dellyTable$FN, mantaTable$FN, smooveTable$FN, gridssTable$FN)
PrecTable <- cbind(popdelTable$Samples, popdelTable$Prec., dellyTable$Prec., mantaTable$Prec., smooveTable$Prec., gridssTable$Prec.)
RecTable <- cbind(popdelTable$Samples, popdelTable$Rec., dellyTable$Rec., mantaTable$Rec., smooveTable$Rec., gridssTable$Rec.)
F1Table <- cbind(popdelTable$Samples, popdelTable$F1, dellyTable$F1, mantaTable$F1, smooveTable$F1, gridssTable$F1)
GTPrecTable <- cbind(popdelTable$Samples, popdelTable$"GT-Prec.", dellyTable$"GT-Prec.", mantaTable$"GT-Prec.", smooveTable$"GT-Prec.", rep(NA, 28))
GTRecTable <- cbind(popdelTable$Samples, popdelTable$"GT-Rec.", dellyTable$"GT-Rec.", mantaTable$"GT-Rec.", smooveTable$"GT-Rec.", rep(NA, 28))
GTF1Table <- cbind(popdelTable$Samples, popdelTable$"GT-F1", dellyTable$"GT-F1", mantaTable$"GT-F1", smooveTable$"GT-F1", rep(NA, 28))
tools <- c("Samples", "PopDel", "Delly", "Manta", "Lumpy", "GRIDSS")
colnames(TPTable) <- tools
colnames(FPTable) <- tools
colnames(FNTable) <- tools
colnames(PrecTable) <- tools
colnames(RecTable) <- tools
colnames(F1Table) <- tools
colnames(GTPrecTable) <- tools
colnames(GTRecTable) <- tools
colnames(GTF1Table) <- tools

write.table(TPTable, "supplement/TP.table.txt", row.names=FALSE)
write.table(FPTable, "supplement/FP.table.txt", row.names=FALSE)
write.table(FNTable, "supplement/FN.table.txt", row.names=FALSE)
write.table(PrecTable, "supplement/Prec.table.txt", row.names=FALSE)
write.table(RecTable, "supplement/Rec.table.txt", row.names=FALSE)
write.table(F1Table, "supplement/F1.table.txt", row.names=FALSE)
write.table(GTPrecTable, "supplement/GT-Prec.table.txt", row.names=FALSE)
write.table(GTRecTable, "supplement/GT-Rec.table.txt", row.names=FALSE)
write.table(GTF1Table, "supplement/GT-F1.table.txt", row.names=FALSE)
