library(readr)
library(plotrix)

################################Functions#####################################
firstColToNames <- function(df)
{
  tmp <- as.data.frame(df)
  n <- strtoi(tmp[,1], 10)
  out <- tmp[,-1]
  row.names(out) <- n
  return(out)
}
##Function for extracting precision and recall
getPrecRec <- function(path, chromosomes)
{
  df <- as.data.frame(t(read.table(gsub("/", paste("/", chromosomes[1], ".",sep=""), path))))
  batches <- df[,1]
  df <- df[,2:4]
  if (length(chromosomes)  > 1)
  {
    for (c in chromosomes[2:length(chromosomes)])
    {
      tmp <- as.data.frame(t(read.table(gsub("/", paste("/", c, ".",sep=""), path))))
      df <- df + tmp[,2:4]
    }
  }  
  df <- cbind(batches, df)
  colnames(df) <- c("SampleNum", "TP","FP","FN")
  df$TotalVar <- df$TP + df$FN
  df$Prec <- df$TP /  (df$TP + df$FP)
  df$Rec <- df$TP / (df$TP + df$FN)
  df$F1 <- 2/(1 / df$Prec + 1 / df$Rec)
  return(df)
}
appendGtInfo <- function(df, path, chromosomes)
{
  t <- read.table(gsub("bed", "GT", gsub("/", paste("/GT/", chromosomes[1], ".",sep=""), path)))
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
  if (length(chromosomes)  > 1)
  {
    for (c in chromosomes[2:length(chromosomes)])
    {
      tmp <- as.data.frame(read.table(gsub("bed", "GT", gsub("/", paste("/GT/", c, ".",sep=""), path))))
      t[,2:9] <- t[,2:9] + tmp[,2:9]
    }
  }
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
load_all <- function(path, chromosomes)
{
  return(appendGtInfo(getPrecRec(path, chromosomes), path, chromosomes))
}
###################################################################
############### Options, settings and paths #######################

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

##Paths to the different files,
popdelCountsPath <- "popdel/popdel.bed.counts"
dellyCountsPath  <- "delly/delly.bed.counts"
lumpyCountsPath  <- "lumpy/lumpy.bed.counts"
mantaCountsPath  <- "manta/manta.bed.counts"

batches<-c(seq(1,9), seq(10,100,10), seq(200,500,100))
############ Plots #############

############## TP-FP-FN ###############

## Loading the data
chromosomes <- c("chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
popdelCounts <- load_all(popdelCountsPath, chromosomes)
dellyCounts <-  load_all(dellyCountsPath, chromosomes)
lumpyCounts <-  load_all(lumpyCountsPath, chromosomes)
mantaCounts <-  load_all(mantaCountsPath, chromosomes)

## Precision

  pdf("figures/chr17-22.G1ksim.pdf", width=8, height=2.5)
  par(mfrow=c(1,2), mar=c(3.5,4.2,0.5,1))
  
  #formatC(x, digits = 1, format = "f")
  plot(x=popdelCounts$SampleNum, y=popdelCounts$Prec, type="l", col=NULL, ylim= c(0.8,1), log="x",
       xlab = "",
       ylab = "",
       yaxt = "n",
       xaxt = "n",
       #main = "Precision for growing number of simulated samples (chr21)",
       lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty = "n", pch=pdchp)
  xticks <- c(seq(1, 10, 1), seq(20, 100, 10), seq(200, 500, 100))
  axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
  axis(side = 1, at = c(1, 10, 100, 500), cex.axis = cexaxis)
  yticks = seq(0.8, 1, 0.1)
  labelsY=c(0.8, 0.9, 1)
  axis(side=2, at=yticks, cex.axis = cexaxis)
  mtext(side=1, line=2.5, "Number of genomes", cex=cexlab)
  mtext(side=2, line=2.5, "Precision", cex=cexlab)
  lines(x=dellyCounts$SampleNum, y=dellyCounts$Prec, type="l", col=dellycol, lwd=lwd, pch=dellychp)
  lines(x=lumpyCounts$SampleNum, y=lumpyCounts$Prec, type="l", col=lumpycol, lwd=lwd, pch=lumpychp)
  lines(x=mantaCounts$SampleNum, y=mantaCounts$Prec, type="l", col=mantacol, lwd=lwd, pch=mantachp)
  lines(x=popdelCounts$SampleNum, y=popdelCounts$Prec, type="l", col=pdcol, lwd=lwd, pch=pdchp)
  
  legend("bottomleft" , inset=0.01, legend=c("PopDel", "Delly", "Lumpy", "Manta"),
         col = c(pdcol, dellycol, lumpycol, mantacol),
         cex=cexLegend, lwd=lwd, bty="n")
  
  ## Recall
  plot(x=popdelCounts$SampleNum, y=popdelCounts$Rec, type="l", col=NULL, ylim= c(0.6,1.0), log="x",
       xlab ="",
       ylab = "",
       xaxt = "n",
       #main = "Recall for growing number of simulated genomes (chr21)",  
       lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty="n")
  axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
  mtext(side=1, line=2.5, "Number of genomes", cex=cexlab)
  axis(side = 1, at = c(1, 10, 100, 500), cex.axis = cexaxis)
  mtext(side=2, line=2.5, "Recall", cex=cexlab)
  lines(x=popdelCounts$SampleNum, y=popdelCounts$Rec, type="l", col=pdcol, lwd=lwd, pch=pdchp)
  lines(x=dellyCounts$SampleNum, y=dellyCounts$Rec, type="l", col=dellycol, lwd=lwd, pch=dellychp)
  lines(x=lumpyCounts$SampleNum, y=lumpyCounts$Rec, type="l", col=lumpycol, lwd=lwd, pch=lumpychp)
  lines(x=mantaCounts$SampleNum, y=mantaCounts$Rec, type="l", col=mantacol, lwd=lwd, pch=mantachp)
  
  dev.off()
  
  
## Calculate GT statistics
getUndercall <- function(t)
{
  return(t$FN1 + 2*t$FN2 + t$FN3)
}
getOvercall <- function(t)
{
  return(t$FP1 + 2*t$FP2 + t$FP3 )
}
getMisscall <- function(t)
{
  return(getUndercall(t) + getOvercall(t))
}
getTotalCalls <- function(t)
{
  return(2 * (t$TP1 + t$TP2 + t$FN1 + t$FN2 + t$FN3 + t$FP1 + t$FP2))
}
getGTstats <- function(t)
{
  df <- as.data.frame(cbind(t$SampleNum, getUndercall(t), getOvercall(t), getMisscall(t), getTotalCalls(t)))
  colnames(df) <- c("SampleNum", "undercall", "overcall", "misscall", "total")
  return(df)
}

popdelGTstats <- getGTstats(popdelCounts)
dellyGTstats <- getGTstats(dellyCounts)
lumpyGTstats <- getGTstats(lumpyCounts)
mantaGTstats <- getGTstats(mantaCounts)

lwd<-1.25
cex<-1
cexLegend<-1
cexlab<-1  
pdf("G1k_GT.pdf", width=6, height=4)
par(mfrow=c(1,1), mar=c(3.5,4,0.5,0.5))
plot(x=popdelGTstats$SampleNum, y=1 - (popdelGTstats$misscall / popdelGTstats$total), type="l", col=NULL, ylim= c(0.7, 1), log="x",
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n",
     lwd=lwd, cex.axis = cexaxis, cex.lab = cexlab, bty = "n", pch=pdchp)
axis(side = 1, at = xticks, labels=FALSE, cex.axis = cexaxis)
mtext(side=1, line=2.5, "Number of genomes", cex=cexlab)
axis(side = 1, at = c(1, 10, 100, 500), cex.axis = cexaxis)
yticks = seq(0.7, 1, 0.1)
axis(side=2, at=yticks, cex.axis = cexaxis)
mtext(side=2, line=2.5, "Ratio of correct alleles", cex=cexlab)
lines(x=popdelGTstats$SampleNum,y=1 - (popdelGTstats$misscall / popdelGTstats$total), type="l", col=pdcol, lwd=lwd, pch=pdchp)
lines(x=dellyGTstats$SampleNum, y=1 - (dellyGTstats$misscall / dellyGTstats$total), type="l", col=dellycol, lwd=lwd, pch=dellychp)
lines(x=lumpyGTstats$SampleNum, y=1 - (lumpyGTstats$misscall / lumpyGTstats$total), type="l", col=lumpycol, lwd=lwd, pch=lumpychp)
lines(x=mantaGTstats$SampleNum, y=1 - (mantaGTstats$misscall / mantaGTstats$total), type="l", col=mantacol, lwd=lwd, pch=mantachp)
legend("bottomright" , inset=0.01, legend=c("PopDel", "Delly", "Lumpy", "Manta"),
       col = c(pdcol, dellycol, lumpycol, mantacol),
       cex=cexLegend, lwd=lwd, bty="n")
dev.off()

## Create tables for output
getTable <- function(t)
{
  table <- as.data.frame(cbind(t$SampleNum, t$TP, t$FP, t$FN, t$Prec, t$Rec, t$F1, t$precGT, t$recGT, t$F1GT))
  colnames(table) <- c("Samples", "TP", "FP", "FN", "Prec.", "Rec.", "F1", "GT-Prec.", "GT-Rec.", "GT-F1")
  return(table)
}
popdelTable <- getTable(popdelCounts)
dellyTable <- getTable(dellyCounts)
mantaTable <- getTable(mantaCounts)
lumpyTable <- getTable(lumpyCounts)

TPTable <- cbind(popdelTable$Samples, popdelTable$TP, dellyTable$TP, mantaTable$TP, lumpyTable$TP)
FPTable <- cbind(popdelTable$Samples, popdelTable$FP, dellyTable$FP, mantaTable$FP, lumpyTable$FP)
FNTable <- cbind(popdelTable$Samples, popdelTable$FN, dellyTable$FN, mantaTable$FN, lumpyTable$FN)
PrecTable <- cbind(popdelTable$Samples, popdelTable$Prec., dellyTable$Prec., mantaTable$Prec., lumpyTable$Prec.)
RecTable <- cbind(popdelTable$Samples, popdelTable$Rec., dellyTable$Rec., mantaTable$Rec., lumpyTable$Rec.)
F1Table <- cbind(popdelTable$Samples, popdelTable$F1, dellyTable$F1, mantaTable$F1, lumpyTable$F1)
GTPrecTable <- cbind(popdelTable$Samples, popdelTable$"GT-Prec.", dellyTable$"GT-Prec.", mantaTable$"GT-Prec.", lumpyTable$"GT-Prec.")
GTRecTable <- cbind(popdelTable$Samples, popdelTable$"GT-Rec.", dellyTable$"GT-Rec.", mantaTable$"GT-Rec.", lumpyTable$"GT-Rec.")
GTF1Table <- cbind(popdelTable$Samples, popdelTable$"GT-F1", dellyTable$"GT-F1", mantaTable$"GT-F1", lumpyTable$"GT-F1")
tools <- c("Samples", "PopDel", "Delly", "Manta", "Lumpy")
colnames(TPTable) <- tools
colnames(FPTable) <- tools
colnames(FNTable) <- tools
colnames(PrecTable) <- tools
colnames(RecTable) <- tools
colnames(F1Table) <- tools
colnames(GTPrecTable) <- tools
colnames(GTRecTable) <- tools
colnames(GTF1Table) <- tools

write.table(TPTable, "supplement/G1kSim.TP.table.txt", row.names=FALSE)
write.table(FPTable, "supplement/G1kSim.FP.table.txt", row.names=FALSE)
write.table(FNTable, "supplement/G1kSim.FN.table.txt", row.names=FALSE)
write.table(PrecTable, "supplement/G1kSim.Prec.table.txt", row.names=FALSE)
write.table(RecTable, "supplement/G1kSim.Rec.table.txt", row.names=FALSE)
write.table(F1Table, "supplement/G1kSim.F1.table.txt", row.names=FALSE)
write.table(GTPrecTable, "supplement/G1kSim.GT-Prec.table.txt", row.names=FALSE)
write.table(GTRecTable, "supplement/G1kSim.GT-Rec.table.txt", row.names=FALSE)
write.table(GTF1Table, "supplement/G1kSim.GT-F1.table.txt", row.names=FALSE)

