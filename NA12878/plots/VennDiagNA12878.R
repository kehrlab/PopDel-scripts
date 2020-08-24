library(VennDiagram)

#t <- read.table("eval/counts.personalis.500-10000.0.5.tsv")     ##select one
t <- read.table("eval/counts.pacbio.500-10000.0.5.tsv")     ##select one


c <- as.data.frame(t(t$V1))
names(c) <- t$V2
rm(t)

#tname <- "Illumina"  ##Adapt according to the reference set
tname <- "PacBio"  ##Adapt according to the reference set

truthcol <- "blue"
pdcol <- "red"
dellycol <- "green" 
smoovecol <- "yellow"

par(mar = rep(0.1, 4), mfrow=c(1,1))
plot.new()
draw.triple.venn(c$truth, c$popdel, c$delly, c$truth_popdel, c$popdel_delly, c$truth_delly, c$truth_popdel_delly, fill = c(truthcol, pdcol, dellycol), c(tname, "PopDel", "Delly"))
plot.new()
draw.triple.venn(c$truth, c$popdel, c$smoove, c$truth_popdel, c$popdel_smoove, c$truth_smoove, c$truth_popdel_smoove, fill = c(truthcol, pdcol, smoovecol), c(tname, "PopDel", "smoove"))
plot.new()
draw.triple.venn(c$truth, c$popdel, c$manta, c$truth_popdel, c$popdel_manta, c$truth_manta, c$truth_popdel_manta, fill = c(truthcol, pdcol, mantacol), c(tname, "PopDel", "Manta"))
