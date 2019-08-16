library(VennDiagram)

t <- read.table("calls/bed/counts.personalis.500-10000_noCentromeres.0.5.rigorous.tsv")     ##select one
#t <- read.table("calls/bed/counts.personalis.500-10000_noCentromeres.0.8.rigorous.tsv")
#t <- read.table("calls/bed/counts.personalis.500-10000_noCentromeres.0.001.rigorous.tsv")
#t <- read.table("calls/bed/counts.pacbio.500-10000_noCentromeres.0.5.rigorous.tsv")
#t <- read.table("calls/bed/counts.pacbio.500-10000_noCentromeres.0.8.rigorous.tsv")
#t <- read.table("calls/bed/counts.pacbio.500-10000_noCentromeres.0.001.rigorous.tsv")

c <- as.data.frame(t(t$V1))
names(c) <- t$V2
rm(t)

tname <- "PacBio"

names <- c(tname, "PopDel", "Delly", "Lumpy")
cols <- rainbow(4)

par(mar = rep(0.1, 4), mfrow=c(1,1))
plot.new()
draw.pairwise.venn(c$truth, c$popdel, c$truth_popdel, category = c(tname, "PopDel"), fill = cols[1:2],
                   alpha = 0.25, euler.d=TRUE, scaled=TRUE)
plot.new()
draw.pairwise.venn(c$truth, c$delly, c$truth_delly, category = c(tname, "Delly"), fill = c(cols[1], cols[3]),
                   alpha = 0.25, euler.d=TRUE, scaled=TRUE)
par(mar = rep(0.5, 4))
plot.new()
draw.pairwise.venn(c$truth, c$lumpy, c$truth_lumpy, category = c(tname, "Lumpy"), fill = c(cols[1], cols[4]),
                   alpha = 0.25, euler.d=TRUE, scaled=TRUE)

plot.new()
draw.quad.venn(c$truth, c$popdel, c$delly, c$lumpy, c$truth_popdel, c$truth_delly, c$truth_lumpy, c$popdel_delly,
               c$popdel_lumpy, c$delly_lumpy, c$truth_popdel_delly, c$truth_popdel_lumpy, c$truth_delly_lumpy, 
               c$popdel_delly_lumpy, c$truth_popdel_delly_lumpy, category = names, fill = cols, alpha = 0.3,
               euler.d=TRUE, scaled=TRUE)