library(readr)
cexLegend <- 1.25

calls <- as.data.frame(t(read_delim("pca/popdel.calls", "\t", escape_double = FALSE, trim_ws = TRUE)))
#calls <- as.data.frame(t(read_delim("pca/delly.calls", "\t", escape_double = FALSE, trim_ws = TRUE)))
#calls <- as.data.frame(t(read_delim("pca/lumpy.calls", "\t", escape_double = FALSE, trim_ws = TRUE)))

perSampleVar <- as.data.frame(read.table("popdel.NoCentromeres.perSampleVariants.GTfilter.txt"))
tmp <- perSampleVar[,1]
perSampleVar <- as.data.frame(perSampleVar[,2])
row.names(perSampleVar) <- tmp

groups <- as.data.frame(read_delim("pca/ancestry.csv", "\t", escape_double = FALSE, col_names = c("sample", "ancestry"), trim_ws = TRUE))
tmp <- groups$sample
groups <- as.data.frame(groups$ancestry)
row.names(groups) <- tmp
colnames(groups) <- c("ancestry")

perSampleVar <- cbind(perSampleVar, groups[rownames(perSampleVar),])
colnames(perSampleVar) <- c("variants", "ancestry")

# Check for duplicate rows/samples
anyDuplicated(calls) # Returns false, so everything is fine.

totalVar <- ncol(calls)
# Get rid of duplicate columns/deletions to remove complete linkage and/or uninformative columns
tcalls <- t(calls)
noDups <- !duplicated(tcalls)
calls <- calls[which(noDups)]
rm(tcalls)
dupVar <- totalVar - ncol(calls)

# Get rid of variants that always have the same cygocity.
vars <- ncol(calls)
calls <- calls[, colSums(calls) != 0]
calls <- calls[,apply(calls,2,function(calls) !all(calls==1))]
calls <- calls[,apply(calls,2,function(calls) !all(calls==2))]
constVar <- vars - ncol(calls)

## Find significant correlations (~LD) of neighboring variants and remove one of the pairs.
iterations <- 0
vars <- ncol(calls)
print("Removing variants in LD...")
while (TRUE)
{
    cors <- rep(TRUE, ncol(calls))
    for (i in 2:ncol(calls))
    {
        if ((cor.test(calls[,i-1], calls[,i], method = "spearman"))$p.value <= 0.05)
            cors[i] <- FALSE
    }
    if (any(!cors))
    {
        iterations <- iterations + 1
        calls <- calls[which(cors)]
    }
    else
        break
}
links <- vars - ncol(calls)
print(paste0("Removed ", totalVar - ncol(calls), " variants in total (", (totalVar - ncol(calls)) / totalVar * 100, "%). ", ncol(calls), " variants remaining."))
print(paste0(dupVar, " variants where perfect duplicates of other variants."))
print(paste0(constVar, " variants had the same genotype in all samples (compared to themselves)."))
print(paste0(links, " variants in LD were removed in ", iterations, " iterations."))

# Shuffle by rows and columns, just to be safe
calls <- calls[sample(nrow(calls)),]
calls <- calls[,sample(ncol(calls))]

#Order the groups so that they correspond to the order in calls
groups <- groups[rownames(calls),]

# Normalize all deletions to have mean 0 and equal variance. Do PCA.
calls.pca <- prcomp(calls, center = TRUE, scale = TRUE)

# Plot the first two principle components and colors the samples by ancestry.
colors <- c("red2", "green3", "blue2")

summary(calls.pca) ## To get the variance of the first two PC's

#pdf("pca.pdf", width=3, height=3)
par(mfrow=c(1,1), mar=c(3.5,4,0.5,0.5))
plot(calls.pca$x[,1:2], col = colors[groups], bty="n", frame=FALSE, pch = 20, cex=0.7, xlim=c(-30,30), cex.axis = 0.66,
     ylim=c(-20,20), xlab = "" , ylab = "", mgp = c(3, 0.5, 0))

mtext(side=1, line=1.5, "PC1 (??% explained variance)", cex = 0.75) ## IMPORTANT UPDATE explained var when changing something!
mtext(side=2, line=1.5, "PC2 (??% explained variance)", cex = 0.75)
legend("bottomright" , inset=0.01, legend=c("AFR", "EAS", "EUR"),
       col = colors, pch = 20, cex=0.8, lwd=-1, bty="n", lty="47")
#dev.off()

#pdf("boxplot.pdf", width=2.5, height=3)
boxplot(perSampleVar$variants ~ perSampleVar$ancestry, notch=TRUE, ylab = "", bty = "n", frame=FALSE,
       boxwex = 1/4,
       outline = FALSE,
       ylim=c(1300,2600),
       lwd = 1.5,
       cex.lab =0.75,
       cex.axis = 0.66,
       par(mgp = c(3, 0.5, 0)))
mtext(side=2, line=1.5, "Deletions per sample", cex = 0.75)
stripchart(perSampleVar$variants ~ perSampleVar$ancestry, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = colors, cex=0.3)
#dev.off()


#Count HET and HOM variants
hetHom <- read_delim("pca/popdel.NoCentromeres.GTfilter.hetHom.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
mean(hetHom$HET_PL27)
mean(hetHom$HOM_PL27)
