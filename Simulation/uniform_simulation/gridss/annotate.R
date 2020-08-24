#!/usr/bin/env Rscript

## Adapted from: https://github.com/PapenfussLab/gridss/blob/master/example/simple-event-annotation.R

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
    stop("Input and output must be provided.", call.=FALSE)
}

library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
#' Simple SV type classifier
simpleEventType <- function(gr) {
	return(ifelse(seqnames(gr) != seqnames(partner(gr)), "CTX", # inter-chromosomosal
								ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
											 ifelse(strand(gr) == strand(partner(gr)), "INV",
											 			 ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
											 			 			 "DUP")))))
}

vcf <- readVcf(args[1])
if (length(vcf) > 0)
{
    info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
	    row.names=c("SIMPLE_TYPE"),
	    Number=c("1"),
	    Type=c("String"),
	    Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))
    gr <- breakpointRanges(vcf)
    svtype <- simpleEventType(gr)
    info(vcf)$SIMPLE_TYPE <- NA_character_
    info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
    info(vcf[gr$sourceId])$SVLEN <- gr$svLen
} else
{
    print("Input VCF does not contain any records. Writing empty output VCF anyway.")
}
writeVcf(vcf, args[2])
