#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
 
parser$add_argument("vcf",
                    help="Input vcf file")
parser$add_argument("output",
                    help="Output file")

args <- parser$parse_args()

suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(BSgenome.Scerevisiae.UCSC.sacCer3))
suppressPackageStartupMessages(library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene))

vcf <- readVcf(args$vcf, "sacCer3")

sc <- BSgenome.Scerevisiae.UCSC.sacCer3
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

coding <- data.frame(predictCoding(vcf, txdb, seqSource=sc))
coding$ALT = sapply(coding$ALT, as.character)
write.table(coding, args$output, sep='\t', row.names=F, quote=F)
