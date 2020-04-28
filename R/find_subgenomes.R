#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)


if (length(args)!=1) {
  stop("Exactly one argument must be supplied (read lengths txt file).n", call.=FALSE)

}

#libraries
library(ggplot2)


src = c("~/github/npTranscript/R" )
source(.findFile(src, "fxns.find_subgenomes.R"))


#constants
n_bins <- 200
#constants for thresholding algorithm
lag=4
threshold=3
influence=0



length_data <- read.csv(args[1],  header=F, col.names = c('Read_Length'))


pdf('leader_read_hist.pdf')
bins <- hist(length_data$Read_Length, breaks=n_bins, plot=FALSE)
plot(bins, main = 'Length histogram of reads containing 5\' leader sequence', ylab='Count', xlab='Read Length')
invisible(dev.off())

#smoothed z-score algorithm, and flatten non-1 signals
zs <- ThresholdingAlgo(bins$density, lag=lag,threshold=threshold, influence=influence)
flat_scores <- flatten_counts(bins, zs)
find_peaks(flat_scores , bins) #will print out bins


#OLD METHOD. deprecated
#thresh_bins <- bins$breaks[which(ThresholdingAlgo(bins$density, lag=lag,threshold=threshold, influence=influence)$signals ==1)]
#sg_bins <- bin_thresholds(thresh_bins) #output as list
