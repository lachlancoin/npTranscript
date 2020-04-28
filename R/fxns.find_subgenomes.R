##Smoothed z-score algorithm
##	Taken from https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/54507329#54507329
ThresholdingAlgo <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}


###
flatten_counts <- function(hist_out, threshold_out) {

flattened_counts <- vector(mode = 'numeric', length = length(hist_out$counts))

for (i in 1:length(flattened_counts)) {
	flattened_counts[i] <- ifelse(threshold_out$signals[i] ==1, hist_out$counts[i], 0)
}
	return(flattened_counts)
}


find_peaks <- function(flattened_counts, hist_out) {
	n_bins = 0
	for (i in 1:length(flattened_counts)) {
		prior <- flattened_counts[i -1]
		current <- flattened_counts[i]
		follow <- flattened_counts[i+1]
		change <- loc_maxima(prior, current, follow)
		
		if ( isTRUE(change) ) {
			cat(prior, follow, '\n') #print out read lengths of i-1 and i+1 to catch reads at peak within 50bp either side
			n_bins=n_bins+1
		}
		}
	message(sprintf('found %d bins', length(n_bins)))
	}

loc_maxima <- function(x, y, z) {

grad_1 = y-x
grad_2 = z-y
result = ifelse(grad_1 * grad_2 < 0 && grad_1> grad_2, TRUE, FALSE )
return (result)
}



#OLD METHOD (DEPRECATED)
# i=1
# bin_limits <- list()

# prev_len = 0
# cont_bin = FALSE
# current_bin = c(0,0)
#get length bins
# for (len in sig_lengths) {
	# if (len - prev_len == bin_size) {
		# if (cont_bin == FALSE) {
			# cont_bin = TRUE
			# current_bin = c(prev_len, len)
			# } else if (cont_bin == TRUE) {
			# current_bin[2] = len
			# }
		# } 
		# else if (cont_bin == TRUE) {
		# cont_bin = FALSE
		# bin_limits[[i]] <- current_bin
		# cat(current_bin,'\n')
		# i=i+1
		# remove(current_bin)
		# }
	# prev_len = len
	# }
# if (cont_bin == TRUE) {
# bin_limits[[i]] <- current_bin
# cat(current_bin,'\n')
# i=i+1
# }
# message(sprintf('found %d bins', length(bin_limits)))
# return(bin_limits)
# }
