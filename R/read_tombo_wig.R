#Author: Daniel Rawlinson

#add this into R with 'source(read_tombo_wig.R)' then run the function with 'read_wig(your file)'
#TO plot for a specific transcript, run 'plot_meth(name of data, your transcript of interest)'


setwd('C:\\Users\\DRAWLINSON\\Documents\\Manuscripts\\early_evolution_and_genomics')
wig_files <- list.files(pattern='wig')

library(ggplot2)


#read in wig file from tombo. each file is from one experiment, and each experiment can contain multiple transcripts
read_wig <- function(wig_file) {
	wig_data <- list()
	meth_list <- list()
	#open for counting lines
	open_wig <- file(wig_file, 'r')
	i=0
	f_size <- length(readLines(open_wig))
	close(open_wig)
	pb <- txtProgressBar(min=i, max = f_size, initial = 0, title = "Reading methylation Data...")
	#reopen for actual reading
	open_wig <- file(wig_file, 'r')
	while ( TRUE ) {
		line = readLines(open_wig, n=1)
		if ( length(line) == 0 ) {
		  #save last list
		  current <- as.data.frame(do.call(rbind, wig_data))
		  colnames(current) <- c('Base','Fraction')
		  current$Base <- as.integer(as.character(current$Base))
		  current$Fraction <- as.numeric(as.character(current$Fraction))
		  meth_list[[unlist(strsplit(transcript,';')) [2]]] <- current	
			message('Finished!')
			break
			}
		else if ( startsWith(line, 'variableStep') ) {
		  if (exists('transcript')) {
		    current <- as.data.frame(do.call(rbind, wig_data))
		    colnames(current) <- c('Base','Fraction')
		    current$Base <- as.integer(as.character(current$Base))
		    current$Fraction <- as.numeric(as.character(current$Fraction))
		    meth_list[[unlist(strsplit(transcript,';')) [2]]] <- current
		    wig_data <- list() #remove all read data
		  }
			transcript <- gsub('chrom=', '',strsplit(line, ' ')[[1]][2])
			getTxtProgressBar(pb)
			}
		else if (startsWith(line, 'track')) {
			next
			}
		else { 
			i=i+1
			wig_data[[i]] <- c(strsplit(line, ' ')[[1]])
			setTxtProgressBar(pb, i)

			
			}
		}

		close(open_wig)
		close(pb)
		message('Imported methylation data for ', length(meth_list), ' transcripts')
		return(meth_list)
}

#add missing bases, set to 0 methylation. applies to object: transcrript within experiment
fill_missing <- function(meth_list_element) {
  all <- seq.int(1:max(meth_list_element$Base))
  w_missing <- meth_list_element$Fraction[match(all,meth_list_element$Base)]
  w_missing[is.na(w_missing)] <- 0
  missing_added <- as.data.frame(cbind(all, w_missing))
  colnames(missing_added) <- colnames(meth_list_element)
  
  return(missing_added)
 }

#set as 0,1,2 for degree of methyation at each base. applies to object: transcript within experiment
meth_encode <- function(meth_list_element, low_for_0 = 0.5, high_for_2 = 0.9) {
  
  coded_meth <- cut(meth_list_element$Fraction, c(0, low_for_0, high_for_2, 1.1), labels = c(0,1,2), right=FALSE)
  return(coded_meth)
}

#add experiment names to list. applies to object: all meth data in a list
name_by_experiment <- function(metameth_list, filenames) {
  name_split <- strsplit(filenames, '\\.')
  exps <- sapply(name_split, function(x) x[1])
  names(metameth_list) <- exps
 return(metameth_list)
} 

#add in missing bases and codify. applies to object: experiment
adjust_experiment <- function(meth_list) {
  meth_list  <- lapply(meth_list, fill_missing)
  coded_meth <- lapply(meth_list, meth_encode)
  for (i in 1:length(meth_list)) {
    meth_list[[i]]$coded <- coded_meth[[i]]
  }
  return(meth_list)
}

get_coded_columns <- function(metameth_list, idxs, transcript) {
  transcript_dfs <- lapply(metameth_list[idxs], `[`, transcript)
  codes <- sapply(unlist(transcript_dfs, recursive = FALSE), `[[`, 'coded')
  lens <- sapply(codes, length)
  #check if lengths are different between experiments
  if (!(all(lens == mean(lens)))) codes <- equalize_code_lengths(codes)
  return(as.data.frame(codes))
  }

equalize_code_lengths <- function(codes) {
  m <- max(sapply(codes, length))
  for (i in 1:length(codes)) {
    if (length(codes[[i]]) < m) {
      codes[[i]] <- as.character(c(as.vector(codes[[i]]), rep(0, m-length(codes[[i]]))))
    }
  }
  return(codes)
}

#merge coded methylation data. applies to object: all meth data in list
matrixify <- function(metameth_list) {
  meth_matrix <- list()
  tn <- lapply(metameth_list, names)
  tn_u <- unique(unlist(tn))
  tn_link <- lapply(tn_u, grep, tn)
  names(tn_link) <- tn_u
  for (t in tn_u) {
    to_take <- tn_link[[t]]
    codes <- get_coded_columns(metameth_list, to_take, t)
    names(codes) <- names(tn)[to_take]
    meth_matrix[[t]] <- codes

  }
  return(meth_matrix)
}

#barcode plot for a specific transcript. applies to object: experiment + transcript name
plot_meth <- function(meth_df, transcript_name) {
	meth_50 <- subset(meth_df[[transcript_name]], Fraction > 0.50)
	meth_90 <- subset(meth_df[[transcript_name]], Fraction > 0.90)
	ggplot() + geom_vline(xintercept = meth_50$Base, alpha = meth_50$Fraction, color = '#94bdc9') + 
	geom_vline(xintercept = meth_90$Base, color = '#2e66a9')+ labs(x = 'Genome Position (bp)', subtitle = transcript_name) +
	theme_classic()
	}

#get differential methyation sites from coded_meth_matrix. applies to object: transcript target from meth_matrix
return_diff <- function(meth_matrix_element, ignore = NULL) {
  if (!is.null(ignore)) {
    meth_matrix_element <- meth_matrix_element[,-ignore]
  }
  diff_idx <- apply(meth_matrix_element, MARGIN = 1, FUN=check_diff_row)
  return(meth_matrix_element[diff_idx,])
  
  }
  


#check if methylation is different at base for each experiment. applies to object: row vector of coded methylation data.
check_diff_row <- function(row) {
  u <- unique(row)
  return(length(u) != 1)
  
}


methpipe <-function(wig_files) {
  meth <- lapply(wig_files, read_wig)
  meth <- name_by_experiment(meth, wig_files) 
  meth <- lapply(meth, adjust_experiment)
  meth_matrix <- matrixify(meth)
  
  return(meth_matrix)
}

#Not useful
# #return n transcripts with highest mean fraction modified score. Not sure mean is most useful. Maybe change it to greatest # of sites > 0.5?
# get_most_meth <- function(meth_df, n =10) {
# 	mean_fracs <- aggregate(meth_df[,'Fraction'], FUN=mean, by = list(meth_df$Transcript))
# 	colnames(mean_fracs) <- c('Transcript', 'Mean_Fraction')
# 	mean_fracs <- mean_fracs[order(mean_fracs$Mean_Fraction, decreasing = T), c(1,2)]
# 	return(head(mean_fracs, n =n))
# }
