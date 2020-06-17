#Author: Daniel Rawlinson

#add this into R with 'source(read_tombo_wig.R)' then run the function with 'read_wig(your file)'
#TO plot for a specific transcript, run 'plot_meth(name of data, your transcript of interest)'


library(ggplot2)

read_wig <- function(wig_file) {
	wig_data <- list()
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
			message('Finished!')
			break
			}
		else if ( startsWith(line, 'variableStep') ) {
			transcript <- gsub('chrom=', '',strsplit(line, ' ')[[1]][2])
			getTxtProgressBar(pb)
			}
		else if (startsWith(line, 'track')) {
			next
			}
		else { 
			i=i+1
			#base_frac <- strsplit(line, ' ')[[1]] added in-line below
			wig_data[[i]] <- c(strsplit(line, ' ')[[1]], transcript)
			setTxtProgressBar(pb, i)
			#print(line)
			#print(line)
			#print(line)
			
			}
		}
		meth_df <- data.frame(t(sapply(wig_data,c)))
		colnames(meth_df) <- c('Base','Fraction','Transcript')
		meth_df$Base <- as.integer(as.character(meth_df$Base))
		meth_df$Fraction <- as.numeric(as.character(meth_df$Fraction))
		close(open_wig)
		close(pb)
		message('Imported methylation data for ', length(levels(meth_df$Transcript)), ' transcripts')
		return(meth_df)
		}

plot_meth <- function(meth_df, transcript_name) {
	meth_50 <- subset(meth_df[which(meth_df$Transcript == transcript_name),], Fraction > 0.50)
	meth_90 <- subset(meth_50[which(meth_df$Transcript == transcript_name),], Fraction > 0.90)
	ggplot() + geom_vline(xintercept = meth_50$Base, alpha = meth_50$Fraction, color = '#94bdc9') + 
	geom_vline(xintercept = meth_90$Base, color = '#2e66a9')+ labs(x = 'Genome Position (bp)', subtitle = transcript_name) +
	theme_classic()
	}

#return n transcripts with highest mean fraction modified score. Not sure mean is most useful. Maybe change it to greatest # of sites > 0.5?
get_most_meth <- function(meth_df, n =10) {
	mean_fracs <- aggregate(meth_df[,'Fraction'], FUN=mean, by = list(meth_df$Transcript))
	colnames(mean_fracs) <- c('Transcript', 'Mean_Fraction')
	mean_fracs <- mean_fracs[order(mean_fracs$Mean_Fraction, decreasing = T), c(1,2)]
	return(head(mean_fracs, n =n))
}
