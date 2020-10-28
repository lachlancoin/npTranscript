# Developed by Josie Gleeson (GitHub: josiegleeson) 2020, based on the 'deseq2-analysis-template.R' script by Stephen Turner https://gist.github.com/stephenturner/f60c1934405c127f09a6
# Execute as: Rscript DESeq_covid_script.R [count_file.csv] [output_prefix] [time1] time2]

suppressPackageStartupMessages({
  library(DESeq2)
  library(RColorBrewer)
  library(gplots)
  library("pheatmap")
  library("pcaExplorer")
  library("ggplot2")
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggrepel)
    })

loadDE <- function() {
  
  # Import data 
  message(paste('importing DE Data'))
  count_files <- list.files(path = 'DE', recursive = T, full.names = T)
  cell_types <- lapply(strsplit(x = count_files, split = '/'), rev) %>% sapply('[', 2) 
  print(cell_types)
  
  print(count_files)
  
  countdata <- lapply(count_files, FUN = 
                        function(x) {read.table((gzfile(x)), header = T, row.names = 1) %>%
                            subset(select = -c(1:5)) %>% 
                            as.matrix() } )
  names(countdata) <- cell_types
  return(countdata)
}

#PCA Function
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main=paste('PCA', cell1, time1, 'vs', cell2, time2), textcx=1, ...) {
  print(main)
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1: ",as.character(pc1var),"% variance")
  pc2lab <- paste0("PC2: ",as.character(pc2var),"% variance")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}


.volcano<-function(df, logFCthresh = 0.5, top=20, prefix="", useadj=TRUE,exclude=NULL, main){
  if(dim(df)[[1]]==0) return(NULL)
  if(!is.null(exclude)){
    df = df[!(df$grp %in% exclude),,drop=F]
  }
  if(useadj){
    pthresh = sort(df$padj)[top]
    ggp<-ggplot(df, aes(x =log2FoldChange, y = -log10(padj),color = ifelse(abs(log2FoldChange)>0.6,"red","grey")))  + theme_bw()
    ggp<-ggp+  geom_point() +  xlab(expression("Fold Change, Log"[2]*"")) +  ylab(expression("Adjusted P value, Log"[10]*"")) 
  }else{
    pthresh = sort(df$pvalue)[top]
    ggp<-ggplot(df, aes(x =log2FoldChange, y = -log10(pvalue),color = ifelse(abs(log2FoldChange)>0.6,"red","grey"))) + theme_bw()
    ggp<-ggp+  geom_point() +  xlab(expression("Fold Change, Log"[2]*"")) +  ylab(expression(" P value, Log"[10]*"")) 
  }
  print(paste("pthresh",pthresh))
  ggp<-ggp+  geom_vline(
    xintercept = c(-0.6,0.6),
    col = "red",
    linetype = "dotted",
    size = 1) 
  ggp<-ggp+  geom_hline(
    yintercept = c(-log10(0.01),-log10(0.05)),
    col = "red",
    linetype = "dotted",
    size = 1)
  ggp<-ggp+  theme(
    plot.title = element_text(size=20)
  )
  ggp<-ggp+  theme(legend.position = "none")+
    scale_colour_manual(values = c("grey", "red")) 
  if(useadj){
    ggp<-ggp+  geom_text_repel(data=subset(df,abs(log2FoldChange) >= logFCthresh & padj < pthresh),
                               aes(log2FoldChange, -log10(padj), label = Row.names),size = 3, color="steelblue")
  }else{
    ggp<-ggp+  geom_text_repel(data=subset(df,abs(log2FoldChange) >= logFCthresh & pvalue < pthresh),
                               aes(log2FoldChange, -log10(pvalue), label = Row.names),size = 3, color="steelblue")
  }
  ggp<-ggp+ggtitle(main)
  ggp
}



# Volcano Plot
volcanoplot <- function (res, lfcthresh=0.5, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=T, textcx=0.8, remove_spurious = F, ...) {
  print(remove_spurious)
  if (remove_spurious == T) {
    res <- subset(res, spurious == F)
  }
  
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Row.names, cex=textcx, ...))
  }
  #legend(legendpos, xjust=1, yjust=1, legend=c(paste("p-adj<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}


runDE <- function(count_list, cell1, cell2, time1, time2, thresh=50, plot_params= NULL) {
  
  #Trim to relevant columns
  print(paste('DE', cell1, cell2, time1, time2))
  print(colnames(count_list[[cell1]]))
  print(colnames(count_list[[cell2]]))
  first_cond_idx <- which(grepl(time1, x= colnames(count_list[[cell1]])))
  second_cond_idx <- which(grepl(time2, x= colnames(count_list[[cell2]])))
  
  count_trim <- merge( x = count_list[[cell1]][,first_cond_idx], 
                       y = count_list[[cell2]][, second_cond_idx], 
                       by = "row.names", sort = FALSE) %>% 
                transform(row.names=Row.names, Row.names=NULL) %>%
              as.matrix()

if(!is.null(plot_params))  count_trim=.subsetFCFile(count_trim,plot_params)
  count_trim=count_trim[apply(count_trim,1,sum)>thresh*dim(count_trim)[2],,drop=F] 
 # print(head(count_trim))
  #uncomment below if all = T in merge
  #count_trim[is.na(count_trim)] <- 0

  # Assign conditions
  (condition <- factor(c(rep(paste('1', time1, cell1, sep = '_'), length(first_cond_idx)), rep(paste('2',time2, cell2, sep = '_'), length(second_cond_idx)))))
  
  # Make DESeq dataset
  (coldata <- data.frame(row.names=colnames(count_trim), condition))
  dds <- DESeqDataSetFromMatrix(countData=count_trim, colData=coldata, design=~condition)
  
  keep <- rowSums(counts(dds)) >= 5
  dds <- dds[keep,]
  
  
  
  
  # Run DESeq2 pipeline and save the results
  dds <- tryCatch( {
    DESeq(dds, fitType = 'parametric') },
    error = function(cond) {
      message(paste('Parametric fit failure :', cond))
      message(paste('attempting DE with fitType = mean'))
      DESeq(dds, fitType = 'mean') }
    
  )
                  
  res <- DESeq2::results(dds)
  
##For flagging of spurious results
count_results_melt <- counts(dds, normalized=TRUE) %>%
    `colnames<-`(condition) %>%
    reshape2::melt(measure.vars = cols)  %>%
    group_by(Var1, Var2)
#find num of zero counts per group
num_zeros <- mutate(count_results_melt, is_zero = ifelse(value ==0, 1, 0)) %>% 
      mutate(value = NULL) %>% 
      dcast(Var1 ~ Var2, fun.aggregate = sum, value.var = 'is_zero')
#get difference in mean counts
count_results_melt <-   dcast(count_results_melt, Var1 ~ Var2, fun.aggregate = mean, ) %>% 
      mutate(diff_mean =  .[[3]]- .[[2]])
#merge flag data
count_results_with_mean <- merge(x = count_results_melt, y = num_zeros, by = 'Var1', suffixes = c('.mean', '.zero_counts')) %>% 
  transform(row.names = Var1, Var1=NULL) %>%
  merge(x = as.data.frame(counts(dds, normalized=TRUE)), y = ., by = 'row.names') %>%
  transform(row.names = Row.names, Row.names = NULL)
  


  
  ## Order by adjusted p-value
  res <- res[order(res$padj), ]
  ## Merge with normalized count data
  resdata <- merge(as.data.frame(res), count_results_with_mean, by="row.names", sort=FALSE) %>% 
    mutate(rev_direction =  ifelse(diff_mean * log2FoldChange < 0, T, F)  , 
           spurious = ifelse(rev_direction == T & .[[length(.)-1]] + .[[length(.)-2]] >1 , T, F ))
  
  
  ## Write results
  #write.csv(resdata, paste0("diff_expr_", output, ".csv"))
  
  # RLD for viewing
  rld <- DESeq2::rlogTransformation(dds)
  
  #Set colours for plotting
  mycols <- brewer.pal(8, "Accent")[1:length(unique(condition))]
  
  # PCA
  #pdf(paste0("pca_1_2_", output, ".pdf"), 18, 18, pointsize=50)
  rpca_command <- list(rld = rld, colors=mycols, intgroup="condition", xlim=c(-10, 10), ylim=c(-10, 10), main = paste('PCA', cell1, time1, 'vs', cell2, time2))
  #dev.off()
  
  #Volcano
  #pdf(paste0("volcanoplot_", output, ".pdf"), 18, 18, pointsize=20)
  vp_command <- list(res = resdata, lfcthresh=0.5, sigthresh=0.05, textcx=0.8, xlim=c(-10, 10), legendpos="topright", main = paste('DESeq2', cell1, time1, 'vs', cell2, time2))
  vp_ggp_command <- list(df = resdata, main = paste('DESeq2', cell1, time1, 'vs', cell2, time2))
  #dev.off()
  
  return(list(data = resdata, rld_pca_params = rpca_command, volcano_params = vp_command, volcano_ggp_params = vp_ggp_command))
}

#DE_out <- main()


