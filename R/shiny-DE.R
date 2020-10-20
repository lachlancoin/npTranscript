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
    })

# loadDE <- function() {
# 
#   # Import data
#   message(paste('importing DE Data'))
#   count_files <- list.files(path = 'DE', recursive = T, full.names = T)
#   cell_types <- lapply(strsplit(x = count_files, split = '/'), rev) %>% sapply('[', 2) 
#   print(count_files)
#   countdata <- lapply(count_files, FUN =
#                         function(x) {read.table((gzfile(x)), header = T, row.names = 1) %>%
#                               subset(select = -c(1:5)) %>%
#                                 as.matrix() } )
#   names(countdata) <- cell_types
#   return(countdata)
# }

#PCA Function
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main=paste('PCA', cell, time1, 'vs', cell, time2), textcx=1, ...) {
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

# Volcano Plot
volcanoplot <- function (res, lfcthresh=0.5, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("p-adj<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}


runDE <- function(count_list, cell, time1, time2) {
  #Trim to relevant columns
  countdata = count_list[[cell]]
  first_cond_idx <- which(grepl(time1, x= colnames(countdata)))
  second_cond_idx <- which(grepl(time2, x= colnames(countdata)))
  count_trim <- countdata[,c(first_cond_idx, second_cond_idx)]
  
  # Assign conditions
  (condition <- factor(c(rep(time1, length(first_cond_idx)), rep(time2, length(second_cond_idx)))))
  
  # Make DESeq dataset
  (coldata <- data.frame(row.names=colnames(count_trim), condition))
  dds <- DESeqDataSetFromMatrix(countData=count_trim, colData=coldata, design=~condition)
  
  keep <- rowSums(counts(dds)) >= 5
  dds <- dds[keep,]
  
  # Run DESeq2 pipeline and save the results
  dds <- DESeq(dds)
  res <- DESeq2::results(dds)
  
  ## Order by adjusted p-value
  res <- res[order(res$padj), ]
  ## Merge with normalized count data
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  ## Write results
  #write.csv(resdata, paste0("diff_expr_", output, ".csv"))
  
  # RLD for viewing
  rld <- DESeq2::rlogTransformation(dds)
  
  #Set colours for plotting
  mycols <- brewer.pal(8, "Accent")[1:length(unique(condition))]
  
  # PCA
  #pdf(paste0("pca_1_2_", output, ".pdf"), 18, 18, pointsize=50)
  rpca <- rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-1, 1), ylim=c(-1, 1))
  #dev.off()

  #pdf(paste0("volcanoplot_", output, ".pdf"), 18, 18, pointsize=20)
  vp <- volcanoplot(resdata, lfcthresh=0.5, sigthresh=0.05, textcx=0.8, xlim=c(-10, 10), legendpos="topright", main = paste('DESeq2', cell, time1, 'vs', cell, time2))
  
  #dev.off()
  
  return(list(data = resdata, rld_pca = rpca, volcano = vp))
}

#DE_out <- main()


