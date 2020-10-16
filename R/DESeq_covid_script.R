<<<<<<< HEAD
#Developed by Josie Gleeson 2020, 
#based on deseq2-analysis-template.R script by Stephen Turner: 
#https://gist.github.com/stephenturner/f60c1934405c127f09a6
=======
#Developed by Josie Gleeson 2020, based on the 'deseq2-analysis-template.R' script by Stephen Turner https://gist.github.com/stephenturner/f60c1934405c127f09a6
>>>>>>> 2b938b2a09458e98be5f7064e604b46fcf9775c2

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  counts <- args[1]
  output <- args[2]
  time1 <- args[3]
  time2 <- args[4]
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
  
  # Import data 
  countdata <- read.table(gzfile(counts), header=TRUE )
  countdata <- as.matrix(countdata)
  
  # Assign conditions
  (condition <- factor(c(rep(time1, 3), rep(time2, 3))))
  
  # Make DESeq dataset
  (coldata <- data.frame(row.names=colnames(countdata), condition))
  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
  
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
  write.csv(resdata, paste0("diff_expr_", output, ".csv"))
  
  # RLD for viewing
  rld <- DESeq2::rlogTransformation(dds)
  
  #Set colours for plotting
  mycols <- brewer.pal(8, "Accent")[1:length(unique(condition))]
  
  # PCA
  rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="Principal Component Analysis", textcx=1, ...) {
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
  pdf(paste0("pca_1_2_", output, ".pdf"), 18, 18, pointsize=50)
  rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-1, 1), ylim=c(-1, 1))
  dev.off()
  
  
  # Volcano Plot
  volcanoplot <- function (res, lfcthresh=1, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
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
  pdf(paste0("volcanoplot_", output, ".pdf"), 18, 18, pointsize=20)
  volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=0.8, xlim=c(-10, 10), legendpos="topright")
  dev.off()
  
  
}

main()

