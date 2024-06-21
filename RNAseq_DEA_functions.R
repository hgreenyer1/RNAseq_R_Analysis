#Helper functions for conducting Differential Expression analysis from salmon quant outputs 



#libraries ------------
library(RVenn)
library(openxlsx)

#Functions---------------

#' construct annotation df from input file names 
#'
#' @param s_dir path to directory containing salmon_quant files 
#' @param out_dir string - option to specify output directory
#' @param out_pref prefix of outputs annotation file to export 
#'
#' @export annotation option to export annotation dataframe  
#'
#' @examples make_annotation('/path/to/salmon/files/', out_dir='NULL'/home/out_figures', out_pref='MCF10a_degron')
get_anno <- function(s_dir, out_dir=NULL, out_pref=NULL){
  
  #note that this function assumes a very specific file name format:
  #cell_conditon_timepoint_replicate_etc..
  #example: MCF10A_DTAG_3H_rep1_otherStuff.....salmon_quant
  
  #parse file names 
  sample <- list.files(s_dir, pattern = ".salmon_quant") %>% str_replace_all( ".salmon_quant", "")
  cell = sapply(strsplit(sample, '_'), `[`, 1)
  condition = sapply(strsplit(sample, '_'), `[`, 2)
  timepoint = sapply(strsplit(sample,'_'), `[`,3)
  replicate = sapply(strsplit(sample, '_'), `[`, 4)
  directory = list.files(s_dir, pattern = ".salmon_quant")
  
  annotation <- data.frame(sample = sample, directory= directory, cell = cell, condition = condition, timepoint = timepoint, replicate = replicate, stringsAsFactors = TRUE)
  
  #option to export file annotation 
  if(!(is.null(out_dir))){
    write.csv(annotation, paste(out_dir,out_pref,'_anno.csv', sep = ''))
  }
  
  return(annotation)
  
}

#' get txi.salmon object from input directories and annotations 
#'
#' @param ref_dir path to annotation gtf
#' @param s_dir path to directory containing salmon_quant files 
#' @param annotation string - option to specify output directory
#'
#' @examples make_annotation(ref_dir, s_dir)
get_txi <- function(ref_dir, s_dir, annotation){
  
  #function assumes some annotations 
  
  txdb <- makeTxDbFromGFF(ref_dir)
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
  
  files <- file.path(s_dir,annotation$directory, "quant.sf")
  names(files) <- paste0(annotation$sample)
  txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
  
  return(txi.salmon)
  
}

#' get DESeq object from sample annotation and txi.salmon object 
#'
#' @param annotation sample annotations 
#' @param txi.salmon salmon reference object 
#' @param frmt format from which DESeq object will be constructed
#' @param export_norm_counts option to export normalized count table from dds object
#' @param out_dir output directory to export norm counts to 
#' @param out_pref prefix for output file name 
#' 
#' @return returns dds object
#'
#' @examples get_dds(MCF10A_anno)
get_dds <- function(annotation, quant_data,  frmt = 'txi', dsgn = ""){
  
  #define column data - assumes sample annotation
  coldata = annotation
  rownames(coldata) = annotation$sample
  
  
  if(frmt == 'txi'){
  
    #get DESeq object from txi 
    dds <- DESeqDataSetFromTximport(txi = quant_data, 
                                     colData = coldata, 
                                    design = ~replicate + condition)
  }
  else if(frmt == 'mat'){
    dds <- DESeqDataSetFromMatrix(countData  = quant_data, 
                                    colData = coldata, 
                                    design = ~replicate + condition)
  }
  
  dds <- DESeq(dds)
  
  return(dds)
  
}

#' construct annotation df from input file names 
#'
#' @param dds path to directory containing salmon_quant files 
#' @param norm option to normalize counts
#' @param out_dir output directory
#' @param out_pref prefix of output file
#'
#' @export annotation option to export annotation dataframe  
#'
#' @examples 
export_counts <- function(dds, norm = TRUE, out_dir=NULL, out_pref=NULL){
  
  #extract counts from dds - option to normalize
  counts <- as.data.frame(counts(dds, normalized=norm))
    
  #remove row names and add as col
  counts_out <- cbind(rownames(counts),counts)
  rownames(counts_out) <- NULL
  names(counts_out)[1] <- 'ID'
  
  #version IDs
  counts_out <- counts_out %>%  
    mutate(ID = str_replace(ID, ".[0-9]+$", ""))
    
  #get gene symbols 
  g_symbol <- gconvert(counts_out$ID, organism = "hsapiens", target = "ENSG", numeric_ns = "", mthreshold = Inf, filter_na = F)
  counts_out <- cbind(g_symbol, counts_out)
  names(counts_out)[1] <- 'Gene_Symbol'
  
  #write to csv 
  if(norm == TRUE){
    write_csv(counts_out, paste(out_dir,out_pref,'_NormCounts.csv', sep = ''))
  }  
  else{
    write_csv(counts_out, paste(out_dir,out_pref,'_RawCounts.csv', sep = ''))
  }


}


#' get full DESeq object from annotation and txi.salmon object 
#'
#' @param dds DESeq object
#'
#' @examples get_dds(annotation, txi.salmon)
get_vst <- function(dds){
  
  #this may not be a necessary function...
  vst <- varianceStabilizingTransformation(dds, blind=TRUE)
  vst_mat <- assay(vst)
  return(vst_mat)
}


#' outputs a bunch of bulk RNAseq QC figures 
#'
#' @param annotation sample annotations 
#' @param txi.salmon salmon reference object 
#' @param dds DESeq object
#' @param vst_mat variance stabilized dds matrix
#'
#' @examples get_dds(annotation, txi.salmon)
RNAseq_QC_figs <- function(annotation, txi.salmon, dds, vst_mat){
  
  #barplot of library sizes
  librarySizes <- colSums(txi.salmon$counts)
  barplot(librarySizes, 
          names=annotation$name, 
          las=2, 
          main="Barplot of library sizes")
  
  plot_total_counts(dds)
  
  logcounts <- log2(vst_mat + 1)
  # make a color vector
  statusCol <- as.numeric(factor(annotation$condition)) + 1
  # Check distributions of samples using boxplots
  boxplot(logcounts, 
          xlab=annotation$name, 
          ylab="Log2(counts)",
          las=2,
          col=statusCol)
  
  plot_library_complexity(dds)
  plot_gene_detection(dds)
  plot_biotypes(dds)
  
  #number of genes that have non–zero counts in all samples
  GeneCounts <- counts(dds)
  idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
  sum(idx.nz)
  
  #estimate the size factors
  DESeq2Table <- estimateSizeFactors(dds)
  sizeFactors(DESeq2Table)
              
  #densities of counts
  multidensity( counts(DESeq2Table, normalized = T)[idx.nz ,],
                xlab="mean counts", xlim=c(0, 1000))
  
  #empircal cumulative distribution functions (ECDFs), essential integrals of the densities and give the probability of observing a certain number of counts equal to x or less given the data.
  multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
             ylab="estimated probability(Fn(x))", xlab="mean counts", xlim=c(0, 1000))
  
} 


#' plot PCA 
#'
#' @param vst vector stabilized dds object 
#' @param color annotation variable mapped to color 
#' @param shape annotation variable mapped to shape
#' @param plot_title string - title for PCA plot 
#'
#' @examples plot_PCA(annotation, txi.salmon)
plot_PCA <- function(vst, color=timepoint, shape=condition, plot_title=NULL){
  
  # Principal components analysis
  pcaData <- plotPCA(vst, intgroup=c('condition','timepoint','replicate'), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar")) 
  
  ggplot(pcaData, aes(x=PC1, y=PC2, color=timepoint, shape=condition)) + 
    geom_point(size =3, aes(fill=timepoint)) + 
    ggtitle(plot_title) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "% variance"))  
  
    #should add a ggsave here 
}

#' extract samples corresponding to a certain annotation 
#'
#' @param txi.salmon large list from all samples
#' @param annotation annotation with all samples 
#' @param sub string of characteristic to subset for
#' 
#' @return returns a DESeq object subset to contain only desired samples 
#'
#' @examples subset_dds(MCF10A_txi.salmon, MCF10A_anno, sub= '3H')
subset_dds <- function(txi.salmon, annotation, sub_col = NULL, sub = NULL){
  
  #get subset of count data as matrix from txi  
  df.salmon_A = as.data.frame(txi.salmon$counts) 
  
  s = annotation$sample[annotation[[sub_col]] == sub]
  
  df.A <- df.salmon_A[,c(s)]
  
  #reformat dataframe 
  df.A <- data.frame(rownames = rownames(df.A), df.A, check.names = FALSE)
  names(df.A)[1] <- "Gene"
  
  salmon.df <- df.A[,c(2:ncol(df.A))] 
  
  mat.A <- as.matrix(as_tibble(salmon.df))
  rownames(mat.A) <- df.A$Gene
  
  #subset annotation
  anno_sub = annotation %>% dplyr::filter(annotation[,sub_col] == sub)
  
  coldata = anno_sub
  rownames(coldata) = anno_sub$sample
  
  #call dds function with subset input 
  dds_sub <- get_dds(annotation = coldata, quant_data = round(mat.A),  frmt = 'mat')
  
  #filter low abundance
  keep <- rowSums(counts(dds_sub) >= 10) >= 3
  
  dds_sub <- dds_sub[keep ,]
  
  
  return(dds_sub)
}


#' extract samples corresponding to a certain annotation 
#'
#' @param dds DESeq object
#'
#' @examples get_DEGs(annotation, txi.salmon)
get_DEGs <- function(dds, ctrst = c('condition','DTAG','DMSO'), export_degs = FALSE, out_dir = NULL, out_pref = NULL, l2fc_cutoff = 1, padj_cutoff = 0.05){

  #get DEGs
  resA <- results(dds, contrast=ctrst)
  
  #make df from results and normalized count data 
  resA <- resA[order(resA$padj), ]
  resA <- as.data.frame(resA)
  resA_data <- data.frame(rownames = rownames(resA), resA, check.names = FALSE)
  
  names(resA_data)[1] <- "Gene"
  
  #remove the version id
  resA_data2 <- resA_data %>%  
    mutate(Gene = str_replace(Gene, ".[0-9]+$", ""))
  
  #get gene symbols 
  resA_symbol <- gconvert(resA_data2$Gene, organism = "hsapiens", target = "ENSG", numeric_ns = "", mthreshold = Inf, filter_na = T)
  names(resA_symbol)[2] <- "Gene"
  
  #merge gene symbols and DESeq results 
  res_data <- merge(resA_symbol, resA_data2, by="Gene", all.x=T, sort=F)
  
  #extract statistically significant DEGs
  deg_df <- res_data %>% dplyr::filter(abs(log2FoldChange) >= l2fc_cutoff & padj <= padj_cutoff)
  
  #option to export degs
  if(export_degs == TRUE){
    
    #export unfiltered DESeq results 
    write_csv(res_data, paste(out_dir,out_pref,'_DESeqRes.csv',sep=''))
    
    #export filtered results 
    tableDEGs <- deg_df[order(deg_df$log2FoldChange), ] %>% dplyr::select("Gene", "name", "log2FoldChange", "padj")
    rownames(tableDEGs) = NULL
    write_csv(tableDEGs, paste(out_dir,out_pref,'_DEGs.csv',sep=''))
  }
  
  #return unfiltered results 
  return(res_data)
  
}


#' extract samples corresponding to a certain annotation 
#'
#' @param dds DESeq object
#'
#' @examples get_DEGs(annotation, txi.salmon)
bulkQC_figs <- function(dds, txi, anno){
  
  #get vst objects 
  vst <- varianceStabilizingTransformation(dds, blind=TRUE)
  vst_mat <- assay(vst)
  
  #barplot of library sizes
  librarySizes <- colSums(txi$counts)
  barplot(librarySizes, 
          names=names(librarySizes), 
          las=2, 
          main="Barplot of library sizes")
  
  plot_total_counts(dds)
  
  logcounts <- log2(vst_mat + 1)
  # make a color vector
  statusCol <- as.numeric(factor(anno$condition)) + 1
  # Check distributions of samples using boxplots
  boxplot(logcounts, 
          xlab="", 
          ylab="Log2(counts)",
          las=2,
          col=statusCol)
  
  plot_library_complexity(dds)
  plot_gene_detection(dds)
  #plot_biotypes(dds)
  
  #number of genes that have non–zero counts in all samples
  GeneCounts <- counts(dds)
  idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
  sum(idx.nz)
  
  #estimate the size factors
  DESeq2Table <- estimateSizeFactors(dds)
  sizeFactors(DESeq2Table)
  #densities of counts
  multidensity(counts(DESeq2Table, normalized = T)[idx.nz ,],
                xlab="mean counts", xlim=c(0, 1000))
  
  #empircal cumulative distribution functions (ECDFs), essential integrals of the densities and give the probability of observing a certain number of counts equal to x or less given the data.
  multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
             ylab="estimated probability(Fn(x))", xlab="mean counts", xlim=c(0, 1000))
}


#' export paired overlap lists from a list of DEGs
#'
#' @param deg_list a named list of DEGs across different condition pairs
#' @param out_pref prefix for the output file name 
#'
#' @examples export_olaps(MCF10a_upreg_DEGs, 'MCF10A_upregulated')
export_olaps <- function(deg_list, out_pref = '', out_dir=''){
  
  #venn diagram object from list of DEGs
  venn_obj <- RVenn::Venn(deg_list)
  
  #get list of overlapping DEGs from venn 
  venn_olap <- overlap_pairs(venn_obj, slice = 'all')
  
  #export xlxs 
  write.xlsx(venn_olap, paste(out_dir,out_pref,'_olaps.xlsx', sep = ''))

}



