source('scDataAnalysis_Utilities_simp.R')
library(igraph)
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(readxl)

library(Rcpp)
sourceCpp(code = '
          #include <Rcpp.h>
          using namespace Rcpp;
          // get overlap between two data frame, output a vector, in which 
          // the ith component to be 1 if the ith row of data1 has overlap in data2, zero otherwise
          // [[Rcpp::export]]
          IntegerVector getOverlaps_C1(DataFrame dat1, DataFrame dat2) {
          CharacterVector chr1 = dat1["chr"];
          CharacterVector chr2 = dat2["chr"];
          NumericVector start1 = dat1["start"];
          NumericVector end1 = dat1["end"];
          NumericVector start2 = dat2["start"];
          NumericVector end2 = dat2["end"];
          
          int n1 = chr1.size(), n2 = chr2.size();
          NumericVector midP1(n1), len1(n1), len2(n2), midP2(n2);
          IntegerVector over1(n1);
          
          len1 = (end1 - start1)/2;
          midP1 = (end1 + start1)/2;
          
          len2 = (end2 - start2)/2;
          midP2 = (end2 + start2)/2;
          
          for(int i=0; i<n1; i++){
          over1[i] = 0;
          for(int j=0; j<n2; j++){
          if((chr2[j] == chr1[i]) && (fabs(midP1[i] - midP2[j]) <= max(NumericVector::create(len1[i], len2[j])))){
          over1[i] = 1;
          break;
          }
          }
          }
          
          return(over1);
          }
')


## virids style, downsample, match_cell, colomn color by sample/cluster
plot_enrich_tf1 <- function(sele.tfs, zscore.all, bc_clusters,
                            cluster.levels = NULL, 
                            cluster.col = NULL,
                            up.qt = 0.95, low.qt = 0.05,
                            ndownsample = 2000, match_cell = F,
                            reScale = F, cluster.rows = F,
                            color_style = 'virid'){
  sele.zscores = zscore.all[sele.tfs, ]
  
  if(reScale) sele.zscores = t(scale(t(sele.zscores), center = T, scale = T))
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample and match_cell
  ncell.cl = min(table(bc_clusters$cluster))
  set.seed(2020)
  bc_clusters.down = bc_clusters
  if(match_cell){
    bc_clusters.down = NULL
    for(cl0 in unique(bc_clusters$cluster)){
      tmp = bc_clusters[cluster == cl0]
      if(nrow(tmp) > ncell.cl) tmp = tmp[sort(sample((1:nrow(tmp)), ncell.cl)), ]
      bc_clusters.down = rbind(bc_clusters.down, tmp)
    }
  }
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = bc_clusters$cluster,
                          'barcode' = bc_clusters$barcode,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  nc = length(unique(bc_clusters$cluster))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(is.null(cluster.col)){
    if(nc >= 3) color_cluster = getPalette(nc)
    if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
    names(color_cluster) = sort(unique(bc_clusters$cluster))
  }else{
    color_cluster = cluster.col
  }
  
  
  ## order 
  ann_column = ann_column[order(factor(ann_column$cluster,
                                       levels = sort(unique(ann_column$cluster)))), ]
  if(is.null(cluster.levels)) {
    
    ann_column = ann_column[order(factor(ann_column$cluster,
                                         levels = sort(unique(ann_column$cluster)))), ]
    
  }else{
    ann_column = ann_column[order(factor(ann_column$cluster,
                                         levels = cluster.levels)), ]
    
    color_cluster = color_cluster[cluster.levels]
  }
  
  
  ann_colors = list('cluster' = color_cluster)
  
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode <- NULL
  
  color_fun = viridis(100)
  if(color_style == 'purple-yellow') color_fun = PurpleAndYellow()
  ph <- pheatmap::pheatmap(sele.zscores, cluster_cols = F, 
                           cluster_rows = cluster.rows, 
                           show_colnames = F, fontsize = 10,
                           annotation_col = ann_column, 
                           color = color_fun,
                           annotation_colors = ann_colors, 
                           fontsize_row = 12)
  return(ph)
  
  
  
  
}


## virids style, downsample, match_cell, column color by sample & cluster
plot_enrich_tf2 <- function(sele.tfs, zscore.all, bc_clusters,
                            cluster.levels = NULL, sample.levels = NULL,
                            up.qt = 0.95, low.qt = 0.05,
                            ndownsample = 2000, 
                            match_cell = F, reScale = F, 
                            cluster.rows = F,
                            color_style = 'virid',
                            order.within = 'sample'){
  sele.zscores = zscore.all[sele.tfs, ]
  
  if(reScale) sele.zscores = t(scale(t(sele.zscores), center = T, scale = T))
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample and match_cell
  ncell.cl = min(table(bc_clusters$sample))
  set.seed(2020)
  bc_clusters.down = bc_clusters
  if(match_cell){
    bc_clusters.down = NULL
    for(sample0 in unique(bc_clusters$sample)){
      tmp = bc_clusters[sample == sample0]
      if(nrow(tmp) > ncell.cl) tmp = tmp[sort(sample((1:nrow(tmp)), ncell.cl)), ]
      bc_clusters.down = rbind(bc_clusters.down, tmp)
    }
  }
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = bc_clusters$cluster,
                          'barcode' = bc_clusters$barcode,
                          'sample' = bc_clusters$sample,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  nc = length(unique(bc_clusters$cluster))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(nc >= 3) color_cluster = getPalette(nc)
  if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
  names(color_cluster) = sort(unique(bc_clusters$cluster))
  
  nsample = length(unique(bc_clusters$sample))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(nsample <= 3) color_sample = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nsample]
  if(nsample > 3) color_sample = getPalette(nsample)
  names(color_sample) = sample.levels
  
  ## order sample by age
  ann_column = ann_column[order(factor(ann_column$sample,
                                       levels = sample.levels[sample.levels %in% ann_column$sample])), ]
  ann_column = ann_column[order(factor(ann_column$cluster,
                                       levels = sort(unique(ann_column$cluster)))), ]
  
  if(order.within == 'sample'){
    if(!is.null(cluster.levels)) {
      ann_column = ann_column[order(factor(ann_column$cluster,
                                           levels = cluster.levels)), ]
      
      color_cluster = color_cluster[cluster.levels]
    }
    
    if(!is.null(sample.levels)){
      ann_column = ann_column[order(factor(ann_column$sample,
                                           levels = sample.levels)), ]
      
      color_sample = color_sample[sample.levels]
    }
    
  }else{
    if(!is.null(sample.levels)){
      ann_column = ann_column[order(factor(ann_column$sample,
                                           levels = sample.levels)), ]
      
      color_sample = color_sample[sample.levels]
    }
    
    if(!is.null(cluster.levels)) {
      ann_column = ann_column[order(factor(ann_column$cluster,
                                           levels = cluster.levels)), ]
      
      color_cluster = color_cluster[cluster.levels]
    }
    
    
    
  }
  
  ann_colors = list('cluster' = color_cluster,
                    'sample' = color_sample)
  
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode <- NULL
  
  color_fun = viridis(100)
  if(color_style == 'purple-yellow') color_fun = PurpleAndYellow()
  
  ph <- pheatmap::pheatmap(sele.zscores, cluster_cols = F, 
                           cluster_rows = cluster.rows, 
                           show_colnames = F, fontsize = 10,
                           annotation_col = ann_column, 
                           color = color_fun,
                           annotation_colors = ann_colors, 
                           fontsize_row = 12)
  return(ph)
  
  
  
  
}
