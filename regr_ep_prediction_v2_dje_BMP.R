source('scDataAnalysis_Utilities_simp.R')
source('regr_ep_prediction_Utilities.R')
`%notin%` = Negate(`%in%`)

library(tidyr)
library(data.table)

#options(error = recover)

## get a binary matrix indicates the gene-peak affinity
## gene.list are data.table including column name gene_name 
## gene_ann should include gene_name,chr,start,end



# combine all coembed results from a single patient ####
seurat.atac <- readRDS('/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds')


final_matching = smooth.rna = smooth.atac = list()
sampleNames = unique(seurat.atac$sample.name)

## regression ####
print("loading combined")
final_matching = readRDS('Coembed_Results2/combined_final_matching.rds')
smooth.rna = readRDS('Coembed_Results2/combined_scRNA_smoothed.rds')
smooth.atac = readRDS('Coembed_Results2/combined_scATAC_smoothed.rds')


print("adding regression analysis")

## filter peaks that accessible in less than 5% of all cell type
if(T){
  
  peaks.mean.ctype <- sapply(unique(seurat.atac$ETP), function(x){
    rowMeans(seurat.atac@assays$ATAC@data[, seurat.atac$ETP == x] > 0)
  })
  rmax = apply(peaks.mean.ctype, 1,  max)
  print("rmax")
  print(rmax)
  summary(rmax)
  filtered.peaks = names(which(rmax > 0.05))
  filtered.peaks = lapply(filtered.peaks, function(x) unlist(strsplit(x, ','))[1])
  filtered.peaks = do.call('c', filtered.peaks)
}


cand.peaks = colnames(smooth.atac)
cand.peaks = sapply(cand.peaks, function(x) unlist(strsplit(x, ','))[1])
#write.csv(cand.peaks,'cand.peaks.txt')
names(cand.peaks) = NULL
colnames(smooth.atac) = cand.peaks
all(filtered.peaks %in% cand.peaks)

atac.cnames = smooth.atac$cell_name
smooth.atac = subset(smooth.atac, select = filtered.peaks)
acces.frac = sapply(filtered.peaks, function(x) mean(smooth.atac[[x]] > 0))


final.peaks = filtered.peaks[acces.frac > 0.01]
smooth.atac = subset(smooth.atac, select = final.peaks)

## focus on selected genes (or degs)
#degs_list = all.de.list <- readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/16_construct-TRN/all.de.list.rds")
#degs = unique(c(degs_list[[1]]$gene, degs_list[[2]]$gene, degs_list[[3]]$gene))
#surface.de = readRDS("/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/16_construct-TRN/de.surface.markers/de.surface.markers")
#degs = unique(c(degs), unique(surface.de[surface.de$group == 'de.genes',]$gene))

degs <- read.csv("degs_BMP.txt")
colnames(degs) <- "gene_name"

## gene-peak-affinity map
## construct gene peak affinity binary matrix
gene_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/gene_ann_hg38.txt')
gene_ann[, 'Tss' := ifelse(strand == '+', start, end)]

names(degs) <- "gene_name"

final.peaks <- as.data.frame(final.peaks)
colnames(final.peaks) <- "peak_name"


#modifying ...
gene2peak.map <- get_gene2peak_map(gene.list = data.table('gene_name' = degs),
                                   peak_names = final.peaks, 
                                   gene_ann = gene_ann,
                                   distal_dist = 5e6)

saveRDS(gene2peak.map, file = 'EP_Prediction2/test_gene2peak.map.rds')

peaks.used = colnames(gene2peak.map)
degs = degs[degs$gene_name %in% rownames(gene2peak.map),]
write.csv(degs, file = 'EP_Prediction2/test_degs_filter.txt')


## do regression gene by gene
smooth.rna = subset(smooth.rna, select = degs)
saveRDS(smooth.rna, file = 'EP_Prediction2/test_subset_smooth.rna.rds')

smooth.atac = subset(smooth.atac, select = peaks.used)
saveRDS(smooth.atac, file = 'EP_Prediction2/test_subset_smooth.atac.rds')

stime = Sys.time()
regr.list = list()
degs <- unique(degs)
for(x in degs){
  exprs <- as.numeric(smooth.rna[[x]])
  names(exprs) <- NULL
  peaks = names(which(gene2peak.map[x, ] == 1))
  covrs = as.matrix(subset(smooth.atac, select = peaks))
  
  rdata <- data.frame(cbind(exprs, covrs))
  res <- coef(summary(lm(exprs ~ ., data = rdata)))
  colnames(res)[4] <- 'P_value'
  regr.list[[x]] <- res
}
etime = Sys.time()
etime-stime
names(regr.list) = degs
saveRDS(regr.list, "EP_Prediction2/regrRes4ep_prediction_5e6_BMP.rds")  #"EP_Prediction2/regrRes4ep_prediction_5e6.rds"



## summarize/filter loops ####
regr.list = readRDS("EP_Prediction2/regrRes4ep_prediction_5e6_BMP.rds") #EP_Prediction2/regrRes4ep_prediction_5e6.rds
regr.sum <- lapply(degs, function(t){
  x = regr.list[[t]]
  x = x[, c(1, 4)]
  x = data.frame(x)
  x = data.table(x, keep.rownames = T)
  x$gene_name = t
  return(x)
})

regr.sum = do.call('rbind', regr.sum)

regr.sum[, 'p_val_adj' := pmin(1, P_value*nrow(regr.sum))]
regr.sum$fdr = p.adjust(regr.sum$P_value, method = 'fdr')

#regr.filtered = regr.sum[fdr < 0.05 & abs(Estimate) > 0.3 & grepl(rn, pattern = '^chr')]
regr.filtered = regr.sum[fdr < 0.05 & abs(Estimate) > 0.10 & grepl(rn, pattern = '^chr')]

regr.filtered$peak_name = sapply(regr.filtered$rn, function(x) gsub('.', '-', x, fixed = T)  )
regr.filtered$rn <- NULL
#saveRDS(regr.filtered, "EP_Prediction2/regr.filtered_fdrthresh_5e6_0.3_v2_nonvis.rds")
saveRDS(regr.filtered, "EP_Prediction2/regr.filtered_fdrthresh_5e6_0.10_0.05_v2_nonvis_BMP.rds") #"EP_Prediction2/regr.filtered_fdrthresh_5e6_0.10_0.05_v2_nonvis.rds"

## to visualize on ucsc genome browser
## (promoter side: closest peak to TSS -- gene level)
gene_ann.deg = gene_ann[gene_name %in% regr.filtered$gene_name, ]

gene_ann.deg[, 'promoter_start' := Tss - 1000]
gene_ann.deg[, 'promoter_end' := Tss + 1000]

setkey(gene_ann.deg, gene_name)
regr.filtered[, 'start' := gene_ann.deg[J(regr.filtered$gene_name)]$promoter_start]
regr.filtered[, 'end' := gene_ann.deg[J(regr.filtered$gene_name)]$promoter_end]
regr.filtered[, 'chr' := gene_ann.deg[J(regr.filtered$gene_name)]$chr]
regr.filtered[, 'promoter_pos' := paste(chr, start, end, sep = '-')]
#saveRDS(regr.filtered, "EP_Prediction2/regr.filtered_fdrthresh_5e6_0.3_v2.rds")
saveRDS(regr.filtered, "EP_Prediction2/regr.filtered_fdrthresh_5e6_0.10_0.05_v2_BMP.rds") #'EP_Prediction2/regr.filtered_fdrthresh_5e6_0.10_0.05_v2.rds

## filter otherend not overlapping with promoters
tss_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/transcript_ann_hg38.txt')
tss_ann = tss_ann[gene_biotype %in% c('protein_coding', 'lincRNA', 'miRNA')]
tss_ann[, 'Tss' := ifelse(strand == '+', start, end)]

# any peak within promoter region
peak.ann = annPeak2Gene(peaks.used, tss_ann, 2000)
setkey(peak.ann, peak_name)
peaks.nprom = peak.ann[nchar(gene_name) == 0]$peak_name
peaks.prom = peak.ann[nchar(gene_name) > 0]$peak_name

regr.filtered.ep = regr.filtered[peak_name %in% peaks.nprom]

## assign nearest peak to promoter
gene_list <- subset(regr.filtered.ep, 
                    select = c(gene_name, chr, start, end)) %>%
  .[!duplicated(.)]

gene_list2peak = geneOverlapPeak(gene_list, peak_names = peaks.prom,
                                 mid_dist = 1000)
gene_list2peak = gene_list2peak[peak_name != 'Not_Found']
setkey(gene_list2peak, gene_name)

regr.filtered.ep = regr.filtered.ep[gene_name %in% gene_list2peak$gene_name]
regr.filtered.ep[, 'promoter_peak' := gene_list2peak[J(regr.filtered.ep$gene_name)]$peak_name]


regr.filtered.ep = subset(regr.filtered.ep, select = c(gene_name, promoter_pos, 
                                                       promoter_peak, peak_name,
                                                       P_value, p_val_adj, fdr, Estimate))
names(regr.filtered.ep)[4] = 'enhancer_peak'


regr.filtered.ep[, 'chr1' := unlist(strsplit(promoter_peak, '-'))[1], by = promoter_peak]
regr.filtered.ep[, 'start1' := as.integer(unlist(strsplit(promoter_peak, '-'))[2]), by = promoter_peak]
regr.filtered.ep[, 'end1' := as.integer(unlist(strsplit(promoter_peak, '-'))[3]), by = promoter_peak]

regr.filtered.ep[, 'chr2' := unlist(strsplit(enhancer_peak, '-'))[1], by = enhancer_peak]
regr.filtered.ep[, 'start2' := as.integer(unlist(strsplit(enhancer_peak, '-'))[2]), by = enhancer_peak]
regr.filtered.ep[, 'end2' := as.integer(unlist(strsplit(enhancer_peak, '-'))[3]), by = enhancer_peak]

regr.filtered.ep[, 'ep_dist' := abs(start1 + end1 - start2 - end2)/2]
regr.filtered.ep = subset(regr.filtered.ep, select = c(gene_name, promoter_pos, 
                                                       promoter_peak, enhancer_peak, ep_dist,
                                                       P_value, p_val_adj, fdr, Estimate))



fwrite(regr.filtered.ep, file = 'EP_Prediction2/regrRes4_EP_overall_5e6_0.10_fdr0.05_BMP.txt',
       sep = '\t')  #EP_Prediction2/regrRes4_EP_overall_5e6_0.10_fdr0.05.txt


