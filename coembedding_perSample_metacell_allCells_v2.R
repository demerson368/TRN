## transfer label from scRNA to scATAC using Seurat 

source('scDataAnalysis_Utilities_simp.R')
print("sourced scDataAnalysis_Utilities_simp.R")

library(gridExtra)
library(patchwork)
`%notin%` = Negate(`%in%`)

args = commandArgs(T)
sampleID0 = args[1]

seuratAtacPath =  '/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds' #'Seurat_Objects/seurat_atac_non_binary_vap20000.rds'
seuratRNAPath = '/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scRNA/42_add_finalized_cMeta_to_X01_samples/t.all.40.rds'
reduction2use = 'pca' # pca or harmony

## load scRNA data and cell type annotation ####
#seurat.rna <- readRDS('/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/scrna/seurat/NB/filtered/integrated.dynamic.filtered.max_mit_pct_10.min_expressed_genes_1000.max_expressed_genes_10000.min_UMI_2000.max_UMI_50000/after_malig_calling/integrated_seurat_object.with_malignancy_labels.rds')
seurat.rna <- readRDS(seuratRNAPath)
seurat.rna <- UpdateSeuratObject(object = seurat.rna)
seurat.rna = subset(seurat.rna, sample.name == sampleID0)

seurat.atac = readRDS(seuratAtacPath)
seurat.atac <- UpdateSeuratObject(object = seurat.atac)
seurat.atac = subset(seurat.atac, sample.name == sampleID0)
seurat.atac <- NormalizeData(seurat.atac)

print("subset performed")

## recreate seurat object using the sample only
mtx = seurat.atac@assays$ATAC@counts
npc = 30
rs = Matrix::rowMeans(mtx > 0)
filtered.pks = names(which(rs < 0.005))
print("filtered peaks")

#print("filtered.pks: ")
#print(filtered.pks)

nvap = 10000
inherited.mdata <- seurat.atac@meta.data
print("collect meta")
#seurat.atac = R.utils::doCall(doBasicSeurat_atac_updated(mtx, npc = npc, 
#                                         norm_by = 'logNormalize',
#                                         top.variable = nvap,
#                                         regressOnPca = TRUE,
#                                         reg.var = 'nFeature_ATAC', 
#                                         excludePks.fromVAP = filtered.pks, 
#                                         meta.data =inherited.mdata))



#seurat.atac = R.utils::doCall(doBasicSeurat_atac_updated(mtx, npc = npc,
#                                         norm_by = 'logNormalize',
#                                         top.variable = nvap,
#                                         regressOnPca = TRUE,
#                                         reg.var = 'nFeature_ATAC',
#                                         meta.data =inherited.mdata))

# Create Seurat object
print("create seurat")
seurat.atac <- CreateSeuratObject(counts = mtx, assay = 'ATAC', project = 'ATAC')
# Normalize data
print("normalize")
seurat.atac <- NormalizeData(seurat.atac, normalization.method = 'LogNormalize')
print("Identify variable features")
seurat.atac <- FindVariableFeatures(seurat.atac, selection.method = 'vst', nfeatures = nvap)
print("Scale data")
seurat.atac <- ScaleData(seurat.atac, features = VariableFeatures(seurat.atac))
print("Run PCA")
seurat.atac <- RunPCA(seurat.atac, features = VariableFeatures(seurat.atac))
print("Regress out unwanted sources of variation")
reg.var <- 'nFeature_ATAC'
seurat.atac <- ScaleData(seurat.atac, vars.to.regress = reg.var)
#print("Add metadata")
#seurat.atac$meta.data <- inherited.mdata

print("doBasicSeurat_atac_updated called !!!!!!!!!!!!!!!")
seurat.atac = RunUMAP(seurat.atac, reduction = 'pca', dim = 1:npc) %>% 
  FindNeighbors(reduction = 'pca', dim = 1:npc) %>% FindClusters(resolution = 0.2)
print("update seurat.atac")

seurat.atac <- UpdateSeuratObject(seurat.atac)
saveRDS(seurat.atac,file = "intermediate_seurat_ver_1.rds")

p0 = DimPlot(seurat.atac, raster = F) + ggtitle('ATAC')



## use GAS = promote + gene body accessibility for label transfer ####
if(all(names(seurat.atac@assays) != 'ACTIVITY')){
  #activity.matrix = readRDS('Seurat_Objects/mtx_gene_activity_signac.rds')
  #activity.matrix = activity.matrix[, colnames(seurat.atac)]
  #seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
  print("running activity")

  activity.matrix = readRDS('Seurat_Objects/mtx_gene_activity_signac.rds')
  activity.matrix = activity.matrix[, colnames(seurat.atac)]
  #seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

  #atac.mtx = seurat.atac@assays$ATAC@counts #seurat.atac0@assays$ATAC@counts
  #activity.matrix = generate_gene_cisActivity('/mnt/isilon/tan_lab/chenc6/Tools/SingleCellAnnotation/GRCh38/genes/genes.gtf',
  #                                          atac.mtx,
  #                                          include_body = T)

  print("add activity")
  seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

  print("running normalize data")
  seurat.atac <- NormalizeData(
    object = seurat.atac,
    assay = 'ACTIVITY',
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat.atac$nCount_ACTIVITY)
  )
  
}
DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac$tech = 'ATAC'
seurat.rna$tech = 'RNA'

seurat.atac <- FindVariableFeatures(seurat.atac, nfeatures = 5000)
vags = VariableFeatures(seurat.atac)

refAssay = 'RNA'

## restrict RNA to the atac region if atac has 1 region but rna has 2
if(length(unique(seurat.rna$sample_name)) > length(unique(seurat.atac$region))){
  reg.atac = unique(seurat.atac$region)
  seurat.rna = subset(seurat.rna, sample_name == paste('NB_7767', sampleID0, reg.atac, sep = '_'))
}

## transfer label 
DefaultAssay(seurat.rna) = refAssay
#seurat.rna = LogNormalize(seurat.rna)
seurat.rna = FindVariableFeatures(seurat.rna, nfeatures = 2200)
vegs = VariableFeatures(seurat.rna)
gfreq = rowMeans(seurat.rna@assays$RNA@counts > 0)
rgenes = names(which(gfreq < 0.0025))
vegs = setdiff(vegs, rgenes)
VariableFeatures(seurat.rna) <- vegs
print("running scaled data rna")
seurat.rna = ScaleData(seurat.rna) %>% RunPCA(verbose = F, npcs = 30) %>% RunUMAP(reduction = 'pca', dim = 1:30)
seurat.rna = FindNeighbors(seurat.rna, reduction = 'pca', dims = 1:30) %>% FindClusters(resolution = 0.5)
p1 <- DimPlot(seurat.rna, group.by = 'seurat_clusters') + ggtitle('RNA')


genes4anchors = vegs

if(length(unique(seurat.rna$sample_name)) > length(unique(seurat.atac$region))) genes4anchors = intersect(vegs, vags)


transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.atac,
                                        features = genes4anchors,
                                        reference.assay = refAssay,
                                        normalization.method = 'LogNormalize',
                                        query.assay = "ACTIVITY",
                                        reduction = 'cca',
                                        k.anchor = 30, k.filter = 200)


## co-embedding ####
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to

refdata <- GetAssayData(seurat.rna, assay = refAssay, 
                        slot = "data")[genes4anchors %in% rownames(seurat.rna), ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
kweight = min(50, nrow(transfer.anchors@anchors) - 5)
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat.atac[["pca"]],
                           dims = 1:ncol(seurat.atac[["pca"]]), k.weight = kweight)

# this line adds the imputed data matrix to the seurat.atac object
seurat.atac[[refAssay]] <- imputation
DefaultAssay(seurat.atac) = refAssay
#seurat.atac[['ATAC']] <- NULL

print("coembed")
coembed <- merge(x = seurat.rna, y = seurat.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
print("running scaled data anchors")
coembed <- ScaleData(coembed, features = genes4anchors, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes4anchors, verbose = FALSE, npcs = 30)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed <- FindNeighbors(coembed, dims = 1:30, reduction = 'pca')
coembed <- FindClusters(coembed, resl = 0.3)

p2 <- DimPlot(coembed, group.by = 'tech', label = T) + ggtitle('Coembedded')
p3 <- DimPlot(coembed, group.by = 'seurat_clusters', label = T) + ggtitle('Coembedded')


message('Co-embedding done, working on metacell calling ... ')


## call metacell on coembedded obj ####
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

coembed <- RenameCells(coembed,  new.names = paste0(coembed$tech, '_', colnames(coembed)))

coembed <- SetupForWGCNA(
  coembed,
  gene_select = "fraction", # the gene selection approach
  #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "test", # the name of the hdWGCNA experiment
  features = genes4anchors
)

# construct metacells  in each group
coembed <- MetacellsByGroups(
  seurat_obj = coembed,
  group.by = c("seurat_clusters"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 3, # maximum number of shared cells between two metacells
  min_cells = 100,
  reduction = 'pca',
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)


# normalize metacell expression matrix:
coembed <- NormalizeMetacells(coembed)
metacell_obj <- GetMetacellObject(coembed)

if(ncol(metacell_obj) < 15){
  coembed <- MetacellsByGroups(
    seurat_obj = coembed,
    group.by = c("seurat_clusters"), # specify the columns in seurat_obj@meta.data to group by
    k = 20, # nearest-neighbors parameter
    max_shared = 3, # maximum number of shared cells between two metacells
    min_cells = 100, # ignore cluster with less than min_cells
    reduction = 'pca',
    ident.group = 'seurat_clusters'
  )
  # normalize metacell expression matrix:
  coembed <- NormalizeMetacells(coembed)
  metacell_obj <- GetMetacellObject(coembed)
  
}

coembed <- ScaleMetacells(coembed, features=VariableFeatures(seurat.rna))
message(paste('#of metacells:', ncol(metacell_obj)))

if(ncol(metacell_obj) > 50){
  coembed <- RunPCAMetacells(coembed, features=VariableFeatures(seurat.rna), npcs = 10)
  #coembed <- RunHarmonyMetacells(coembed, group.by.vars='Sample')
  coembed <- RunUMAPMetacells(coembed, reduction='pca', dims=1:10)
  
  p3 <- DimPlotMetacells(coembed, group.by = 'seurat_clusters') + ggtitle("cluster_metacell")
  
}
  
ggsave(p0 + p1 + p2 + p3 + plot_layout(ncol = 2, heights = c(2, 2)), 
         width = 15, height = 12, 
         filename = paste0('Coembed_Results2_v2/coembed_', sampleID0, '.png'), 
         device = 'png')


message('Call metacells Done! Saving metadata matrices...')

## construct and save metacell rna and atac matrices 
merged_cells = metacell_obj$cells_merged
nrna_cells = sapply(merged_cells, function(x) stringr::str_count(x, pattern = 'RNA_'))
natac_cells = sapply(merged_cells, function(x) stringr::str_count(x, pattern = 'ATAC_'))

write(merged_cells, file = paste0('Coembed_Results2_v2/metacells_names_', sampleID0, '.txt'))


#calculate metacell mtx
sele.metacells = merged_cells[nrna_cells >= 5 & natac_cells >= 5] ## ignore metacells which are mostly comprised of single modality
nrna_cells = nrna_cells[names(sele.metacells)]
natac_cells = natac_cells[names(sele.metacells)]

if(T){
  ## manually
  mask <- sapply(names(nrna_cells), function(x) colnames(coembed) %in%
                   unlist(strsplit(merged_cells[x], ',', fixed = T)))
  rownames(mask) = colnames(coembed)
  rna.mask <- sapply(names(nrna_cells), function(x) mask[, x]/nrna_cells[x])
  rna.mtx <- coembed@assays$RNA@data %*% rna.mask # remember atac cells in RNA assay have 0 account
  
  atac.mask <- sapply(names(natac_cells), function(x) mask[, x]/natac_cells[x])
  atac.mtx <- coembed@assays$ATAC@data %*% atac.mask[colnames(coembed@assays$ATAC@data), ]

}else{
  rna.mtx <- metacell_obj@assays$RNA@data
  atac.mtx <- metacell_obj@assays$ATAC@data
  atac.mtx = atac.mtx[, names(sele.metacells)]
  rna.mtx = rna.mtx[, naems(sele.metacells)]
}

colnames(rna.mtx) = colnames(atac.mtx) = paste0(sampleID0, '_', names(sele.metacells))

saveRDS(rna.mtx, file = paste0('Coembed_Results2_v2/rna_metacell_mtx_v2_noTFIDF_', sampleID0, '.rds'))
saveRDS(atac.mtx, file = paste0('Coembed_Results2_v2/atac_metacell_mtx_v2_noTFIDF_', sampleID0, '.rds'))
message(paste('#of final metacells:', ncol(rna.mtx)))

message('All Done!')
