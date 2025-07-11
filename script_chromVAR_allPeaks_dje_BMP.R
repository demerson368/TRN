library(Seurat)

##reture chromvar result given a seurat obj

source('scDataAnalysis_Utilities_simp.R')
#args = commandArgs(T)
seuratPath = '/mnt/isilon/tan_lab/xuj5/ETP_ALL/Final_RScripts/scATAC/14_import_r_nr_add_metadata/t.all.40.objects/t.all.40.atac.blasts.1146.each.rds'#args[1]
chromvarPath = 'MetaData2/chromVar_output.rds'#args[2]
ndownSample = 150000 #as.integer(args[3])

seurat.obj = readRDS(seuratPath)

mtx = seurat.obj@assays$ATAC@counts

## further filter peaks
rs = Matrix::rowSums(mtx > 0)
mtx = mtx[rs > 10, ]

rn = rownames(mtx)
new_rn = sapply(rn, function(x) unlist(strsplit(x, ','))[1])

print("downsample upper limit: ")
print(ncol(mtx))

rownames(mtx) = new_rn
if(ndownSample < ncol(mtx)) {
    set.seed(2020)
    sele.cells = sample((1:ncol(mtx)), ndownSample)
    mtx = mtx[, sele.cells]
}
chromVAR.obj = run_chromVAR(mtx, genomeName = 'BSgenome.Hsapiens.UCSC.hg38')

saveRDS(chromVAR.obj, file = chromvarPath) 
