## get a binary matrix indicates the gene-peak affinity
## gene.list are data.table including column name gene_name 
## gene_ann should include gene_name,chr,start,end


get_gene2peak_map <- function(gene.list, peak_names,
                              gene_ann, distal_dist = 2e05){


  colnames(peak_names) <- "peak_name"
  colnames(gene.list) <- "gene_name"

  print(head(peak_names))
  peaks <- separate(
     data.frame(peak_names),
     col = peak_name,
     into = c("chr", "start", "end"),
     remove = FALSE,
     sep = "-"
  )
  #peaks <- tidyr::separate(data.tab1le(peak_names), 
  #                         col = peak_name, remove = F,
  #                         into = c('chr', 'start', 'end'))
  class(peaks$start) = class(peaks$end) = 'integer'
  setkey(gene_ann, gene_name)
  peaks <- as.data.table(peaks)
  setkey(peaks, peak_name)
  gene.list = gene.list[gene_name %in% gene_ann$gene_name]
  gene.list$chr = gene_ann[gene.list$gene_name]$chr
  gene.list$start = gene_ann[gene.list$gene_name]$start
  gene.list$end = gene_ann[gene.list$gene_name]$end

   ## for each gene, get the corresponding peaks
  gene2peaks = lapply(gene.list$gene_name, function(x) {

    gene_info = gene_ann[.(x), on = .(gene_name)]

    chr0 = gene_info$chr #gene_ann[x]$chr
    start0 = gene_info$start #gene_ann[x]$start
    end0 = gene_info$end #gene_ann[x]$end

    peaks0 = peaks[chr == chr0]
    peaks0 = peaks0[abs(start/2 + end/2 - start0/2 - end0/2) <= distal_dist]
    return(peaks0$peak_name)
  } )

  ## pool all peaks relate to one gene #for some reason sort does not work on a list?
  gene2peaks.u <- lapply(sort(unlist(unique(gene.list$gene_name))), function(x){
    id = which(gene.list$gene_name == x)
    tmp_peak <- do.call('c', lapply(id, function(x) gene2peaks[[x]]))
    return(tmp_peak)
  })
  #print(gene.list)

  names(gene2peaks.u) <- sort(unlist(unique(gene.list$gene_name)))  #for some reason sort does not work on a list?
  lens = sapply(gene2peaks.u, length)

  genes.f <- names(which(lens > 0))
  lens = lens[lens > 0]

  ## construct overlap matrix
  gene2peaks.dt <- data.table('gene' = rep(genes.f, lens),
                              'peak' = do.call('c', lapply(genes.f,
                                                           function(x) gene2peaks.u[[x]])))
  upeaks = sort(unique(gene2peaks.dt$peak))
  gene2peaks.dt[, 'id1' := which(genes.f == gene), by = gene]
  gene2peaks.dt[, 'id2' := which(upeaks == peak), by = peak]
  gene2peak.map <- sparseMatrix(i = gene2peaks.dt$id1,
                                j = gene2peaks.dt$id2,
                                dimnames = list(genes.f, upeaks))
  gene2peak.map = gene2peak.map * 1

  return(gene2peak.map)
}

## annotate peaks with gene +/- 5kb of its TSS
# input peak_coords with chr-start-end, format
annPeak2Gene <- function(peak_coords, gene_ann, proxim_dist = 5e+03){
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  peaks = tidyr::separate(data.table(x = peak_coords),
                          col = x,
                          into = c('chr', 'start', 'end'))


  peaks$peak_name = peak_coords
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'

  chrs = unique(peaks$chr)
  peaks <- as.data.table(peaks)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    if (chr0 == "chrM"){
      genes0 = gene_ann[chr == "chrMT"]
    } else {
      genes0 = gene_ann[chr == chr0]
    }

    peaks0$gene_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= (peaks0$end[i] + proxim_dist) & tss >= (peaks0$start[i] - proxim_dist)]
      if(nrow(tss0) > 0 ) {
        peaks0$gene_name[i] = paste(unique(tss0$gene_name), collapse = ',')
      }
    }

    peaks_ann = rbind(peaks_ann, peaks0)
    }

    peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1,
                                        paste0(peak_name, ',', gene_name), peak_name)]


    setkey(peaks_ann, peak_name)

    return(peaks_ann)
}


## map gene to overlapping atac peak
## gene_list with genename, chr, start, end
geneOverlapPeak <- function(gene_list, peak_names, mid_dist = 1000){
  # should include tss information in gene_list
  peaks = tidyr::separate(data = data.table('peak_name' = peak_names),
                          col = peak_name, into = c('chr', 'start', 'end'),
                          remove = F)
  class(peaks$chr) = 'character'
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'

  chrs = unique(gene_list$chr)
  gene_new = NULL
  peaks <- as.data.table(peaks)
  peaks[, 'midP' := start/2 + end/2]
  for(chr0 in chrs){
    gene0 = gene_list[chr == chr0, ]
    gene0$peak_name = 'Not_Found'
    peaks0 = peaks[chr == chr0]
    gene0[, 'peak_id0' := any( abs(peaks0$midP -start) < mid_dist | abs(peaks0$midP - end) < mid_dist),
          by = gene_name]
    gene1 = gene0[peak_id0 == FALSE]
    gene2 = gene0[peak_id0 == TRUE]
    gene2[, 'peak_id' :=  which.min(abs(peaks0$midP - start - 1000)), by = gene_name]

    gene2[, 'peak_name' :=  peaks0[peak_id]$peak_name, by = gene_name]
    gene2$peak_id = NULL
    gene_new = rbind(gene_new, gene1, gene2)

  }
  gene_new[, c('peak_id0') := NULL]
  return(gene_new)
}



