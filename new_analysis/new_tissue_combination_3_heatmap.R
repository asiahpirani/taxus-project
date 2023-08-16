
args = commandArgs(trailingOnly=TRUE)
input_set = as.numeric(args[1])

source('load_and_save_functions.R')
source('apply_normalization.R')
source('apply_batch_effect.R')
source('load_and_save_functions.R')

library(dplyr)

normalize_and_correct <- function(data, annot, n_method, c_method)
{
  n_data = normalize(data, n_method)
  c_data = batch_correction(n_data, c_method, annot)
  return(c_data)
}


kallisto_raw = read_exp('../data/kallisto_trimmed_est_counts.txt')

ids = read.csv('../data/names_projects_2.csv')

select_subset = function(in_data, in_annot, projects, tissues, outname)
{
  annot = in_annot %>% filter(Project %in% projects &
                                Tissue %in% tissues)
  data  = in_data[, annot$Sample]
  save_data(data, annot, outname)
}

project_coms = list(c('PRJNA661543',                'PRJNA730337'),
                    c('PRJNA661543',                'PRJNA730337'),
                    c('PRJNA661543',                'PRJNA730337'),
                    c('PRJNA661543', 'PRJNA493167', 'PRJNA730337'),
                    c('PRJNA661543', 'PRJNA493167'               ),
                    c('PRJNA661543', 'PRJNA493167'               ),
                    c('PRJNA661543', 'PRJNA493167'               ))

combinations = list(c('Bark', 'Leaf'),         # PRJNA661543             PRJNA730337
                    c('Bark', 'Leaf', 'Root'), # PRJNA661543             PRJNA730337
                    c('Bark', 'Root'),         # PRJNA661543             PRJNA730337
                    c('Leaf', 'Root'),         # PRJNA661543 PRJNA493167 PRJNA730337
                    c('Leaf', 'Root', 'Stem'), # PRJNA661543 PRJNA493167
                    c('Leaf', 'Stem'),         # PRJNA661543 PRJNA493167
                    c('Root', 'Stem'))         # PRJNA661543 PRJNA493167

c_method_names = c('none', 'removeBatchEffect', 'combat', 'zm', 'zm_sd')
n_method_names = c('quantile', 'rlog', 'cpm', 'tmm', 'med', 'uq', 'none')


require(ComplexHeatmap)
require(reshape2)
require(ggplot2)
require(ggpubr)
require(ape)
require(tidyr)
require(dplyr)
require(grid)
require(ComplexHeatmap)
require(circlize)
require(PupillometryR)

require(topGO)
require(ALL)
require(gridExtra)
require(genefilter)
require(dplyr)

require(FactoMineR)
require(factoextra)
require(ggrepel)
require(ggfortify)

require(gplots)

make_pca = function(data, annot, outname, tstr)
{
  X = t(data)
  X = data.frame(X)
  X = cbind(Project = annot$Project, Genus = annot$Genus, Tissue=annot$Tissue, X)

  m1 = apply(X[, -c(1, 2, 3)], 2, min)
  m2 = apply(X[, -c(1, 2, 3)], 2, max)
  idx = c(T, T, T, m1 != m2)
  XX = X[, idx]

  pca.res <- prcomp(XX[, -c(1, 2, 3)],
                    center = TRUE,
                    scale. = TRUE)

  autoplot(pca.res,
           data = XX,
           colour = 'Project') +
    ggtitle(tstr)
  ggsave(filename = paste(outname,'_color_by_project.pdf', sep=''), width = 6, height = 4)

  autoplot(pca.res,
           data = XX,
           colour = 'Tissue') +
    ggtitle(tstr)
  ggsave(filename = paste(outname,'_color_by_tissue.pdf', sep=''), width = 6, height = 4)

  autoplot(pca.res,
           data = XX,
           colour = 'Genus') +
    ggtitle(tstr)
  ggsave(filename = paste(outname,'_color_by_genus.pdf', sep=''), width = 6, height = 4)
}

make_heat = function(data, ks, outname, tstr)
{
  dist_mat   <- dist(data, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'average')

  for (k in ks)
  {
    pdf(file = paste0(outname, '_k', k, '.pdf'), width = 6, height = 8)
    ht = Heatmap(as.matrix(data), name = "mat",
                 row_order = hclust_avg$order,
                 top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 1+(1:k)))),
                 column_km = k,
                 row_labels = rep('', nrow(data)),
                 column_title = tstr)
    draw(ht)
    dev.off()
  }
}


for (i in input_set:input_set)
{
  tissues = combinations[[i]]
  tissues_tstr = paste0(tissues, collapse = '+')
  fname = paste(tissues, collapse = '_')
  print(fname)

  outdir  = paste0('new_tissue_combinations/',fname)
  outname = paste0(outdir,'/kallisto_',fname)

  load(paste0(outname,'_raw.RData'))
  raw_data = data
  m1 = apply(raw_data, 1, min)
  m2 = apply(raw_data, 1, max)

  for (c_method in c_method_names)
  {
    for (n_method in n_method_names)
    {
      print(paste(fname, c_method, n_method))

      suboutdir = paste0(outdir,'/',n_method,'_',c_method)
      load(paste0(suboutdir,'/kallisto_',fname,'_',n_method,'_',c_method,'_log10.RData'))
      nrm_data = data

      outname = paste0(suboutdir,'/kallisto_',fname,'_',n_method,'_',c_method)

      k1 = length(unique(annot$Project))
      k2 = length(unique(annot$Tissue))
      ks = min(k1,k2):(k1*k2)

      sub_nrm_data = nrm_data[m1!=m2, ]
      colnames(sub_nrm_data) = paste(annot$Sample, annot$Project, annot$Tissue, sep='.')

      make_heat(sub_nrm_data, ks, paste0(outname, '_heatmap_all'),
                paste0('Normalize:',n_method,'+Batch:',c_method,', ',tissues_tstr))
    }
  }
}

