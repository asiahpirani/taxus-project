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

for (i in 1:length(combinations))
{
  tissues = combinations[[i]]
  projects = project_coms[[i]]
  fname   = paste(tissues, collapse = '_')
  print(fname)
  outdir  = paste0('new_tissue_combinations/',fname)
  dir.create(file.path(outdir), showWarnings = F)

  outname = paste0(outdir,'/kallisto_',fname)
  select_subset(kallisto_raw, ids,
                projects,
                tissues,
                paste0(outname,'_raw.RData'))

  load(paste0(outname,'_raw.RData'))

  for (c_method in c_method_names)
  {
    for (n_method in n_method_names)
    {
      print(paste(fname, c_method, n_method))

      suboutdir = paste0(outdir,'/',n_method,'_',c_method)
      dir.create(file.path(suboutdir), showWarnings = F)

      nrm_data = normalize_and_correct(round(data), annot, n_method, c_method)
      save_data(nrm_data, annot,
                paste0(suboutdir,'/kallisto_',fname,'_',n_method,'_',c_method,'_log10.RData'))
    }
  }
}

