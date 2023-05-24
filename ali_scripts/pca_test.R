setwd('D:/Projects/taxus/')

require(reshape2)
require(ggplot2)
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

read_exp = function(inname)
{
  tpm = read.csv(inname, sep='\t')
  rownames(tpm) = tpm[,1]
  tpm = tpm[, -1]
  return(tpm)
}

kallisto_tpm = read_exp('kallisto_trimmed.txt')

ids = read.csv('data/names_projects_2.csv')

kallisto_log_tpm = log(1+kallisto_tpm)

inc_idx =  ids$Tissue %in% c('Bark', 'Cone', 'Leaf', 'Root', 'Stem', 'Twig') & 
          !ids$Sample %in% c('SRR17673884', 'SRR17673894', 'SRR17673892')

ids_sub = ids[inc_idx, ]
kallisto_log_tpm_sub = kallisto_log_tpm[, inc_idx]

get_per_proj_zm = function(kallisto_log_tpm, ids)
{
  kallisto_log_tpm_zm = c()
  ids_reorder = c()
  
  uprojects = unique(ids$Project)
  
  for (i in 1:length(uprojects))
  {
    p = uprojects[i]
    temp = kallisto_log_tpm[, ids$Project == p]
    ids_temp = ids[ids$Project == p, ]
    print(dim(temp))
    m = apply(temp, 1, mean)
    m = matrix(rep(m, ncol(temp)), ncol = ncol(temp))
    if (i == 1)
    {
      kallisto_log_tpm_zm = temp - m
      ids_reorder = ids_temp
    }
    else
    {
      kallisto_log_tpm_zm = cbind(kallisto_log_tpm_zm, temp - m)
      ids_reorder = rbind(ids_reorder, ids_temp)
    }
  }
  return(list(data=kallisto_log_tpm_zm, ids=ids_reorder))
}

get_per_type_zm = function(kallisto_log_tpm, ids)
{
  kallisto_log_tpm_zm = c()
  ids_reorder = c()
  
  uprojects = unique(ids[, c('Project', 'Genus', 'Tissue')])
  print(uprojects)
  
  for (i in 1:nrow(uprojects))
  {
    p = uprojects[i, 1]
    g = uprojects[i, 2]
    t = uprojects[i, 3]
    print(c(p, g, t))
    count = sum(ids$Project == p & ids$Genus == g & ids$Tissue == t)
    print(count)
    if (count == 1)
    {
      next
    }
    temp = kallisto_log_tpm[, ids$Project == p & ids$Genus == g & ids$Tissue == t]
    ids_temp = ids[ids$Project == p & ids$Genus == g & ids$Tissue == t, ]
    
    print(dim(temp))
    m = apply(temp, 1, mean)
    m = matrix(rep(m, ncol(temp)), ncol = ncol(temp))
    if (i == 1)
    {
      kallisto_log_tpm_zm = temp - m
      ids_reorder = ids_temp
    }
    else
    {
      kallisto_log_tpm_zm = cbind(kallisto_log_tpm_zm, temp - m)
      ids_reorder = rbind(ids_reorder, ids_temp)
    }
  }
  return(list(data=kallisto_log_tpm_zm, ids=ids_reorder))
}



res = get_per_proj_zm(kallisto_log_tpm_sub, ids_sub)
kallisto_log_tpm_sub_zm_per_proj = res$data
ids_per_proj = res$ids


res = get_per_type_zm(kallisto_log_tpm_sub, ids_sub)
kallisto_log_tpm_sub_zm_per_type = res$data
ids_per_type = res$ids

X1 = as.data.frame(t(kallisto_log_tpm_sub_zm_per_proj))
X1 = cbind(Project=ids_per_proj$Project, 
           Genus=ids_per_proj$Genus, 
           Tissue = ids_per_proj$Tissue, 
           X1)
X1 = X1[X1$Project != 'PRJNA427840', ]

X2 = as.data.frame(t(kallisto_log_tpm_sub_zm_per_type))
X2 = cbind(Project=ids_per_type$Project, 
           Genus=ids_per_type$Genus, 
           Tissue=ids_per_type$Tissue, 
           X2)

X3 = as.data.frame(t(kallisto_log_tpm_sub))
X3 = cbind(Project=ids_sub$Project, 
           Genus=ids_sub$Genus, 
           Tissue=ids_sub$Tissue, 
           X3)
X3 = X3[X3$Project != 'PRJNA427840', ]

dim(X2)

make_pca = function(X, outname, tstr)
{
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

X1_sub_1 = X1 %>% filter(Genus == 'Taxus Chinensis' & Tissue == 'Leaf')
X2_sub_1 = X2 %>% filter(Genus == 'Taxus Chinensis' & Tissue == 'Leaf')
X3_sub_1 = X3 %>% filter(Genus == 'Taxus Chinensis' & Tissue == 'Leaf')

X1_sub_2 = X1 %>% filter(Genus == 'Taxus Chinensis' & Tissue %in% c('Leaf', 'Root'))
X2_sub_2 = X2 %>% filter(Genus == 'Taxus Chinensis' & Tissue %in% c('Leaf', 'Root'))
X3_sub_2 = X3 %>% filter(Genus == 'Taxus Chinensis' & Tissue %in% c('Leaf', 'Root'))

X1_sub_3 = X1 %>% filter(Genus == 'Taxus Chinensis')
X2_sub_3 = X2 %>% filter(Genus == 'Taxus Chinensis')
X3_sub_3 = X3 %>% filter(Genus == 'Taxus Chinensis')


make_pca(X1_sub_1, 'new_pca_plots/taxus_chinensis_leaf_zm_per_proj', 'Taxus Chinensis, Leaf (zm per Project)')
make_pca(X2_sub_1, 'new_pca_plots/taxus_chinensis_leaf_zm_per_type', 'Taxus Chinensis, Leaf (zm per Tissue)')
make_pca(X3_sub_1, 'new_pca_plots/taxus_chinensis_leaf_no_zm', 'Taxus Chinensis, Leaf (no zm)')
make_pca(X1_sub_2, 'new_pca_plots/taxus_chinensis_leaf_root_zm_per_proj', 'Taxus Chinensis, Leaf+Root (zm per Project)')
make_pca(X2_sub_2, 'new_pca_plots/taxus_chinensis_leaf_root_zm_per_type', 'Taxus Chinensis, Leaf+Root (zm per Tissue)')
make_pca(X3_sub_2, 'new_pca_plots/taxus_chinensis_leaf_root_no_zm', 'Taxus Chinensis, Leaf+Root (no zm)')
make_pca(X1_sub_3, 'new_pca_plots/taxus_chinensis_zm_per_proj', 'Taxus Chinensis (zm per Project)')
make_pca(X2_sub_3, 'new_pca_plots/taxus_chinensis_zm_per_type', 'Taxus Chinensis (zm per Tissue)')
make_pca(X3_sub_3, 'new_pca_plots/taxus_chinensis_no_zm', 'Taxus Chinensis (no zm)')
make_pca(X1,       'new_pca_plots/all_zm_per_proj', 'All (zm per Project)')
make_pca(X2,       'new_pca_plots/all_zm_per_type', 'All (zm per Tissue)')
make_pca(X3,       'new_pca_plots/all_no_zm', 'All (no zm)')


X2_PRJNA797697 = X2 %>% filter(Project == 'PRJNA797697')
make_pca(X2_PRJNA797697, 'new_pca_plots/PRJNA797697_zm_per_type', 'PRJNA797697 (zm per Tissue)')

PRJNA797697_log_tpm = kallisto_log_tpm[, ids$Project == 'PRJNA797697']
PRJNA797697_cor = cor(PRJNA797697_log_tpm)

pdf('new_pca_plots/PRJNA79769_heatmap.pdf', width = 6, height = 6)
heatmap.2(PRJNA797697_cor, trace = 'none')
dev.off()

PRJNA797697_log_tpm = kallisto_log_tpm[, ids$Project == 'PRJNA797697' & 
                                         !ids$Sample %in% c('SRR17673884', 'SRR17673894', 'SRR17673892')]
PRJNA797697_cor = cor(PRJNA797697_log_tpm)

pdf('new_pca_plots/PRJNA79769_heatmap_excluded.pdf', width = 6, height = 6)
heatmap.2(PRJNA797697_cor, trace = 'none')
dev.off()


class(X1$Project)

table(X1$Project)


X3[1:4, 1:5]

X1 = X1[, -c(2,3)]
X1$Project = as.factor(X1$Project)
X2 = X2[, -c(2,3)]
X2$Project = as.factor(X2$Project)
X3 = X3[, -c(2,3)]
X3$Project = as.factor(X3$Project)

table(X1$Project)

table(ids$Project)


