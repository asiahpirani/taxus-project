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

require(randomForest)
library(e1071)
require(rsample)

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


predict_lbl_rf <- function(X)
{
  acc = c()
  for (i in 1:100)
  {
    split_data <- initial_split(X, prop = 0.5, strata = 'lbl')
    train_data <- training(split_data)
    test_data  <- testing(split_data)

    model = randomForest(lbl ~ ., data=train_data)
    res = predict(model, newdata = test_data)
    acc = c(acc, sum(res == test_data$lbl)/length(test_data$lbl))
  }
  return(acc)
}

predict_lbl_lr <- function(X)
{
  acc = c()
  for (i in 1:100)
  {
    split_data <- initial_split(X, prop = 0.5, strata = 'lbl')
    train_data <- training(split_data)
    test_data  <- testing(split_data)

    model = glm(lbl ~ ., data=train_data, family = 'binomial')
    res = predict(model, newdata = test_data)
    res = res>0.5
    acc = c(acc, sum(res == test_data$lbl)/length(test_data$lbl))
  }
  return(acc)
}

predict_lbl_svm_lin <- function(X)
{
  acc = c()
  for (i in 1:100)
  {
    split_data <- initial_split(X, prop = 0.5, strata = 'lbl')
    train_data <- training(split_data)
    test_data  <- testing(split_data)

    model = svm(lbl ~ ., data=train_data, kernel = "linear")
    res = predict(model, newdata = test_data)
    acc = c(acc, sum(res == test_data$lbl)/length(test_data$lbl))
  }
  return(acc)
}

predict_lbl_svm_rad <- function(X)
{
  acc = c()
  for (i in 1:100)
  {
    split_data <- initial_split(X, prop = 0.5, strata = 'lbl')
    train_data <- training(split_data)
    test_data  <- testing(split_data)

    model = svm(lbl ~ ., data=train_data, kernel = "radial")
    res = predict(model, newdata = test_data)
    acc = c(acc, sum(res == test_data$lbl)/length(test_data$lbl))
  }
  return(acc)
}

make_data_with_label = function(in_data, lbl)
{
  t_in_data <- t(in_data)
  colnames(t_in_data) <- paste0('G', 1:ncol(t_in_data))
  rownames(t_in_data) <- NULL
  data <- cbind(lbl=as.factor(lbl), data.frame(t_in_data))
  return(data)
}

run_one_sample_raw = function(clas_func, sub_raw_proj, sub_raw_tiss, sub_raw_gns)
{
  acc_raw_proj = clas_func(sub_raw_proj)
  acc_raw_tiss = clas_func(sub_raw_tiss)
  acc_raw_gns  = clas_func(sub_raw_gns)

  all_acc = rbind(data.frame(normalized=F, class='Project', accuracy=acc_raw_proj),
                  data.frame(normalized=F, class='Tissue',  accuracy=acc_raw_tiss),
                  data.frame(normalized=F, class='Genus',   accuracy=acc_raw_gns))

  return(all_acc)
}

run_one_sample_nrm = function(clas_func, sub_nrm_proj, sub_nrm_tiss, sub_nrm_gns)
{
  acc_nrm_proj = clas_func(sub_nrm_proj)
  acc_nrm_tiss = clas_func(sub_nrm_tiss)
  acc_nrm_gns  = clas_func(sub_nrm_gns)

  all_acc = rbind(data.frame(normalized=T, class='Project', accuracy=acc_nrm_proj),
                  data.frame(normalized=T, class='Tissue',  accuracy=acc_nrm_tiss),
                  data.frame(normalized=T, class='Genus',   accuracy=acc_nrm_gns))

  return(all_acc)
}

make_one_box = function(all_acc, outname, tstr)
{
  ggplot(all_acc, aes(x=normalized, fill=class, y=accuracy)) +
    geom_boxplot(lwd=0.25, outlier.size = .5) + ggtitle(tstr) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(0.3, 1.4)

  ggsave(filename = paste0(outname, '_accuracy_3.pdf'), width = 4, height = 5)
}


for (i in input_set:input_set)
{
  tissues = combinations[[i]]
  tissues_tstr = paste0(tissues, collapse = '+')
  fname = paste(tissues, collapse = '_')
  print(fname)

  outdir  = paste0('new_tissue_combinations/',fname)
  outname = paste0(outdir,'/kallisto_',fname)

  suboutdir = paste0(outdir,'/','none','_','none')
  load(paste0(suboutdir,'/kallisto_',fname,'_','none','_','none','_log10.RData'))
  raw_data = data
  m1 = apply(raw_data, 1, min)
  m2 = apply(raw_data, 1, max)
  
  rs = rowSums(raw_data)
  # sub_raw_data = raw_data[order(rs, decreasing = T)[1:100], ]
  sub_raw_data = raw_data
  colnames(sub_raw_data) = paste(annot$Sample, annot$Project, annot$Tissue, sep='.')

  sub_raw_proj = make_data_with_label(sub_raw_data, annot$Project)
  sub_raw_tiss = make_data_with_label(sub_raw_data, annot$Tissue)
  sub_raw_gns  = make_data_with_label(sub_raw_data, annot$Genus)

  ## Random Forest
  all_acc_raw_rf = run_one_sample_raw(predict_lbl_rf,
			   sub_raw_proj, sub_raw_tiss, sub_raw_gns)
  ### SVM Linear
  #all_acc_raw_svm_lin = run_one_sample_raw(predict_lbl_svm_lin,
  #		   sub_raw_proj, sub_raw_tiss, sub_raw_gns)
  ### SVM Radial
  #all_acc_raw_svm_rad = run_one_sample_raw(predict_lbl_svm_rad,
  #		   sub_raw_proj, sub_raw_tiss, sub_raw_gns)

  for (c_method in c_method_names)
  {
    for (n_method in n_method_names)
    {
      if (c_method == 'none' & n_method == 'none')
      {
	next
      }
      print(paste(fname, c_method, n_method))

      suboutdir = paste0(outdir,'/',n_method,'_',c_method)
      load(paste0(suboutdir,'/kallisto_',fname,'_',n_method,'_',c_method,'_log10.RData'))
      nrm_data = data
      
      # sub_nrm_data = nrm_data[order(rs, decreasing = T)[1:100], ]
      sub_nrm_data = nrm_data
      colnames(sub_nrm_data) = paste(annot$Sample, annot$Project, annot$Tissue, sep='.')

      sub_nrm_proj = make_data_with_label(sub_nrm_data, annot$Project)
      sub_nrm_tiss = make_data_with_label(sub_nrm_data, annot$Tissue)
      sub_nrm_gns  = make_data_with_label(sub_nrm_data, annot$Genus)

      ## Random Forest
      all_acc = run_one_sample_nrm(predict_lbl_rf,
			       sub_nrm_proj, sub_nrm_tiss, sub_nrm_gns)
      all_acc = rbind(all_acc_raw_rf, all_acc)
      outname = paste0(suboutdir,'/rf_100r_',fname,'_',n_method,'_',c_method)
      make_one_box(all_acc, outname, paste0('RF, Normalize:',n_method,'+Batch:',c_method,', ',tissues_tstr))

      ### SVM Linear
      #all_acc = run_one_sample_nrm(predict_lbl_svm_lin,
      #		       sub_nrm_proj, sub_nrm_tiss, sub_nrm_gns)
      #all_acc = rbind(all_acc_raw_svm_lin, all_acc)
      #outname = paste0(suboutdir,'/svm_lin_100r_',fname,'_',n_method,'_',c_method)
      #make_one_box(all_acc, outname, paste0('SVM (lin), Normalize:',n_method,'+Batch:',c_method,', ',tissues_tstr))

      ### SVM Radial
      #all_acc = run_one_sample_nrm(predict_lbl_svm_rad,
      #		       sub_nrm_proj, sub_nrm_tiss, sub_nrm_gns)
      #all_acc = rbind(all_acc_raw_svm_rad, all_acc)
      #outname = paste0(suboutdir,'/svm_rad_100r_',fname,'_',n_method,'_',c_method)
      #make_one_box(all_acc, outname, paste0('SVM (rad), Normalize:',n_method,'+Batch:',c_method,', ',tissues_tstr))
    }
  }
}
