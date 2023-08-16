
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
## from here:                                                                                        ##
## https://github.com/alexisvdb/rnaseq_coexpression/blob/main/scripts/Rfunctions_batch_correction.R  ##
## from this paper:                                                                                  ##
## https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0263344                         ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

library(limma)
library(sva)

batch_correction = function(dat, method = "none", annotation){
  
  # c("none", "removeBatchEffect", "combat", "combat-seq")
  # batch correction
  dat <- switch (method,
                 removeBatchEffect = batch_effect_removeBatchEffect(dat, annotation),
                 combat            = batch_effect_combat(dat, annotation),
                 zm                = batch_effect_zm(dat, annotation),
                 zm_sd             = batch_effect_zm_std(dat, annotation),
                 # combat_seq        = batch_effect_combat_seq(dat, study.sets, annotation),
                 none              = dat,
  )
  
  # return result
  dat
}

batch_effect_zm = function(dat, annotation){
  
  cell.types <- annotation$Tissue
  cell.type.n <- length(unique(cell.types))
  
  batches <- annotation$Project
  batches.n <- length(unique(batches))
  
  subset.size <- length(unique(batches))
  
  indices <- 1:subset.size
  names(indices) <- unique(batches)
  batch.indices <- indices[batches]
  
  corrected.data <- dat
  
  for (b in unique(batches))
  {
    temp = dat[, batches == b]
    m = apply(temp, 1, mean)
    m = matrix(rep(m, ncol(temp)), ncol = ncol(temp))
    corrected.data[, batches == b] = temp - m
  }
  
  corrected.data
}

batch_effect_zm_std = function(dat, annotation){
  
  cell.types <- annotation$Tissue
  cell.type.n <- length(unique(cell.types))
  
  batches <- annotation$Project
  batches.n <- length(unique(batches))
  
  subset.size <- length(unique(batches))
  
  indices <- 1:subset.size
  names(indices) <- unique(batches)
  batch.indices <- indices[batches]
  
  corrected.data <- dat
  
  for (b in unique(batches))
  {
    temp = dat[, batches == b]
    m = apply(temp, 1, mean)
    m = matrix(rep(m, ncol(temp)), ncol = ncol(temp))
    
    s = apply(temp, 1, sd)
    s[s==0] = 1
    
    s = matrix(rep(s, ncol(temp)), ncol = ncol(temp))
    
    corrected.data[, batches == b] = (temp - m)/s
  }
  
  corrected.data
}


batch_effect_removeBatchEffect = function(dat, annotation){
  
  cell.types <- annotation$Tissue
  cell.type.n <- length(unique(cell.types))
  
  batches <- annotation$Project
  batches.n <- length(unique(batches))
  
  subset.size <- length(unique(batches))
  
  indices <- 1:subset.size
  names(indices) <- unique(batches)
  batch.indices <- indices[batches]
  
  # there are several different cases, requiring different processing:
  # - subset of samples including > 1 cell type
  # - subset of samples of only 1 batch and 1 cell type -> no batch correction needed
  # - subset of samples of > 1 batch but all of the same 1 cell type 
  
  if(cell.type.n > 1){ # if there is more than 1 cell type in the data
    
    mod = model.matrix(~as.factor(Tissue), data=annotation)
    corrected.data = removeBatchEffect(dat, batch = batches, design = mod)
    
  } else if (cell.type.n == 1 & batches.n == 1){ # if there is only 1 cell type in just 1 batch
    
    corrected.data = dat # nothing to be done
    
  } else if (cell.type.n == 1){ # if there is only 1 cell type in the data
    
    mod = model.matrix(~1, data=annotation)
    corrected.data = removeBatchEffect(dat, batch = batches, design = mod)
  }
  
  # return the result
  corrected.data
}

batch_effect_combat = function(dat, annotation){

  cell.types <- annotation$Tissue
  cell.type.n <- length(unique(cell.types))
  
  batches <- annotation$Project
  batches.n <- length(unique(batches))
  
  subset.size <- length(unique(batches))
  
  indices <- 1:subset.size
  names(indices) <- unique(batches)
  batch.indices <- indices[batches]

  # there are several different cases, requiring different processing:
  # - subset of samples including > 1 cell type
  # - subset of samples of only 1 batch and 1 cell type -> no batch correction needed
  # - subset of samples of > 1 batch but all of the same 1 cell type 
  
  if(cell.type.n > 1){ # if there is more than 1 cell type in the data
    
    mod = model.matrix(~as.factor(Tissue), data=annotation)
    corrected.data = ComBat(dat=as.matrix(dat), batch=batch.indices, mod=mod, par.prior=TRUE, prior.plots=FALSE)
    
  } else if (cell.type.n == 1 & batches.n == 1){ # if there is only 1 cell type in just 1 batch
    
    corrected.data = dat # nothing to be done
    
  } else if (cell.type.n == 1){ # if there is only 1 cell type in the data
    
    mod = model.matrix(~1, data=annotation)
    corrected.data = ComBat(dat=as.matrix(dat), batch=batch.indices, mod=mod, par.prior=TRUE, prior.plots=FALSE)
  }
  
  # return the result
  corrected.data
}
