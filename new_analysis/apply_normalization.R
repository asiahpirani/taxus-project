
#######################################################################################################
#######################################################################################################
##                                                                                                   ##
## from here:                                                                                        ##
## https://github.com/alexisvdb/rnaseq_coexpression/blob/main/scripts/Rfunctions_normalization.R     ##
## from this paper:                                                                                  ##
## https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0263344                         ##
##                                                                                                   ##
#######################################################################################################
#######################################################################################################

library(limma)
library(DESeq2)
library(edgeR)

normalize = function(dat, method = "none"){
  
  # normalize
  dat <- switch (method,
                 quantile = normalization_quantile(dat), 
                 rlog     = normalization_rlog(dat),
                 cpm      = normalization_cpm(dat),
                 tmm      = normalization_tmm(dat),
                 med      = normalization_med(dat),
                 uq       = normalization_uq(dat),
                 none     = dat,
  )
  
  # set a small pseudo count
  pseudo <- quantile(dat[dat > 0], 0.01)
  
  # and convert to log10 values
  dat.log10 <- log10(dat+(10*pseudo))
  
  # return result
  dat.log10
}

normalization_uq = function(dat){

  # get the UQ of the NON-ZERO values
  uqs <- apply(dat, 2, function(x) quantile(x[x!=0],0.75))
  # the geometric mean
  mean.uq <- exp(mean(log(uqs))) 
  correction.factors <- mean.uq/uqs
  
  # make a DGEList object from the counts
  dat <- DGEList(counts=dat)
  dat$samples[,3] <- correction.factors
  
  # return the normalized tag counts
  cpm(dat)
}

normalization_med = function(dat){

  # get the median of the NON-ZERO values
  meds <- apply(dat, 2, function(x) median(x[x!=0]))
  # the geometric mean
  mean.med <- exp(mean(log(meds))) 
  correction.factors <- mean.med/meds
  
  # make a DGEList object from the counts
  dat <- DGEList(counts=dat)
  dat$samples[,3] <- correction.factors
  
  # return the normalized tag counts
  cpm(dat)
}

normalization_tmm = function(dat){
  
  # make a DGEList object from the counts
  dat <- DGEList(counts=dat)
  
  # normalization 
  # by default this is the TMM (trimmed mean of M-values)
  dat <- calcNormFactors(dat)
  cpm(dat)
}

normalization_cpm = function(dat){
  
  total_tag_count <- apply(dat,2,sum)
  for(c in 1:ncol(dat)){
    dat[,c] <- dat[,c]*1e6/total_tag_count[c]
  }
  # return result
  dat
}

normalization_rlog = function(dat){
  
  # # remove ridiculously large read counts in the ComBat-seq output (some >>> 1e30!!)
  # # DESeq2 can not handle them.
  # dat[dat>1e9] <- 1e9
  
  coldata <- matrix(NA, nrow=ncol(dat), ncol=1)
  dds <- DESeqDataSetFromMatrix(countData = dat,
                                colData = coldata,
                                design= ~ 1)
  # run this or DESeq() first
  dds <- estimateSizeFactors(dds) 
  
  # return the normalized tag counts
  counts(dds, normalized=T)
} 

normalization_quantile = function(dat){
  
  normalizeQuantiles(dat)
}

