
read_exp = function(inname)
{
  
  tpm = read.csv(inname, sep='\t')
  rownames(tpm) = tpm[,1]
  tpm = tpm[, -1]
  return(tpm)
}

save_data <- function(data, annot, fname)
{
  save(data, annot, file=fname)
}
