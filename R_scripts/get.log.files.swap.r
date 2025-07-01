get.log.files.swap <- function(MainDirectory, Folder, scheme, burnin_fraction,nberOfSamples) {
  output <- list()
  ListFilesToAnalyze <- list.files(path = paste0(MainDirectory, Folder, "/dta_tipSwap/"), pattern = "*tipSwap.log")
  
  if (scheme == 500) {
    ListFilesToAnalyze <- ListFilesToAnalyze[!grepl("150", ListFilesToAnalyze)]
  }
  
  if (scheme == 150) {
    ListFilesToAnalyze <- ListFilesToAnalyze[!grepl("500", ListFilesToAnalyze)]
  }
  
  for (filename in ListFilesToAnalyze) {
    name <- gsub("*_tipSwap.log", "", filename)
    log <- read_tsv(file = paste0(MainDirectory, Folder, "/dta_tipSwap/", filename), 
                    comment = "#", col_names = T, show_col_types = F)
    burnIn <- round(burnin_fraction * dim(log)[1]) +1
    log <- log[burnIn:dim(log)[1],]
    #subsample 100
    index1 <- dim(log)[1]
    interval <- floor(index1/nberOfSamples)
    subsampled_log <- log[seq(1, index1, by = interval), ]
    subsampled_log <- subsampled_log[1:nberOfSamples, ] #make sure there are 100
    
    log1 <- subsampled_log[, grepl("regions.indicators.", colnames(subsampled_log))]
    output[[name]] <- log1
  }
  return(output)
}