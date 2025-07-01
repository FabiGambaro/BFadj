inferred.root <- function(folders, input_dir, input_dir_tipSWap, scheme_500,BF_threshold,nberOfSamples) {
  
  result <- list() 
  result.true.phylo <- list() 
  final.results <- list()
  
  for (h in 1:length(folders)) {
    print(paste0("Working on ", folders[[h]]))
    MCC.tipSwap.to.analyze <- list.files(paste0(folders[h], "/", input_dir_tipSWap), pattern = "MCC.tree")
    
    if (scheme_500) {
      MCC.tipSwap.to.analyze <- MCC.tipSwap.to.analyze[grepl("500", MCC.tipSwap.to.analyze)]
      nb.seq <- 500
    } else {
      MCC.tipSwap.to.analyze <- MCC.tipSwap.to.analyze[grepl("150", MCC.tipSwap.to.analyze)]
      nb.seq <- 150
    }
    
    for (j in 1:length(MCC.tipSwap.to.analyze)) {
      name <- gsub("_tipSwap.MCC.tree", "", MCC.tipSwap.to.analyze[[j]])
      print(paste0("Working on ", name))
      TP <- 0; FP <- 0; adjTP <- 0; adjFP <- 0; adjTN <- 0; TN <- 0; FN <- 0; adjFN <- 0
      
      log <- read_delim(paste0(folders[h], "/", input_dir, name, ".log.txt"), delim = "\t",
                        escape_double = FALSE, trim_ws = TRUE, skip = 3, progress = F, show_col_types = F)
      log1 <- log[,grepl("regions.indicators.", colnames(log))]
      D <- (1 + sqrt(1 + 4 * ncol(log1))) / 2
      
      tree <- readAnnotatedNexus(paste0(folders[h], "/", input_dir, name, ".mcc.tree"))
      df <- data.frame(region = unlist(tree$root.annotation$regions.states.set), 
                       prob = unlist(tree$root.annotation$regions.states.set.prob))
      
      tree.swap <- readAnnotatedNexus(paste0(folders[h], "/", input_dir_tipSWap, name, "_tipSwap.MCC.tree"))
      df.swap <- data.frame(region = unlist(tree.swap$root.annotation$region.set),
                            prob.swap = unlist(tree.swap$root.annotation$region.set.prob))
      
      df <- dplyr::left_join(df, df.swap, by = c("region" = "region"))
      df$prob.swap[is.na(df$prob.swap)] <- 1/nberOfSamples
      #if any of p values is equal to one, modify to avoid an error or undefined result
      df$prob[df$prob==1] <- (nberOfSamples - 1) / nberOfSamples
      df$prob.swap[df$prob.swap==1] <- (nberOfSamples - 1) / nberOfSamples
      
      q = (1/D)
      #q <- (log(2) + (D - 1)) / (D * (D - 1) / 1)
      df$BF <- df$prob / (1 - df$prob) / (q / (1 - q))
      df$adjBF <- df$prob / (1 - df$prob) / (df$prob.swap / (1 - df$prob.swap))
      #df$diff <- df$adjBF - df$BF
      
      result[[name]] <- df
      
      true.phylo <- readAnnotatedNexus(paste0(folders[h], "/files/", name, "_modi.nex"))
      true.root <- true.phylo$root.annotation
      df.true.phylo <- data.frame(root = true.root)
      
      result.true.phylo[[name]] <- df.true.phylo
      
      inferred.root <- df$region[which(df$prob == max(df$prob))]
      simulated.root <- unique(df.true.phylo$regions.states)
      
      if (df.true.phylo$regions.states %in% df$region) {
        if (simulated.root == inferred.root) {
          if (!is.na(df$BF[which(df$region == inferred.root)]) & df$BF[which(df$region == inferred.root)] >= BF_threshold) {
            TP <- 1
          } else {
            FN <- 1
          }
        } else {
          if (!is.na(df$BF[which(df$region == inferred.root)]) & df$BF[which(df$region == inferred.root)] >= BF_threshold) {
            FP <- 1
          } else {
            TN <- 1
          }
        }
        
        if (simulated.root == inferred.root) {
          if (!is.na(df$adjBF[which(df$region == inferred.root)]) & df$adjBF[which(df$region == inferred.root)] >= BF_threshold) {
            adjTP <- 1
          } else {
            adjFN <- 1
          }
        } else {
          if (!is.na(df$adjBF[which(df$region == inferred.root)]) & df$adjBF[which(df$region == inferred.root)] >= BF_threshold) {
            adjFP <- 1
          } else {
            adjTN <- 1
          }
        }
      } else {
        if (!is.na(df$BF[which(df$region == inferred.root)]) & df$BF[which(df$region == inferred.root)] >= BF_threshold) {
          FP <- 1
        } else {
          TN <- 1
        }
        if (!is.na(df$adjBF[which(df$region == inferred.root)]) & df$adjBF[which(df$region == inferred.root)] >= BF_threshold) {
          adjFP <- 1
        } else {
          adjTN <- 1
        }
      }
      
      dt <- data.frame(TP = TP, adjTP = adjTP, FP = FP, adjFP = adjFP, TN = TN, adjTN = adjTN, FN = FN, adjFN = adjFN)
      final.results[[name]] <- dt
    }
  }
  
  return(list(result = result, result.true.phylo = result.true.phylo, final.results = final.results))
}

