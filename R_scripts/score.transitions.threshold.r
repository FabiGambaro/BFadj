score.transitions.threshold <- function(mylist,BF_threshold) {
  tab_transitions <- list()
  
  for (i in 1:length(mylist)) {
    simulation <- names(mylist)[i]
    df <- mylist[[i]]
  
    df_TP <- 0; df_TN <- 0; df_FP <- 0; df_FN <- 0
    df_TPadj <- 0; df_TNadj <- 0; df_FPadj <- 0; df_FNadj <- 0
    
    for (j in 1:nrow(df)) {
      count <- df$count[j]
      mean_count <- df$mean[j]
      BF <- df$BF[j]
      BFadj <- df$BFadj[j]
      
      #BF
      if (!is.na(count)) {
        if (count == mean_count && !is.na(BF) && BF >= BF_threshold) {
          df_TP <- df_TP + count
        } else if (mean_count > count && !is.na(BF) && BF >= BF_threshold) {
          df_TP <- df_TP + min(count, mean_count)
          df_FP <- df_FP + (mean_count - count)
        } else if (mean_count < count && !is.na(BF) && BF >= BF_threshold) {
          df_TP <- df_TP + mean_count
          df_FN <- df_FN + (count - mean_count)
        } else if (is.na(BF) || BF < BF_threshold) {
          df_FN <- df_FN + count
        }
      } else {
        if (is.na(BF) || BF < BF_threshold) {
          df_TN <- df_TN + 1
        } else if (!is.na(BF) && BF >= BF_threshold) {
          df_FP <- df_FP + mean_count
        }
      }
      
      # BFadj
      if (!is.na(count)) {
        if (count == mean_count && !is.na(BFadj) && BFadj >= BF_threshold) {
          df_TPadj <- df_TPadj + count
        } else if (mean_count > count && !is.na(BFadj) && BFadj >= BF_threshold) {
          df_TPadj <- df_TPadj + min(count, mean_count)
          df_FPadj <- df_FPadj + (mean_count - count)
        } else if (mean_count < count && !is.na(BFadj) && BFadj >= BF_threshold) {
          df_TPadj <- df_TPadj + mean_count
          df_FNadj <- df_FNadj + (count - mean_count)
        } else if (is.na(BFadj) || BFadj < BF_threshold) {
          df_FNadj <- df_FNadj + count
        }
      } else {
        if (is.na(BFadj) || BFadj < BF_threshold) {
          df_TNadj <- df_TNadj + 1
        } else if (!is.na(BFadj) && BFadj >= BF_threshold) {
          df_FPadj <- df_FPadj + mean_count
        }
      }
    }
    
    # store results in a vector
    result <- c(TP = df_TP, TN = df_TN, FP = df_FP, FN = df_FN,
                TPadj = df_TPadj, TNadj = df_TNadj, FPadj = df_FPadj, FNadj = df_FNadj)

    tab_transitions[[simulation]] <- result
  }
  
  return(tab_transitions)
}
