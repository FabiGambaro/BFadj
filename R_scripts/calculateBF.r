calculateBF <- function(logFile1, logFile2, nberOfSamples) {
  
  colnames(logFile1) <- gsub("regions.indicators.", "", colnames(logFile1))
  D <- (1 + sqrt(1 + 4 * ncol(logFile1))) / 2
  q <- (log(2) + (D - 1)) / (D * (D - 1) / 1) #asymmetric q = (log(2)+K-1)/(K*(K-1))
  p <- colMeans(logFile1)
  
  #if any of p values is equal to one, modify to avoid an error or undefined result
  p[p == 1] <- (nberOfSamples - 1) / nberOfSamples
  
  BF <- p / (1 - p) / (q / (1 - q)) #Bayes Factor as the ratio of the posterior odds and prior odds
  
  colnames(logFile2) <- gsub("regions.indicators.", "", colnames(logFile2))
  p2 <- colMeans(logFile2)
  
  #if any of p values is equal to one, modify to avoid an error or undefined result
  p2[p2 == 1] <- (nberOfSamples - 1) / nberOfSamples
  
  BF2 <- p/(1-p) /(p2/(1-p2))
  
  result_df <- data.frame(p.standard =p,p.swap=p2, BF = BF, BFadj = BF2)

  return(result_df)
}

