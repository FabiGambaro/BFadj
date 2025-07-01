require(readr)
require(dplyr)
require(plyr)
require(treeio)
require(seraphim)
require(HDInterval)
  
MainDirectory = "~/Dropbox/BF_adjusted/SIMULATIONS/"

scheme = 150
scheme = 500
folders = list.files(pattern = "^simulation[0-9]{1,2}$")

source(paste0(MainDirectory, "scripts/trans.events.extractions.postTrees.r"))
source(paste0(MainDirectory, "scripts/trans.events.simTrees.r"))
source(paste0(MainDirectory, "scripts/get.log.files.r"))
source(paste0(MainDirectory, "scripts/get.log.files.swap.r"))
source(paste0(MainDirectory, "scripts/calculateBF.r"))
source(paste0(MainDirectory, "scripts/score.transitions.threshold.r"))

for (h in 1:length(folders)) {
  cat("Processing folder:", folders[h], "\n")
  inferred.transitions <- trans.events.extractions(MainDirectory, folders[h], scheme)
}

BF_thresholds <- c(3, 10, 20)
for (h in 2:length(folders)) {
  cat("Processing folder:", folders[h], "\n")
  
  LogFilesToAnalyze <- get.log.files(MainDirectory, folders[h], scheme, burnin_fraction = 0.1, nberOfSamples = 100)
  LogFilesToAnalyzeSwap <- get.log.files.swap(MainDirectory, folders[h], scheme, burnin_fraction = 0.1, nberOfSamples = 100)
  
  BF_BFadj_list <- list()
  options(scipen = 999); options(digits = 4)
  #estimate BF and adjBF using log files pre analyze before
  for (i in 1:length(LogFilesToAnalyze)) {
    log.name <- names(LogFilesToAnalyze)[i]
    log1 <- LogFilesToAnalyze[[log.name]]
    log2 <- LogFilesToAnalyzeSwap[[log.name]]
    BF_BFadj_list[[log.name]] <- calculateBF(log1, log2,nberOfSamples = 100)
  }
  
  simulated.transitions <- trans.events.simTrees(MainDirectory, folders[h], scheme)
  inferred.transitions <- list()
  for (i in 1:length(simulated.transitions)) {
    name <- names(simulated.transitions)[i]
    if(file.exists(paste0(MainDirectory, folders[h],"/dta/",name,"_ext/transitions.events.tsv"))){
      inferred.transitions[[name]] <- read_tsv(paste0(MainDirectory, folders[h],"/dta/",name,"_ext/transitions.events.tsv"), 
                                               show_col_types = F, progress = F) 
    }
}
  
  #merge the info from the inferred, simulated phylo and the BF estimates
  merged_list <- list()
  if(!file.exists(paste0(MainDirectory, folders[h], "/outputs/"))){
    dir.create(paste0(MainDirectory, folders[h], "/outputs/"))
    }
  
  for (i in 1:length(inferred.transitions)){
    df1_name <- names(inferred.transitions)[i]
    df1 <- inferred.transitions[[df1_name]]
    df2 <-  BF_BFadj_list[[df1_name]]
    df3 <- simulated.transitions[[df1_name]]
    df1$count <- df3$count[match(df1$transition, df3$transition)]
    df1$p.standard <- df2$p.standard[match(df1$transition, rownames(df2))]
    df1$BF <- df2$BF[match(df1$transition, rownames(df2))]
    df1$p.swap <- df2$p.swap[match(df1$transition, rownames(df2))]
    df1$BFadj <- df2$BFadj[match(df1$transition, rownames(df2))]
    write_tsv(df1, paste0(MainDirectory, folders[h], "/outputs/", df1_name, ".tsv"))
    merged_list[[df1_name]] <- df1
  }
  
  
  #on the final data frame estimate TP, TN, FP, FN considering different thresholds
  
  for (k in 1:length(BF_thresholds)){
    BF_threshold = BF_thresholds[k]
    
    if (!dir.exists(paste0(MainDirectory, folders[h], "/BF_results/"))) {
      dir.create(paste0(MainDirectory, folders[h], "/BF_results/"))
    }
    
    results <- score.transitions.threshold(merged_list,BF_threshold)
    results.all <- ldply(results, rbind)
    write_tsv(results.all, paste0(MainDirectory, folders[h], "/BF_results/BF_transitions_",scheme,
                                   "_",BF_threshold,".tsv"))
    #saveRDS(results_all,paste0("Transitions_events/all.transitions_", 
    #                            scheme,"_",BF_threshold,".rds"))
  }

}
  


