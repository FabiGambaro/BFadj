scripts_dir = "~/Dropbox/seraphim/"
source(paste0(scripts_dir,"DTA_tree_extraction_regions.r")) 

trans.events.extractions <- function(MainDirectory, Folder, scheme) {
  print(paste0("Working on folder: ", Folder))
  
  ListTreesToAnalyze <- list.files(path = paste0(MainDirectory, Folder, "/dta/"), pattern = "*_select100.trees")
  
  ListTreesToAnalyze <- list.files(path = paste0(
    MainDirectory,Folder,"/dta/"), pattern = "*_select100.trees")
  
  if (scheme == 500) {
    ListTreesToAnalyze <- ListTreesToAnalyze[!grepl("150", ListTreesToAnalyze)]
  }else{
    ListTreesToAnalyze <- ListTreesToAnalyze[!grepl("500", ListTreesToAnalyze)]
  }
  
  transitions_mean_list <- list()
  
  # load and extract info from each tree
  for (tree_file in ListTreesToAnalyze) {
    tree_name <- gsub("_select100.trees", "", tree_file)
    localTreesDirectory <- paste0(MainDirectory, Folder, "/dta/", tree_name, "_ext")
    dir.create(localTreesDirectory, showWarnings = FALSE)
    
    trees <- readAnnotatedNexus(paste0(MainDirectory, Folder, "/dta/", tree_file))
    if (length(trees) == 0) stop("No trees were loaded from ", tree_file)
    
    # extract unique regions 
    regions <- unique(unlist(lapply(trees[[1]]$annotations, function(x) x$regions.states)))
    regions <- sort(regions)
    
    ## initialize transition matrices
    transitions = matrix(0, nrow = length(regions), ncol = length(regions))
    row.names(transitions) = regions; colnames(transitions) = regions
    transitions_95HPDhigh <- transitions_95HPDlow <- transitions
    
    ## process trees
    tabs = list()
    for (k in 1:length(trees)) {
      cat("Processing tree:", k, "\n")
      tree = trees[[k]]
      tab = DTA_tree_extraction1(tree, mostRecentSamplingDatum = 0)
      write.csv(tab, paste0(localTreesDirectory,"/TreeExtractions_",k,".csv"), row.names=F, quote=F)
      tabs[[k]] = tab
    }
    
    ## compute transition counts
    for (tab in tabs) {
      transition_counts <- table(tab$startLoc, tab$endLoc)
      transitions[rownames(transition_counts), colnames(transition_counts)] <- 
        transitions[rownames(transition_counts), colnames(transition_counts)] + transition_counts
    }
    
    # Compute mean transitions
    transitions_mean <- transitions / length(trees)
    transitions_mean_df <- as.data.frame(as.table(transitions_mean))
    colnames(transitions_mean_df) <- c("from", "to", "mean")
    transitions_mean_df$transition <- paste0(transitions_mean_df$from, ".", transitions_mean_df$to)
    
    # Compute 95% HPD intervals
    
    for (i in 1:dim(transitions_mean)[1]){
      for (j in 1:dim(transitions_mean)[2]){
        if (i != j){
          vS = rep(NA, length(tabs))
          for (h in 1:length(tabs)){
            vS[h] = length(which((tabs[[h]]["startLoc"]==regions[i])&(tabs[[h]]["endLoc"]==regions[j])))
          }
          transitions_95HPDhigh[i,j] = HDInterval::hdi(vS)[2]
          transitions_95HPDlow[i,j] = HDInterval::hdi(vS)[1]
        }
      }
    }
    
    # Convert HPD matrices to dataframes and merge
    transitions_95HPDhigh_df <- as.data.frame(as.table(transitions_95HPDhigh))
    transitions_95HPDlow_df <- as.data.frame(as.table(transitions_95HPDlow))
    colnames(transitions_95HPDhigh_df) <- c("from", "to", "HPDhigh")
    colnames(transitions_95HPDlow_df) <- c("from", "to", "HPDlow")
    transitions_95HPDhigh_df$transition <- paste0(transitions_95HPDhigh_df$from, ".", transitions_95HPDhigh_df$to)
    transitions_95HPDlow_df$transition <- paste0(transitions_95HPDlow_df$from, ".", transitions_95HPDlow_df$to)
    
    transitions_mean_df$HPDhigh <- transitions_95HPDhigh_df$HPDhigh[match(transitions_mean_df$transition, transitions_95HPDhigh_df$transition)]
    transitions_mean_df$HPDlow <- transitions_95HPDlow_df$HPDlow[match(transitions_mean_df$transition, transitions_95HPDlow_df$transition)]
    
    # filter out self-transitions
    transitions_mean_df <- transitions_mean_df[, c("from", "to", "transition", "mean", "HPDlow","HPDhigh")]
    transitions_mean_df <- transitions_mean_df[which(transitions_mean_df$from != transitions_mean_df$to),] 
    write_tsv(transitions_mean_df, paste0(localTreesDirectory, "/transitions.events.tsv"))

    transitions_mean_list[[tree_name]] <- transitions_mean_df
    print(paste0("Done with file: ", tree_name))
  }
  return(transitions_mean_list)
}
