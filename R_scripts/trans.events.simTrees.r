scripts_dir = "~/Dropbox/seraphim/"
source(paste0(scripts_dir,"DTA_tree_extraction_regions.r")) 

trans.events.simTrees <- function(MainDirectory, Folder, scheme){
  transitions.simTrees <- list()
  print(paste0("Working on folder:", " ", Folder))
  
  SimTreesToAnalyze <- list.files(path = paste0(
    MainDirectory,Folder,"/files/"), pattern = "*_modi.nex")
  
  if (scheme == 500) {
    SimTreesToAnalyze <- SimTreesToAnalyze[!grepl("150", SimTreesToAnalyze)]
  }else{
    SimTreesToAnalyze <- SimTreesToAnalyze[!grepl("500", SimTreesToAnalyze)]
  }
  
  for (i in 1:length(SimTreesToAnalyze)){
    tree_file <- SimTreesToAnalyze[i]
    simTree <- readAnnotatedNexus(paste0(MainDirectory, Folder, "/files/", tree_file))
    tree_name = gsub("*_modi.nex", "", tree_file)
    
    # extract unique regions 
    sim.regions <- unique(unlist(lapply(simTree$annotations, function(x) x$regions.states)))
    sim.regions <- sort(sim.regions)
    
    ## initialize transition matrices
    simTransitions = matrix(0, nrow=length(sim.regions), ncol=length(sim.regions))
    row.names(simTransitions) = sim.regions; colnames(simTransitions) = sim.regions
    
    #extract and compute transition counts
    tab = DTA_tree_extraction1(simTree, mostRecentSamplingDatum=0)
    #write.csv(tab, paste0(MainDirectory, Folder, "/files/", tree_name, "_ext.csv"))
    
    transition_counts <- table(tab$startLoc, tab$endLoc)
    simTransitions[rownames(transition_counts), colnames(transition_counts)] <- 
      simTransitions[rownames(transition_counts), colnames(transition_counts)] + transition_counts
    
    #convert matrices in dataframes and re-format
    simTransitions.df <- as.data.frame(as.table(simTransitions))
    colnames(simTransitions.df) <- c("from", "to", "count")
    simTransitions.df$transition <- paste0(simTransitions.df$from, ".", simTransitions.df$to)
    simTransitions.df <- simTransitions.df[which(simTransitions.df$from != simTransitions.df$to),]
    simTransitions.df <- simTransitions.df[which(simTransitions.df$count != 0),]
    #write_tsv(simTransitions.df, paste0(MainDirectory, Folder, "/files/", tree_name, "_trans.events.tsv"))
    
    transitions.simTrees[[tree_name]] <- simTransitions.df
    
  }
  
  return(transitions.simTrees)
  
}