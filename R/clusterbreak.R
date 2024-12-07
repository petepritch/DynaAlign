.onLoad <- function(libname, pkgname) {
  bioc_packages <- c("Biostrings", "DECIPHER")
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing missing Bioconductor package: %s", pkg))
      BiocManager::install(pkg, ask = FALSE)
    }
  }
}


#' Modified Louvain Clustering function
#'
#' @param gin Input network.
#' @param res Louvain algorithm sensitivity.
#' @param res_range Percent of sensitivity for the function loop over for best modularity.
#' @param res_step Number of even breaks over the sensitivity range.
#' @param itr Number of iteration over each break point. 
#' @return A clustered graph network optimized for modularity
#' @importFrom igraph graph_from_adjacency_matrix E cluster_louvain V layout_with_fr
#' @export
louvain_mod<-function(gin,res,res_range_perc=0,res_step=0,itr=3){
  
  res=seq(res - res_range_perc*res, res + res_range_perc*res, by = res_step)# set range of resolution parameters
  
  best_modularity <- NULL
  best_clusters <- NULL
  best_resolution<-NULL
  for (j in seq_along(res)){
    if (j==1){
      best_resolution<-res[j]
    }
    for (i in 1:itr) {#iteration loop
      gin.cluster<-igraph::cluster_louvain(gin,weights=igraph::E(gin)$weight,resolution=res[j]) # cluster using weights
      gin.mod<-  igraph::modularity(gin.cluster)#obtain modularity
      
      if(i==1){#initial loop
        best_modularity <- gin.mod
        best_clusters <- gin.cluster
      }
      if (i>1 & gin.mod > best_modularity) {#replace initial values if higher modularity
        best_modularity <- gin.mod
        best_clusters <- gin.cluster
        best_resolution<-res[j]
      }
    }
  }
  
  return(list(cluster=best_clusters,
              resolution=best_resolution,
              modularity=best_modularity))
}

#' Generate cluster id using graph network and louvain method
#' 
#' @param pepmat A square or upper triangular similarity matrix
#' @param igraph_mode Mode settings for igraph (default: "upper")
#' @param igraph_weight Weight setting (TRUE/NULL) for igraph (default: TRUE)
#' @param louvain_resolution Resolution setting for Louvain algorithm (default: 1.05)
#' @importFrom igraph graph_from_adjacency_matrix E cluster_louvain
#' @return A vector containing cluster assignment
#' @export
netcluster<-function(pepmat,igraph_mode="upper",igraph_weight=TRUE,louvain_resolution=1.05,louvain_range_perc=0,louvain_step=0,louvain_itr=3){
  if(nrow(pepmat)!=ncol(pepmat)) {
    stop("Input must be a square pairwise similarity matrix")
  }
  network<-igraph::graph_from_adjacency_matrix(pepmat,mode=igraph_mode,weighted=igraph_weight)
  # cluster.n<-igraph::cluster_louvain(network,weights=igraph::E(network)$weight,resolution=louvain_resolution)
  cluster.n<-louvain_mod(network,res=louvain_resolution,res_range_perc=louvain_range_perc,res_step=louvain_step,itr=louvain_itr)
  return(cluster.n$cluster$membership)
}

#' Generate clusters with specified sizes using graph network and louvain method
#' 
#' @param pep A vector of peptide sequences
#' @param thresh Threshold similarity score used to remove edge between two sequences if not similar enough
#' @param k_size k-mer size for MinHash algorithm (default: k_size = 2)
#' @param hash_size Hash function size for MinHash algorithm (default: hash_size = 50)
#' @param size_max Maximum size of cluster desired (default: size_max = 10)
#' @param size_min Minimum size of cluster desired (default: size_min = 3)
#' @param sens Resolution setting for Louvain algorithm (default: sens = 1.05)
#' @param max_itr Maximum function calls wanted before halting function execution (default: max_itr = 500)
#' @importFrom igraph graph_from_adjacency_matrix E cluster_louvain
#' @return A nx2 matrix with a column containing n peptide sequences and their corresponding cluster assignment
#' @export
#' 
clusterbreak <- function(pep, thresh=0.8,k_size=2, hash_size=50, size_max=10, size_min=3, sens=1.05, max_itr=10000) {
  if (size_max <= size_min) {
    stop("size_max must be greater than size_min")
  }
  if (length(pep) == 0) {
    stop("empty input sequence vector")
  }
  
  # Create new state environment
  state <- new.env()
  state$out.df <- matrix(nrow = 0, ncol = 2)
  state$itr <- 1
  state$clusters_processed <- 0
  state$convergence <- 1
  
  cluster_recursive <- function(pep) {
    
    # Create helper function for logging adapted from Claude AI prompt
    log_message <- function(msg, level="INFO") {
      timestamp <- format(Sys.time(), "%H:%M:%S")
      cat(sprintf("[%s] %s: %s\n", timestamp, level, msg))
    }
    
    if (state$itr > max_itr) {
      log_message("Maximum function calls reached", "WARNING")
      state$convergence<-0 # change status to not converged since max itr reached
      return(state$out.df)
      break
    }
    
    pep.sim <- minhash(pep, k_size, hash_size)  #minhash similarity matrix
    pep.sim <- pep.sim$dist_matrix
    pep.sim[pep.sim<thresh] <- 0 # remove edges from nodes with similarity below threshold
    c.index <- netcluster(pep.sim, louvain_resolution=sens) #cluster id
    pep.ref <- cbind(pep, c.index) # combine cluster id with sequences
    c.size <- tabulate(c.index) # count each cluster size
    id.itr <- which(c.size > size_max) # cluster id above max size
    id.rm <- which(c.size < size_min) # cluster id below min size
    
    # stop or recursion conditions for output
    if (length(id.itr) == 0) {  
      out.pep <- pep.ref[!pep.ref[,2] %in% id.rm,] # filter out those under minimum length
      if (nrow(out.pep) > 0) { #ensure non-empty output
        out.pep[,2] <- paste0(state$itr, ".", out.pep[,2]) #combine iteration and cluster to ensure uniqueness
        state$out.df <- rbind(state$out.df, out.pep) # combine output df
        state$clusters_processed <- state$clusters_processed + nrow(out.pep) # sum number of clusters processed
      }
    } else {
      pep.out <- pep.ref[(!pep.ref[,2] %in% id.rm) & (!pep.ref[,2] %in% id.itr),] # recursion condition
      if (nrow(pep.out) > 0) {
        pep.out[,2] <- paste0(state$itr, ".", pep.out[,2])#combine iteration and cluster to ensure uniqueness
        state$out.df <- rbind(state$out.df, pep.out)# combine output df
        state$clusters_processed <- state$clusters_processed + nrow(pep.out)# sum number of clusters processed
      }
      
      # filter remaining clusters
      pep.new <- pep.ref[pep.ref[,2] %in% id.itr,]#clusters with size > max
      loop.lvl <- unique(pep.new[,2])#set levels for loop
      
      # loop for remaining qualified clusters
      for (i in seq_along(loop.lvl)) {
        sub.df <- pep.new[pep.new[,2] == loop.lvl[i], 1] #filter ith cluster
        state$itr <- state$itr + 1 #increment for iteration/function call
        cluster_recursive(sub.df) #recursive call
      }
    }
    
    return(state$out.df)
  }
  
  # run recursive clustering
  result <- cluster_recursive(pep)
  
  # Final status report adapted from claude AI output
  if (state$convergence==1){
    cat(sprintf("\nClustering complete:\n"))
  }else{
    cat(sprintf("\nClustering incomplete, consider adjusting parameters:\n"))
  }
  cat(sprintf("Total clusters processed: %d\n", state$clusters_processed))
  cat(sprintf("Total function calls: %d\n", state$itr))
  
  rm(state) # reset state global 
  
  return(result)
}

#' Generate consensus sequence
#'
#' @param df Input clusterbreak output df
#' @return df with first column being unique cluster id and second column being the corresponding consensus sequence
#' @importFrom Biostrings AAStringSet
#' @importFrom DECIPHER AlignSeqs ConsensusSequence
#' @export
clusterconsensus<-function(df){
  cluster.id<-unique(df[,2])
  out.df<-matrix(nrow=0,ncol=2)
  for (i in 1:length(cluster.id)){
    df.sub<-df[which(df[,2] == cluster.id[i]),1]
    aa.set<-Biostrings::AAStringSet(df.sub)
    aa.set.align<-DECIPHER::AlignSeqs(aa.set)
    con.seq<-as.character(DECIPHER::ConsensusSequence(aa.set.align))
    out.df<-rbind(out.df,c(cluster.id[i],con.seq))
  }
  return(out.df)
}

#' Plot consensus sequences for each cluster in a clustered network
#'
#' @param df Input clusterconsensus function output
#' @param k_size minhash kmer size
#' @param threshold binary threshold for adjacency matrix
#' @param sens Louvain algorithm sensitivity
#' @return hash_size hash function size
#' @importFrom igraph graph_from_adjacency_matrix E cluster_louvain V layout_with_fr
#' @export
consensusplot<-function(df,k_size=2, hash_size=50,threshold=0.8,sens=1.05,...){
  #similarity matrix and adjacency matrix
  df.hash<-minhash(df[,2],k_size,hash_size)
  df.hash<-df.hash$dist_matrix
  df.hash[df.hash<threshold]<-0
  #plot call
  g<-igraph::graph_from_adjacency_matrix(df.hash,mode="upper",weighted=TRUE) # base plot
  g.weight<-igraph::E(g)$weight # edge weight
  g.cluster<-igraph::cluster_louvain(g,weights=g.weight,resolution=sens) # cluster using weights
  
  g.layout<-igraph::layout_with_fr(g, weights = g.weight) #set layout for plot
  igraph::V(g)$name <- df[,1]#set node names to cluster name
  g.out<-plot(g.cluster,g,layout = g.layout, ...)
}