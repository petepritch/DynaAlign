set.seed(123)
data <- matrix(runif(100), nrow = 10, ncol = 10)

##' Cosine similarity function for generating similarity matrices
##'
##' @title
##' @param x
##' @return
##' @author Pete Pritchard
cosine_similarity <- function(x) {

  x_norm <- sqrt(rowSums(x^2))
  x_normalized <- x / x_norm

  similarity_matrix <- x_normalized %*% t(x_normalized)

  return(similarity_matrix)
}

similarity_matrix <- cosine_similarity(data)

print(similarity_matrix)



##' Plotting...
##'
##' @title
##' @param J
##' @param thres
##' @return
##' @author Pete Pritchard
plotGraph <- function(J, thres) {

  J[J < thres] <- 0

  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

  clusters <- cluster_fast_greedy(g)

  return( plot(clusters, g,
               vertex.size = 15,
               vertex.label = NA,
               edge.width = E(g)$weight * 2,
               main = "Test graph"
  ) )

}
