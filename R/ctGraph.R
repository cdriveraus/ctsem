ctGraphPlot <- function(x,DRIFT=TRUE, DIFFUSION=TRUE){
  if('ctStanModel' %in% class(x)){
    model <- x
    parmats <- listOfMatrices(x$pars)
    A <- parmats$DRIFT
    G <- parmats$DIFFUSION
    
    A <- matrix(as.integer(!A %in% '0'),nrow(A),ncol(A))
    G <- matrix(as.numeric(!G %in% '0'),nrow(G),ncol(G))
  }
  if('ctStanFit' %in% class(x)){
    model <-x$ctstanmodelbase
    A <- x$stanfit$transformedparsfull$DRIFT[1,,]
    G <- x$stanfit$transformedparsfull$DIFFUSION[1,,]
  }
  
  dimnames(A) <- list(model$latentNames,model$latentNames)
  dimnames(G) <- list(model$latentNames,model$latentNames)
  
  
  # Load required libraries
  library(tidygraph)
  library(ggraph)
  library(dplyr)
  library(igraph)  # Only needed for matrix handling
  

  
  # Create node data frame from rownames
  nodes_df <- data.frame(name = rownames(G), stringsAsFactors = FALSE)
  
  # Create edge data frame for G edges (undirected; use only lower triangle to avoid duplicates)
  edges_G <- do.call(rbind, lapply(2:nrow(G), function(i) {
    do.call(rbind, lapply(1:(i - 1), function(j) {
      if (G[i, j] != 0)
        data.frame(from = rownames(G)[j],
          to   = rownames(G)[i],
          weight = G[i, j],
          edge_type = "G",
          stringsAsFactors = FALSE)
    }))
  }))
  edges_G <- edges_G[!is.na(edges_G$from), ]
  
  # Create edge data frame for A edges (directed; include all nonzero entries)
  edges_A <- do.call(rbind, lapply(1:nrow(A), function(i) {
    do.call(rbind, lapply(1:ncol(A), function(j) {
      if (A[i, j] != 0)
        data.frame(from = rownames(A)[i],
          to   = rownames(A)[j],
          weight = A[i, j],
          edge_type = "A",
          stringsAsFactors = FALSE)
    }))
  }))
  edges_A <- edges_A[!is.na(edges_A$from), ]
  
  # Combine the two edge data frames and add a column for layout weights (absolute weight)
  edges_all <- bind_rows(edges_G, edges_A) %>% 
    mutate(layout_weight = abs(weight))
  
  # Create a tidygraph object
  graph_obj <- tbl_graph(nodes = nodes_df, edges = edges_all, directed = TRUE)
  
  # Compute the layout using the positive weights (Fruchterman-Reingold)
  # (Alternatively, you can use create_layout() as below if you want a layout object.)
  layout_obj <- create_layout(graph_obj, layout = "fr", weights = layout_weight)
  
  # Plot the network using ggraph.
  # We use the filter aesthetic (with .data) to choose which edges to plot in each layer.
  ggraph(graph_obj, layout = "fr", weights = layout_weight) +
    # Plot undirected G edges as dashed green arcs (no arrowheads)
    geom_edge_arc(aes(filter = .data$edge_type == "G"),
      edge_colour = "green",
      linetype = "dashed",
      edge_width = 0.8,
      edge_alpha = 0.8) +
    # Plot directed A edges using geom_edge_parallel to offset overlapping edges.
    geom_edge_parallel(aes(filter = .data$edge_type == "A",
      edge_colour = ifelse(.data$weight > 0, "red", "blue")),
      arrow = arrow(length = unit(3, 'mm'), type = "closed"),
      end_cap = circle(3, 'mm'),
      edge_width = 0.8) +
    # Plot nodes and labels.
    geom_node_point(size = 5, color = "black") +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    theme_void() +
    guides(
      edge_alpha = "none",
      edge_width = "none",
      edge_colour = "none"
    )
  
}
