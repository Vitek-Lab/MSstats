
#' Create visualization of networks in plotly
#' 
#' @param input list of Uniprot protein IDs
#' @importFrom hgnc import_hgnc_dataset latest_archive_url filter_by_keyword
#' @importFrom jsonlite toJSON
#' @importFrom httr POST, add_headers
#' 
#' @export
#' 
#' @NoRd
#'
#'               
visualizeNetworks = function(input) {
  # library(igraph)
  # library(plotly)
  input = c("O14558", "P04792", "P23141", "P01833", "Q6E0U4", "P60985", "P17661")
  url = "https://discovery.indra.bio/api/indra_subnetwork_relations"
  dt = import_hgnc_dataset(file = latest_archive_url())
  rows = lapply(input, function(x) filter_by_keyword(dt, x, cols = c("uniprot_ids")))
  hgnc_ids = lapply(rows, function(x) list("HGNC", x$hgnc_id2))
  groundings = list(nodes = hgnc_ids)
  json_body = jsonlite::toJSON(groundings, auto_unbox = TRUE)
  res = POST(url, body = json_body, add_headers("Content-Type" = "application/json"), encode = "raw")
  output = content(res)
  hgnc_ids_2 = lapply(rows, function(x) x$hgnc_id2)
  vertices = data.frame(node = unlist(hgnc_ids_2), protein_id = input)
  edges = data.frame(
      from = lapply(output, function(x) x$source_id), 
      to = lapply(output, function(x) x$target_id),
      evidence_count = lapply(output, function(x) x$data$evidence_count),
      belief = lapply(output, function(x) x$data$belief),
      stmt_type = lapply(output, function(x) x$data$stmt_type)
  )
  g <- graph_from_data_frame(edges, directed=TRUE, vertices=vertices)
  L = layout.circle(g)
  vs <- V(G)
  es <- as.data.frame(get.edgelist(G))
  
  Nv <- length(vs)
  Ne <- length(es[1]$V1)
  Xn <- L[,1]
  Yn <- L[,2]
  
  network <- plot_ly(x = ~Xn, y = ~Yn, mode = "markers", text = vs$label, hoverinfo = "text")
  
  edge_shapes <- list()
  for(i in 1:Ne) {
      v0 <- es[i,]$V1
      v1 <- es[i,]$V2
      
      edge_shape = list(
          type = "line",
          line = list(color = "#030303", width = 0.3),
          x0 = Xn[v0],
          y0 = Yn[v0],
          x1 = Xn[v1],
          y1 = Yn[v1]
      )
      
      edge_shapes[[i]] <- edge_shape
  }
  
  axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  
  fig <- layout(
      network,
      title = 'INDRA Network',
      shapes = edge_shapes,
      xaxis = axis,
      yaxis = axis
  )
  
  fig
  
  # Create a graph
  # g <- graph_from_data_frame(input, directed = FALSE)
  # 
  # # Plot the graph
  # p <- plot_ly(data = input, x = ~x, y = ~y, text = ~name, mode = "markers+text", textposition = "bottom center") %>%
  #   add_trace(data = input, x = ~x, y = ~y, xend = ~xend, yend = ~yend, mode = "lines") %>%
  #   layout(title = "Network Visualization", showlegend = FALSE)
  # 
  # p
}