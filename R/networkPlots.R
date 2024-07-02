
#' Fetch data on a set of proteins from INDRA
#' 
#' @param input groupComparison comparisionResult table
#' @param pvalue_cutoff p-value cutoff for filtering
#' @param stmt_types types of statements to filter, default is "Complex"
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @importFrom r2r hashmap query keys
#' @importFrom dplyr filter
#' 
#' @export
#' 
#' @noRd
#'
#'   
fetchIndraData = function(input, pvalue_cutoff = 0.05, stmt_types = c("Complex")) {
    input = filter(input, adj.pvalue < pvalue_cutoff)
    input = filter(input, is.na(issue))
    hgnc_ids = as.character(input$HgncId)
    uniprot_ids = input$Protein
    gene_id_map = hashmap()
    gene_id_map[input$HgncId] = input$HgncName
    url = "https://discovery.indra.bio/api/indra_subnetwork_relations"
    groundings = lapply(hgnc_ids, function(x) list("HGNC", x))
    groundings = list(nodes = groundings)
    json_body = jsonlite::toJSON(groundings, auto_unbox = TRUE)
    res = POST(url, body = json_body, add_headers("Content-Type" = "application/json"), encode = "raw")
    output = content(res)
    output = Filter(function(x) x$data$stmt_type %in% stmt_types, output)
    
    edge_data = hashmap()
    for (edge in output) {
        key = paste(edge$source_id, edge$target_id, edge$data$stmt_type, sep="_")
        if (key %in% keys(edge_data)) {
            edge_data[[key]]$data$evidence_count = edge_data[[key]]$data$evidence_count + edge$data$evidence_count
        } else {
            edge_data[[key]] = edge
        }
    }
    evidenceList = sapply(keys(edge_data), function(x) 
        paste("https://db.indra.bio/statements/from_agents?subject=", 
              query(gene_id_map, query(edge_data, x)$source_id), "&object=", 
              query(gene_id_map, query(edge_data, x)$target_id), "&type=",
              query(edge_data, x)$data$stmt_type, "&format=html", sep=""))
    
    
    nodes = data.frame(id=hgnc_ids,
                       uniprot_id=uniprot_ids,
                       logFC=input$log2FC, 
                       pvalue=input$adj.pvalue,
                       stringsAsFactors=FALSE)
    edges = data.frame(source=sapply(keys(edge_data), function(x) query(edge_data, x)$source_id),
                       target=sapply(keys(edge_data), function(x) query(edge_data, x)$target_id),
                       interaction=sapply(keys(edge_data), function(x) query(edge_data, x)$data$stmt_type), 
                       evidenceCount=sapply(keys(edge_data), function(x) query(edge_data, x)$data$evidence_count),
                       evidenceLink=evidenceList,
                       stringsAsFactors=FALSE)
    
    return(list(nodes=nodes, edges=edges))
}

#' Create visualization of networks in cytoscape
#' 
#' @param nodes dataframe of nodes
#' @param edges dataframe of edges
#' @importFrom RCy3 createNetworkFromDataFrames mapVisualProperty createVisualStyle setVisualStyle
#' 
#' @export
#' 
#' @noRd
#'
#'               
visualizeNetworks = function(nodes, edges) {
  network_id = createNetworkFromDataFrames(nodes, edges, title="my first network", collection="DataFrame Example")
  arrowShapes = mapVisualProperty(
      'Edge Target Arrow Shape','interaction',
      'd',
      c("Complex", "Activation", "Inhibition"),
      c("Arrow","Arrow","Arrow")
  )
  nodeLabels = mapVisualProperty('Node Label','uniprot_id','p')
  createVisualStyle("Y",
                    list(
                        NODE_FILL_COLOR="lightblue",
                        NODE_SHAPE="ROUNDRECT",
                        NODE_SIZE=50,
                        NODE_LABEL_FONT_SIZE=6,
                        NODE_LABEL_POSITION="center"), 
                    list(
                        nodeLabels, 
                        arrowShapes
                    ))
  setVisualStyle("Y")
  
}