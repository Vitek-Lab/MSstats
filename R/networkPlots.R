
#' Create visualization of networks in plotly
#' 
#' @param input list of hgnc IDs
#' @importFrom jsonlite toJSON
#' @importFrom httr POST, add_headers
#' 
#' @export
#' 
#' @NoRd
#'
#'               
visualizeNetworks = function(input) {
  # Fetch data from INDRA
  input = c("1925", "2092", "3621", "25528", "10289")
  uniprot_ids = c("CHK1_HUMAN", "CLH1_HUMAN", "FCGRN_HUMAN", "NECP2_HUMAN", "RFA1_HUMAN")
  url = "https://discovery.indra.bio/api/indra_subnetwork_relations"
  hgnc_ids = lapply(input, function(x) list("HGNC", x))
  groundings = list(nodes = hgnc_ids)
  json_body = jsonlite::toJSON(groundings, auto_unbox = TRUE)
  res = POST(url, body = json_body, add_headers("Content-Type" = "application/json"), encode = "raw")
  output = content(res)
  output = Filter(function(x) x$data$stmt_type %in% c("Complex", "Activation", "Inhibition"), output)
  
  
  # Construct Cytoscape network
  nodes = data.frame(id=input,
                     uniprot_id=uniprot_ids,
                      logFC=c(0.1,-0.3,-0.7,0.7,-0.6), 
                      pvalue=c(0.99, 0.003, 0.001, 0.00000002, 0.0000000003),
                      stringsAsFactors=FALSE)
  edges = data.frame(source=sapply(output, function(x) x$source_id),
                      target=sapply(output, function(x) x$target_id),
                      interaction=sapply(output, function(x) x$data$stmt_type), 
                      correlation=c(0.2,0.3,0.4,0.4),
                      evidenceCount=sapply(output, function(x) x$data$evidence_count),
                      evidenceList=sapply(output, function(x) jsonlite::toJSON(jsonlite::fromJSON(x$data$stmt_json)$evidence)),
                      stringsAsFactors=FALSE)
  
  
  createNetworkFromDataFrames(nodes,edges, title="my first network", collection="DataFrame Example")
  arrowShapes = mapVisualProperty(
      'Edge Target Arrow Shape','interaction',
      'd',
      c("Complex", "Activation", "Inhibition"),
      c("Arrow","Arrow","Arrow")
  )
  nodeLabels = mapVisualProperty('Node Label','uniprot_id','p')
  edgeWidth = mapVisualProperty('Edge Width','evidenceCount','p')
  createVisualStyle("Z3",
                    list(
                        NODE_SHAPE="ROUNDRECT",
                        NODE_WIDTH=50,
                        NODE_HEIGHT=30,
                        NODE_LABEL_FONT_SIZE=5,
                        NODE_LABEL_POSITION="center"), 
                    list(nodeLabels, arrowShapes, edgeWidth))
  setVisualStyle("Z3")
  
}