
#' Create visualization of networks in plotly
#' 
#' @param input list of hgnc IDs
#' @importFrom jsonlite toJSON
#' @importFrom httr POST, add_headers
#' @importFrom r2r hashmap, query
#' 
#' @export
#' 
#' @NoRd
#'
#'               
visualizeNetworks = function(input) {
  # Fetch data from INDRA
  input = c("1925", "2092", "3621", "25528", "10289")
  # input = c("1103", "1104", "13575", "25528", "10289")
  uniprot_ids = c("CHK1_HUMAN", "CLH1_HUMAN", "FCGRN_HUMAN", "NECP2_HUMAN", "RFA1_HUMAN")
  gene_id_map = hashmap()
  gene_id_map[c("1925", "2092", "3621", "25528", "10289")] = 
      c("CHEK1", "CLTC", "FCGRT", "NECAP2", "RPA1")
  url = "https://discovery.indra.bio/api/indra_subnetwork_relations"
  hgnc_ids = lapply(input, function(x) list("HGNC", x))
  groundings = list(nodes = hgnc_ids)
  json_body = jsonlite::toJSON(groundings, auto_unbox = TRUE)
  res = POST(url, body = json_body, add_headers("Content-Type" = "application/json"), encode = "raw")
  output = content(res)
  output = Filter(function(x) x$data$stmt_type %in% c("Complex", "Activation", "Inhibition"), output)
  # NEXT STEP: COLLAPSE EDGES WITH SAME SOURCE, TARGET, and INTERACTION AND SUM EVIDENCE COUNT - HELPER FUNCTION
  evidenceList = sapply(output, function(x) 
      paste("https://db.indra.bio/statements/from_agents?subject=", 
            query(gene_id_map, x$source_id), "&object=", 
            query(gene_id_map, x$target_id), "&type=",
            x$data$stmt_type, "&format=html", sep=""))
  
  
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
                      evidenceLink=evidenceList,
                      stringsAsFactors=FALSE)
  
  
  network_id = createNetworkFromDataFrames(nodes,edges, title="my first network", collection="DataFrame Example")
  arrowShapes = mapVisualProperty(
      'Edge Target Arrow Shape','interaction',
      'd',
      c("Complex", "Activation", "Inhibition"),
      c("Arrow","Arrow","Arrow")
  )
  nodeLabels = mapVisualProperty('Node Label','uniprot_id','p')
  edgeWidth = mapVisualProperty('Edge Width','evidenceCount','p')
  createVisualStyle("Y",
                    list(
                        NODE_FILL_COLOR="lightblue",
                        NODE_SHAPE="ROUNDRECT",
                        NODE_SIZE=50,
                        NODE_LABEL_FONT_SIZE=6,
                        NODE_LABEL_POSITION="center"), 
                    list(nodeLabels, arrowShapes, edgeWidth))
  setVisualStyle("Y")
  
}