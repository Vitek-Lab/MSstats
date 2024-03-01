#' Get subnetwork from INDRA database with differential analysis results
#' 
#' @param input groupComparison comparisionResult table with additional HGNC ID 
#' and HGNC name columns
#' @param pvalue_cutoff p-value cutoff for filtering. Default is NULL, i.e. no 
#' filtering
#' 
#' @return list of 2 data.frames, nodes and edges
#' 
#' @export
#' 
#' @examples 
#' input = data.table::fread(system.file("tinytest/processed_data/groupComparisonModel.csv", 
#'                              package = "MSstats"))
#' # subnetwork = getSubnetworkFromIndra(input, pvalue_cutoff = 0.05)
#' # head(subnetwork$nodes)
#' # head(subnetwork$edges)
#'
#'   
getSubnetworkFromIndra = function(input, pvalue_cutoff = NULL) {
    input = .filterGetSubnetworkFromIndraInput(input, pvalue_cutoff)
    res = .callIndraCogexApi(input$HgncId)
    nodes = .constructNodesDataFrame(input)
    edges = .constructEdgesDataFrame(res, input)
    return(list(nodes=nodes, edges=edges))
}