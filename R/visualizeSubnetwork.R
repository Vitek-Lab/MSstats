#' Create visualization of subnetwork in cytoscape
#' 
#' @param nodes dataframe of nodes
#' @param edges dataframe of edges
#' @param pvalue_cutoff p-value cutoff for highlighting significant proteins. Default is 0.05
#' @param logfc_cutoff log fold change cutoff for highlighting significant proteins. Default is 0.5
#' @importFrom RCy3 createNetworkFromDataFrames mapVisualProperty createVisualStyle setVisualStyle setNodeColorMapping
#' 
#' @export
#' 
#' @examples 
#' input = data.table::fread(system.file("tinytest/processed_data/groupComparisonModel.csv", 
#'                              package = "MSstats"))
#' # subnetwork = getSubnetworkFromIndra(input)
#' # visualizeSubnetwork(subnetwork$nodes, subnetwork$edges)
#'
#'               
visualizeSubnetwork = function(nodes, edges, pvalue_cutoff = 0.05, logfc_cutoff = 0.5) {
    DEFAULT_VISUAL_STYLE = list(
        # NODE_FILL_COLOR="lightblue",
        NODE_SHAPE="ROUNDRECT",
        NODE_SIZE=50,
        NODE_LABEL_FONT_SIZE=6,
        NODE_LABEL_POSITION="center",
        EDGE_TARGET_ARROW_SHAPE="Arrow")
    VISUAL_STYLE_NAME = "MSstats-Indra Visual Style"
    
    createNetworkFromDataFrames(nodes, edges)
    control.points <- c (-logfc_cutoff, 0.0, logfc_cutoff)
    colors <-  c ('#CC5522',  '#DD8855', '#FFFFFF','#5588DD','#2255CC')
    RCy3::setNodeColorMapping("logFC", control.points, colors)
    
    VISUAL_STYLE_MAPPINGS = list(
        mapVisualProperty('Node Label','id','p')
    )
    createVisualStyle(VISUAL_STYLE_NAME, DEFAULT_VISUAL_STYLE, VISUAL_STYLE_MAPPINGS)
    setVisualStyle(VISUAL_STYLE_NAME)
}