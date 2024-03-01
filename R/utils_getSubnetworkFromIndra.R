#' Call INDRA Cogex API and return response
#' @param hgnc_ids list of hgnc_ids
#' @return list of INDRA statements
#' @importFrom jsonlite toJSON
#' @importFrom httr POST add_headers content
#' @keywords internal
#' @noRd
.callIndraCogexApi = function(hgnc_ids) {
    INDRA_COGEX_URL = "https://discovery.indra.bio/api/indra_subnetwork_relations"
    
    groundings = lapply(hgnc_ids, function(x) list("HGNC", x))
    groundings = list(nodes = groundings)
    groundings = jsonlite::toJSON(groundings, auto_unbox = TRUE)

    res = POST(
        INDRA_COGEX_URL, 
        body = groundings, 
        add_headers("Content-Type" = "application/json"), 
        encode = "raw"
    )
    res = content(res)
    return(res)
}

#' Filter groupComparison result input based on user-defined cutoffs
#' @param input groupComparison result
#' @param pvalue_cutoff p-value cutoff
#' @return filtered groupComparison result
#' @keywords internal
#' @noRd
.filterGetSubnetworkFromIndraInput = function(input, pvalue_cutoff) {
    if (!is.null(pvalue_cutoff)) {
        input = input[input$adj.pvalue < pvalue_cutoff,]
    }
    input = input[is.na(input$issue),]
    return(input)
}

#' Add additional metadata to an edge
#' @param edge object representation of an INDRA statement
#' @param input filtered groupComparison result
#' @return edge with additional metadata
#' @keywords internal
#' @noRd
.addAdditionalMetadataToIndraEdge = function(edge, input) {
    edge$evidence_list = paste(
        "https://db.indra.bio/statements/from_agents?subject=", 
        edge$source_id, "@HGNC&object=", 
        edge$target_id, "@HGNC&type=",
        edge$data$stmt_type, "&format=html", sep="")
    edge$source_uniprot_id = input[input$HgncId == edge$source_id,]$Protein
    edge$target_uniprot_id = input[input$HgncId == edge$target_id,]$Protein
    return(edge)
}


#' Collapse duplicate INDRA statements into a mapping of edge to metadata
#' @param res INDRA response
#' @param input filtered groupComparison result
#' @importFrom r2r hashmap keys
#' @return processed edge to metadata mapping
#' @keywords internal
#' @noRd
.collapseDuplicateEdgesIntoEdgeToMetadataMapping = function(res, input) {
    edge_to_metadata_mapping = hashmap()

    for (edge in res) {
        key = paste(edge$source_id, edge$target_id, edge$data$stmt_type, sep="_")
        if (key %in% keys(edge_to_metadata_mapping)) {
            edge_to_metadata_mapping[[key]]$data$evidence_count = 
                edge_to_metadata_mapping[[key]]$data$evidence_count + 
                    edge$data$evidence_count
        } else {
            edge = .addAdditionalMetadataToIndraEdge(edge, input)
            edge_to_metadata_mapping[[key]] = edge
        }
    }
    
    return(edge_to_metadata_mapping)
}

#' Construct edges data.frame from INDRA response
#' @param res INDRA response
#' @param input filtered groupComparison result
#' @importFrom r2r query keys
#' @return edge data.frame
#' @keywords internal
#' @noRd
.constructEdgesDataFrame = function(res, input) {
    res = .collapseDuplicateEdgesIntoEdgeToMetadataMapping(res, input)
    edges = data.frame(source=sapply(keys(res), function(x) query(res, x)$source_uniprot_id),
                       target=sapply(keys(res), function(x) query(res, x)$target_uniprot_id),
                       interaction=sapply(keys(res), function(x) query(res, x)$data$stmt_type), 
                       evidenceCount=sapply(keys(res), function(x) query(res, x)$data$evidence_count),
                       evidenceLink=sapply(keys(res), function(x) query(res, x)$evidence_list),
                       stringsAsFactors=FALSE)
    return(edges)
}

#' Construct nodes data.frame from groupComparison output
#' @param input filtered groupComparison result
#' @return nodes data.frame
#' @keywords internal
#' @noRd
.constructNodesDataFrame = function(input) {
    nodes = data.frame(id=input$Protein,
                       logFC=input$log2FC, 
                       pvalue=input$adj.pvalue,
                       stringsAsFactors=FALSE)
    return(nodes)
}
