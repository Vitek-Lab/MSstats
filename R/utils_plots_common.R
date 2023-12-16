#' Theme for MSstats plots
#' 
#' @param type type of a plot
#' @param x.axis.size size of text on the x axis
#' @param y.axis.size size of text on the y axis
#' @param legend_size size of the legend
#' @param strip_background background of facet
#' @param strip_text_x size of text on facets
#' @param legend_position position of the legend
#' @param legend_box legend.box
#' @param text_angle angle of text on the x axis (for condition and comparison plots)
#' @param text_hjust hjust parameter for x axis text (for condition and comparison plots)
#' @param text_vjust vjust parameter for x axis text (for condition and comparison plots)
#' @param ... additional parameters passed on to ggplot2::theme()
#' 
#' @import ggplot2
#' @export
#' 
theme_msstats = function(
    type, x.axis.size = 10, y.axis.size = 10, legend_size = 13, 
    strip_background = element_rect(fill = "gray95"),
    strip_text_x = element_text(colour = c("black"), size = 14),
    legend_position = "top", legend_box = "vertical", text_angle = 0, text_hjust = NULL, text_vjust = NULL, 
    ...
) {
    if (type %in% c("CONDITIONPLOT", "COMPARISONPLOT")) {
        ggplot2::theme(
            panel.background = element_rect(fill = 'white', colour = "black"),
            axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
            axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
            axis.ticks = element_line(colour = "black"),
            title = element_text(size = x.axis.size + 8, vjust = 1.5),
            panel.grid.major.y = element_line(colour = "grey95"),
            panel.grid.minor.y = element_blank(),
            axis.text.y = element_text(size = y.axis.size, colour = "black"),
            axis.text.x = element_text(size = x.axis.size, colour = "black",
                                       angle = text_angle, hjust = text_hjust, 
                                       vjust = text_vjust),
            ...
        )
    } else {
        ggplot2::theme(
            panel.background = element_rect(fill = 'white', colour = "black"),
            legend.key = element_rect(fill = 'white', colour = 'white'),    
            panel.grid.minor = element_blank(),
            strip.background = strip_background,
            axis.text.x = element_text(size = x.axis.size, colour = "black"),
            axis.text.y = element_text(size = y.axis.size, colour = "black"),
            axis.title.x = element_text(size = x.axis.size + 5, vjust = -0.4),
            axis.title.y = element_text(size = y.axis.size + 5, vjust = 0.3),
            axis.ticks = element_line(colour = "black"),
            title = element_text(size = x.axis.size + 8, vjust = 1.5),
            strip.text.x = strip_text_x,
            legend.position = legend_position,
            legend.box = legend_box,
            legend.text = element_text(size = legend_size),
            ...
        )
    }
}

#' Get proteins based on names or integer IDs
#' 
#' @param chosen_proteins protein names or integers IDs
#' @param all_proteins all unique proteins
#' 
#' @return character
#' 
#' @export
getSelectedProteins = function(chosen_proteins, all_proteins) {
    if (is.character(chosen_proteins)) {
        selected_proteins = chosen_proteins
        missing_proteins = setdiff(selected_proteins, all_proteins)
        if (length(missing_proteins) > 0) {
            stop(paste("Please check protein name. Dataset does not have this protein. -", 
                       toString(missing_proteins), sep = " "))
        }
    }
    if (is.numeric(chosen_proteins)) {
        selected_proteins <- all_proteins[chosen_proteins]
        if (length(all_proteins) < max(chosen_proteins)) {
            stop(paste("Please check your selection of proteins. There are ",
                       length(all_proteins)," proteins in this dataset."))
        }
    }
    selected_proteins
}


#' Save a plot to pdf file
#' 
#' @inheritParams .saveTable
#' @param width width of a plot
#' @param height height of a plot
#' 
#' @return NULL
#' 
#' @export
#' 
savePlot = function(name_base, file_name, width, height) {
    if (name_base != FALSE) {
        file_path = getFileName(name_base, file_name, width, height)
        file_path = paste0(file_path,".pdf")
        pdf(file_path, width = width, height = height)
    }
    NULL
}

getFileName = function(name_base, file_name, width, height) {
    all_files = list.files(".")
    if(file_name == 'ProfilePlot'){
        num_same_name = sum(grepl(paste0("^", name_base, file_name, "_[0-9]?"), all_files))
    } else {
        num_same_name = sum(grepl(paste0("^", name_base, file_name, "[0-9]?"), all_files))
    }
    if (num_same_name > 0) {
        file_name = paste(file_name, num_same_name + 1, sep = "_")
    }
    file_path = paste0(name_base, file_name)
    return(file_path)
}


#' Save a data table to a file
#' @param input data.table
#' @param name_base path to a folder (or "" for working directory)
#' @param file_name name of a file to save. If this file already exists,
#' an integer will be appended to this name
#' @return NULL
#' @keywords internal
.saveTable = function(input, name_base, file_name) {
    if (name_base != FALSE) {
        all_files = list.files(".")
        num_same_name = sum(grepl(paste0("^", file_name, "_[0-9]?"), all_files))
        if (num_same_name > 0) {
            file_name = paste(file_name, num_same_name + 1, sep = "_")
        }
        file_path = paste0(paste0(name_base, "/"), file_name, ".pdf")
        data.table::fwrite(input, file = file_path)
    }
    NULL
}
