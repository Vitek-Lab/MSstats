#' Plot quality metrics from converter output
#'
#' Visualizes a quality metric column from the output of MSstats converter
#' functions against run order for a single protein. Each
#' PeptideSequence + PrecursorCharge combination is drawn as a distinct
#' coloured line, mirroring the feature-level view in
#' \code{\link[MSstats]{dataProcessPlots}}.
#'
#' @param input data.frame or data.table returned by an MSstatsConvert
#'   converter function (e.g. \code{SpectronauttoMSstatsFormat}).
#' @param metric character, name of the column to plot on the y-axis.
#'   Defaults to \code{"AnomalyScores"}. Must be a column of \code{input}.
#' @param which.Protein character, name of the protein to plot. Required.
#' @param address prefix for the filename used when saving the plot.
#'   If \code{FALSE} (default), the plot is returned without saving.
#'   When \code{isPlotly = FALSE} a PDF is saved; when \code{isPlotly = TRUE}
#'   an HTML file is saved.
#' @param isPlotly logical. If \code{TRUE} returns an interactive
#'   \code{\link[plotly]{plotly}} object (and saves as HTML when
#'   \code{address} is provided). If \code{FALSE} (default) returns a
#'   \code{\link[ggplot2]{ggplot}} object.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object, or a \code{plotly} object
#'   when \code{isPlotly = TRUE}.
#'
#' @details
#' The x-axis order is determined by the factor levels of the \code{Run}
#' column. When \code{runOrder} is passed to the converter the \code{Run}
#' column is automatically set to an ordered factor; otherwise the runs appear
#' in alphabetical order.
#'
#' Metric values are averaged across fragment ions within each
#' PeptideSequence + PrecursorCharge + Run combination before plotting, so
#' each precursor contributes exactly one point per run.
#'
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @importFrom htmltools save_html
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- SpectronauttoMSstatsFormat(
#'   input, calculateAnomalyScores = TRUE,
#'   anomalyModelFeatures = c("FGShapeQualityScoreMS2", "EGDeltaRT"),
#'   anomalyModelFeatureTemporal = c("mean_decrease", "dispersion_increase"),
#'   runOrder = my_run_order
#' )
#' MSstatsQualityMetricsPlot(result, which.Protein = "ProteinA")
#' MSstatsQualityMetricsPlot(result, metric = "EGDeltaRT",
#'                           which.Protein = "ProteinA", isPlotly = TRUE)
#' }
MSstatsQualityMetricsPlot <- function(input, metric = "AnomalyScores",
                                      which.Protein,
                                      address = FALSE, isPlotly = FALSE) {
    if (missing(which.Protein)) {
        stop("'which.Protein' is required. Please specify a protein name.")
    }

    input_df <- as.data.frame(input)

    required_cols <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "Run")
    missing_cols <- setdiff(required_cols, colnames(input_df))
    if (length(missing_cols) > 0) {
        stop(paste0(
            "Required column(s) not found in input: ",
            paste(missing_cols, collapse = ", ")
        ))
    }
    if (!metric %in% colnames(input_df)) {
        stop(paste0(
            "Column '", metric, "' not found in input. ",
            "Available columns: ", paste(colnames(input_df), collapse = ", ")
        ))
    }
    if (!which.Protein %in% input_df$ProteinName) {
        stop(paste0("Protein '", which.Protein, "' not found in input."))
    }

    input_df <- input_df[input_df$ProteinName == which.Protein, ]

    if (!is.factor(input_df$Run)) {
        input_df$Run <- factor(input_df$Run)
    }

    input_df$Precursor <- paste(input_df$PeptideSequence,
                                input_df$PrecursorCharge, sep = "_")

    # Average across fragment ions so each precursor has one value per run
    plot_df <- aggregate(
        input_df[[metric]],
        by  = list(Run = input_df$Run, Precursor = input_df$Precursor),
        FUN = mean, na.rm = TRUE
    )
    colnames(plot_df)[colnames(plot_df) == "x"] <- metric

    # Preserve run factor ordering from the original data
    plot_df$Run <- factor(plot_df$Run, levels = levels(input_df$Run))

    p <- ggplot(plot_df,
                aes(x     = .data[["Run"]],
                    y     = .data[[metric]],
                    color = .data[["Precursor"]],
                    group = .data[["Precursor"]])) +
        geom_line(linewidth = 0.6) +
        geom_point(size = 1.5) +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme_bw() +
        theme(axis.text.x  = element_text(size = 8),
              legend.title = element_text(size = 9),
              legend.text  = element_text(size = 7)) +
        labs(x        = "Run (temporal order)",
             y        = metric,
             title    = paste("Quality Metric:", metric),
             subtitle = which.Protein,
             color    = "Peptide_Charge")

    if (isPlotly) {
        plotly_p <- ggplotly(p)
        if (!identical(address, FALSE)) {
            save_html(plotly_p,
                      file = paste0(address, "QualityMetricsPlot.html"))
        }
        return(plotly_p)
    }

    if (!identical(address, FALSE)) {
        pdf(paste0(address, "QualityMetricsPlot.pdf"))
        print(p)
        dev.off()
    }

    p
}
