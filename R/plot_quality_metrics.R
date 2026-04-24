#' Plot quality metrics from converter output
#'
#' Visualizes a quality metric column from the output of MSstats converter
#' functions against run order. When the \code{Run} column is a factor with
#' levels sorted by run order (as produced automatically when \code{runOrder}
#' is supplied to the converter), the x-axis follows that temporal ordering.
#'
#' @param input data.frame or data.table returned by an MSstatsConvert
#'   converter function (e.g. \code{SpectronauttoMSstatsFormat}).
#' @param metric character, name of the column to plot on the y-axis.
#'   Defaults to \code{"AnomalyScores"}. Must be a column of \code{input}.
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
#' Each point represents a single feature (precursor / fragment) measurement.
#' A boxplot layer summarises the distribution per run, and individual points
#' are overlaid with jitter to avoid over-plotting.
#'
#' The x-axis order is determined by the factor levels of the \code{Run}
#' column. When \code{runOrder} is passed to the converter the \code{Run}
#' column is automatically set to an ordered factor; otherwise the runs appear
#' in alphabetical order.
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
#' MSstatsQualityMetricsPlot(result)
#' MSstatsQualityMetricsPlot(result, metric = "EGDeltaRT")
#' MSstatsQualityMetricsPlot(result, isPlotly = TRUE)
#' }
MSstatsQualityMetricsPlot <- function(input, metric = "AnomalyScores",
                                      address = FALSE, isPlotly = FALSE) {
    input_df <- as.data.frame(input)

    if (!metric %in% colnames(input_df)) {
        stop(paste0(
            "Column '", metric, "' not found in input. ",
            "Available columns: ", paste(colnames(input_df), collapse = ", ")
        ))
    }
    if (!"Run" %in% colnames(input_df)) {
        stop("'Run' column not found in input.")
    }

    if (!is.factor(input_df$Run)) {
        input_df$Run <- factor(input_df$Run)
    }

    p <- ggplot(input_df, aes(x = .data[["Run"]], y = .data[[metric]])) +
        geom_boxplot(outlier.shape = NA, fill = "lightblue",
                     alpha = 0.6, width = 0.5) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 0.8,
                    color = "steelblue") +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 8)) +
        labs(
            x = "Run (temporal order)",
            y = metric,
            title = paste("Quality Metric:", metric)
        )

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
