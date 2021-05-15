#' Visualization for model-based quality control in fitting model
#' 
#' @description To check the assumption of linear model for whole plot inference, 
#' modelBasedQCPlots takes the results after fitting models from function 
#' (\code{\link{groupComparison}}) as input and automatically generate two types 
#' of figures in pdf files as output: 
#' (1) normal quantile-quantile plot (specify "QQPlot" in option type) for checking
#' normally distributed errors.; 
#' (2) residual plot (specify "ResidualPlot" in option type).
#'  
#' @param data output from function groupComparison.
#' @param type choice of visualization. "QQPlots" represents normal quantile-quantile 
#' plot for each protein after fitting models. "ResidualPlots" represents a plot
#' of residuals versus fitted values for each protein in the dataset.
#' @param axis.size size of axes labels. Default is 10.
#' @param dot.size size of points in the graph for residual plots and QQ plots. Default is 3.
#' @param width width of the saved file. Default is 10.
#' @param height height of the saved file. Default is 10.
#' @param which.Protein Protein list to draw plots. List can be names of Proteins 
#' or order numbers of Proteins from levels(testResultOneComparison$ComparisonResult$Protein). 
#' Default is "all", which generates all plots for each protein.
#' @param address name that will serve as a prefix to the name of output file.
#' 
#' @details Results based on statistical models for whole plot level inference are 
#' accurate as long as the assumptions of the model are met. The model assumes that 
#' the measurement errors are normally distributed with mean 0 and constant variance. 
#' The assumption of a constant variance can be checked by examining the residuals from the model.
#' \itemize{
#' \item{QQPlots : a normal quantile-quantile plot for each protein is generated in order to check whether the errors are well approximated by a normal distribution. If points fall approximately along a straight line, then the assumption is appropriate for that protein. Only large deviations from the line are problematic.}
#' \item{ResidualPlots : The plots of residuals against predicted(fitted) values. If it shows a random scatter, then the assumption is appropriate.}
#' }
#' 
#' @return produce a pdf file
#' 
#' @export
#' 
#' @examples
#' QuantData <- dataProcess(SRMRawData, use_log_file = FALSE)
#' head(QuantData$FeatureLevelData)
#' levels(QuantData$FeatureLevelData$GROUP)
#' comparison <- matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
#' row.names(comparison) <- "T7-T1"
#' colnames(comparison) <- unique(QuantData$ProteinLevelData$GROUP)
#' # Tests for differentially abundant proteins with models:
#' # label-based SRM experiment with expanded scope of biological replication.
#' testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData,
#' use_log_file = FALSE)
#' # normal quantile-quantile plots
#' modelBasedQCPlots(data=testResultOneComparison, type="QQPlots", address="")
#' # residual plots
#' modelBasedQCPlots(data=testResultOneComparison, type="ResidualPlots", address="")
#' 
modelBasedQCPlots = function(
    data, type, axis.size = 10, dot.size = 3, width = 10, height = 10, 
    which.Protein = "all", address = ""
) {
    if (length(setdiff(toupper(type), c("QQPLOTS","RESIDUALPLOTS"))) != 0) {
        stop(paste0("Input for type=", type, 
                    ". However,'type' should be one of QQPlots, ResidualPlots."))
    }
    if (address == FALSE) {
        if(all(which.Protein == 'all')){
            stop('** Cannnot generate all plots in a screen. Please set one protein at a time.')
        } else if (length(which.Protein) > 1) {
            stop('** Cannnot generate multiple plots in a screen. Please set one protein at a time.')
        }
    }
    
    fitted_models = data[["FittedModel"]]
    all_proteins = levels(data$ComparisonResult$Protein)
    if (all(which.Protein != "all")) {
        selected_proteins = getSelectedProteins(which.Protein, all_proteins)
        fitted_models = fitted_models[all_proteins %in% selected_proteins]
        all_proteins = all_proteins[all_proteins %in% selected_proteins]
    }
    
    if (toupper(type) == "QQPLOTS") {
        .plotQQ(fitted_models, all_proteins, width, height, address, 
                dot.size, axis.size)
    } else if (toupper(type) == "RESIDUALPLOTS") {
        .plotResiduals(fitted_models, all_proteins, width, height,
                       address, dot.size, axis.size)
    }
}


#' @importFrom stats resid quantile qnorm 
#' @importFrom utils setTxtProgressBar
.plotQQ = function(fitted_models, all_proteins, width, height, address,
                   dot.size, axis.size) {
    residual = NULL
    savePlot(address, "QQPlot", width, height)
    pb = utils::txtProgressBar(min = 0, max = length(fitted_models), style = 3)
    for (i in seq_along(fitted_models)) {	
        sub = fitted_models[[i]]
        if(is.null(sub)){
            next()
        }
        if (is(sub, "lm")) {
            sub.residuals = sub$residuals
        } else {
            sub.residuals = resid(sub)
        }
        sub.residuals.table = data.frame("residual" = sub.residuals)
        ## get slope and intercept for qline
        y = quantile(sub.residuals.table$residual[!is.na(sub.residuals.table$residual)],
                     c(0.25, 0.75))
        x = qnorm(c(0.25, 0.75))
        slope = diff(y) / diff(x)
        int = y[1L] - slope * x[1L]
        
        plot = ggplot(sub.residuals.table, aes(sample = residual)) +
            geom_point(stat = "qq", alpha = 0.8, 
                       shape = 20, size = dot.size) +
            scale_shape(solid=FALSE) + 
            geom_abline(slope = slope, intercept = int, colour = "red") +
            scale_y_continuous("Sample Quantiles") +
            scale_x_continuous("Theoretical Quantiles") +
            labs(title=paste("Normal Q-Q Plot (", all_proteins[i], ")")) +
            theme(
                panel.background = element_rect(fill = "white", colour = "black"),
                panel.grid.major = element_line(colour = "grey95"),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(size = axis.size, colour = "black"),
                axis.text.y = element_text(size = axis.size, colour = "black"),
                axis.ticks = element_line(colour = "black"),
                axis.title.x = element_text(size = axis.size+5, vjust = -0.4),
                axis.title.y = element_text(size = axis.size+5, vjust = 0.3),
                title = element_text(size = axis.size + 8, vjust = 1.5),
                legend.position = "none"
            )
        print(plot)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    if (address != FALSE) {
        dev.off()
    }
}


#' @importFrom stats resid fitted
#' @importFrom utils setTxtProgressBar
.plotResiduals = function(fitted_models, all_proteins, width, height, 
                          address, dot.size, axis.size) {
    savePlot(address, "ResidualPlot", width, height)
    pb = utils::txtProgressBar(min = 0, max = length(fitted_models), style = 3)
    for (i in seq_along(fitted_models)) {	
        fitted_model = fitted_models[[i]]
        if(is.null(fitted_model)) {
            next()
        }
        if (is(fitted_model, "lm")) {
            sub.residuals = fitted_model$residuals
            sub.fitted = fitted_model$fitted.values
        } else {
            sub.residuals = resid(fitted_model)
            sub.fitted = fitted(fitted_model)
        }
        sub.residuals.table = data.frame("residual" = sub.residuals, 
                                         "fitted" = sub.fitted)
        
        plot = ggplot(aes_string(x = "fitted", y = "residual"), 
                       data = sub.residuals.table) +
            geom_point(size = dot.size, 
                       alpha = 0.5) +
            geom_hline(yintercept = 0, 
                       linetype = "twodash", 
                       colour = "darkgrey", 
                       size = 0.6) +
            scale_y_continuous("Residuals") +
            scale_x_continuous("Predicted Abundance") +
            labs(title = all_proteins[i]) +
            theme(
                panel.background = element_rect(fill = "white", colour = "black"),
                panel.grid.major = element_line(colour = "grey95"),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(size = axis.size, colour="black"),
                axis.text.y = element_text(size = axis.size, colour="black"),
                axis.ticks = element_line(colour = "black"),
                axis.title.x = element_text(size = axis.size + 5, vjust = -0.4),
                axis.title.y = element_text(size = axis.size + 5, vjust = 0.3),
                title = element_text(size = axis.size + 8, vjust = 1.5),
                legend.position = "none")
        print(plot)
        setTxtProgressBar(pb, i)
    }
    close(pb)
    if (address != FALSE) {
        dev.off()
    }
}
