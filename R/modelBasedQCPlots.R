#' @export
modelBasedQCPlots = function(
    data, type, axis.size = 10, dot.size = 3, text.size = 7, legend.size = 7,
    width = 10, height = 10, which.Protein = "all", address = ""
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
    
    fitted_models = data[["fittedmodel"]]
    all_proteins = levels(data$ComparisonResult$Protein)
    if (all(which.Protein != "all")) {
        selected_proteins = .getSelectedProteins(which.Protein, all_proteins)
        fitted_models = fitted_models[all_proteins %in% selected_proteins]
        all_proteins = all_proteins[all_proteins %in% selected_proteins]
    }
    
    if (toupper(type) == "QQPLOTS") {
        .plotQQ(fitted_models, all_proteins, width, height, address, 
                dot.size, axis.size)
    } else if (toupper(type) == "RESIDUALPLOTS") {
        .plotResiduals(fitted_models, all_proteins, width, height,
                       address, axis.size)
    }
}

.plotQQ = function(fitted_models, all_proteins, width, height, address,
                   dot.size, axis.size) {
    # normality
    .savePlot(address, "QQPlot", width, height)
    pb = txtProgressBar(1, length(fitted_models))
    for (i in seq_along(fitted_models)) {	
        sub = fitted_models[[i]]
        if(is.null(sub)){
            next()
        }
        if (class(sub) == "lm") {
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
                panel.background = element_rect(fill="white", colour="black"),
                panel.grid.major = element_line(colour="grey95"),
                panel.grid.minor = element_blank(),
                axis.text.x = element_text(size=axis.size, colour="black"),
                axis.text.y = element_text(size=axis.size, colour="black"),
                axis.ticks = element_line(colour="black"),
                axis.title.x = element_text(size=axis.size+5, vjust=-0.4),
                axis.title.y = element_text(size=axis.size+5, vjust=0.3),
                title = element_text(size=axis.size+8, vjust=1.5),
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


.plotResiduals = function(fitted_models, all_proteins, width, height, 
                          address, axis.size) {
    .savePlot(address, "ResidualPlot", width, height)
    pb = txtProgressBar(0, length(fitted_models))
    for (i in seq_along(fitted_models)) {	
        sub = fitted_models[[i]]
        if(is.null(sub)) {
            next()
        }
        if (class(sub)=="lm") {
            sub.residuals = sub$residuals
            sub.fitted = sub$fitted.values
        } else {
            sub.residuals = resid(sub)
            sub.fitted = fitted(sub)
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
