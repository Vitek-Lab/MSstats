
.makeProfilePlot = function(input, is_censored, type) {
    if (type == "TRANSITION") {
        type_color_type = "FEATURE"
    } else {
        type_color = "PEPTIDE"
    }
    
    profile_plot = ggplot(input, aes_string(x='RUN', y='ABUNDANCE', 
                                            color=type_color, linetype='FEATURE')) +
        facet_grid(~LABEL) +
        geom_line(size=0.5)
    
    if (is_censored) {
        profile_plot = profile_plot +
            geom_point(aes_string(x='RUN', y='ABUNDANCE', color=type_color, shape='censored'), data=input, 
                       size=dot.size.profile) +
            scale_shape_manual(values=c(16, 1), 
                               labels=c("Detected data", "Censored missing data"))
    } else {
        profile_plot = profile_plot +
            geom_point(size=dot.size.profile) +
            scale_shape_manual(values=c(16))
    }
    
    if (type == "FEATURE" | (type == "NA" & !is_censored)) {
        profile_plot == profile_plot + scale_colour_manual(values=cbp[s])
    } else if (type == "PEPTIDE" | (type == "NA" & is_censored)) {
        profile_plot = profile_plot + scale_colour_manual(values=cbp[1:length(unique(s))])
    }
    
    profile_plot = profile_plot +
        scale_linetype_manual(values=ss) +
        scale_x_continuous('MS runs', breaks=cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup)) +
        geom_vline(xintercept=lineNameAxis + 0.5, colour="grey", linetype="longdash") +
        labs(title=unique(input$PROTEIN)) +
        geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), 
                  size=text.size, 
                  angle=text.angle, 
                  color="black") +
        theme(
            panel.background=element_rect(fill='white', colour="black"),
            legend.key=element_rect(fill='white', colour='white'),
            panel.grid.minor = element_blank(),
            strip.background=element_rect(fill='gray95'),
            strip.text.x=element_text(colour=c("#00B0F6"), size=14),
            axis.text.x=element_text(size=x.axis.size, colour="black"),
            axis.text.y=element_text(size=y.axis.size, colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
            axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
            title=element_text(size=x.axis.size+8, vjust=1.5),
            legend.position="top",
            legend.text=element_text(size=legend.size))
    
    color_guide = guide_legend(title=paste("# peptide:", nlevels(input$PEPTIDE)), 
                               title.theme = element_text(size=13, angle=0),
                               keywidth=0.1,
                               keyheight = 0.1,
                               default.unit = 'inch',
                               ncol=3)
    linetype_guide = guide_legend(title=paste("# peptide:", nlevels(input$PEPTIDE)), 
                                  title.theme = element_text(size=13, angle=0),
                                  keywidth=0.1,
                                  keyheight = 0.1,
                                  default.unit = 'inch',
                                  ncol=3)
    shape_guide = guide_legend(title=NULL,
                               label.theme = element_text(size=11, angle=0),
                               keywidth=0.1,
                               keyheight = 0.1,
                               default.unit = 'inch')
    if (is_censored) {
        if (type == "FEATURE") {
            profile_plot = profile_plot +
                guides(color = color_guide,
                       linetype = linetype_guide,
                       shape = shape_guide)
        } else if (type == "PEPTIDE") {
            profile_plot = profile_plot + guides(color = color_guide)
        } else {
            profile_plot = profile_plot + guides(shape = shape_guide)
        }
    } else {
        if (type == "FEATURE") {
            profile_plot = profile_plot +
                guides(color = color_guide, 
                       linetype = linetype_guide)
        } else if (type == "PEPTIDE") {
            profile_plot = profile_plot + 
                guides(color = color_guide, shape = shape_guide)
        }
    }
    profile_plot    
}


.makeSummaryProfilePlot = function(input, is_censored) {
    profile_plot = ggplot(aes_string(x='RUN', y='ABUNDANCE', 
                                     color='analysis', linetype='FEATURE', size='analysis'), data=input) +
        facet_grid(~LABEL) +
        geom_line(size=0.5)
    
    if (is_censored) {
        profile_plot = profile_plot +
            geom_point(aes_string(x='RUN', y='ABUNDANCE', 
                                  color='analysis', size='analysis', shape='censored'), data=input) +
            scale_shape_manual(values=c(16, 1), labels=c("Detected data", "Censored missing data"))
    } else {
        profile_plot = profile_plot +         
            geom_point(size=dot.size.profile) +
            scale_shape_manual(values=c(16))
        
    }
    
    profile_plot = profile_plot +
        scale_colour_manual(values=c("lightgray", "darkred")) +
        scale_size_manual(values=c(1.7, 2), guide="none") +
        scale_linetype_manual(values=c(rep(1, times=length(unique(input$FEATURE))-1), 2), guide="none") +
        scale_x_continuous('MS runs',breaks=cumGroupAxis) +
        scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup)) +
        geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash") +
        labs(title=unique(input$PROTEIN)) +
        geom_text(data=groupNametemp, aes(x=RUN, y=ABUNDANCE, label=Name), 
                  size=text.size, 
                  angle=text.angle, 
                  color="black") +
        theme(
            panel.background=element_rect(fill='white', colour="black"),
            legend.key=element_rect(fill='white', colour='white'),
            panel.grid.minor = element_blank(),
            strip.background=element_rect(fill='gray95'),
            strip.text.x=element_text(colour=c("#00B0F6"), size=14),
            axis.text.x=element_text(size=x.axis.size, colour="black"),
            axis.text.y=element_text(size=y.axis.size, colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
            axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
            title=element_text(size=x.axis.size+8, vjust=1.5),
            legend.position="top",
            legend.text=element_text(size=legend.size),
            legend.title=element_blank())
    
    color_guide = guide_legend(order=1,
                               title=NULL,
                               label.theme = element_text(size=10, angle=0))
    shape_guide = shape=guide_legend(order=2, 
                                     title=NULL,
                                     label.theme = element_text(size=10, angle=0))
    if (is_censored) {
        profile_plot = profile_plot +
            guides(color = color_guide, shape = shape_guide)
    } else {
        profile_plot = profile_plot +
            guides(color = color_guide) +
            geom_point(aes_string(x = "RUN", y = "ABUNDANCE", size = "analysis",
                                  color = "analysis"), data = input)
    }
    profile_plot
}


.makeQCPlot = function(input, all_proteins) {
    if (all_proteins) {
        plot_title = "All proteins"
    } else {
        plot_title = unique(input$PROTEIN)
    }
    
    ggplot(input, aes_string(x='RUN', y='ABUNDANCE'))+
        facet_grid(~LABEL)+
        geom_boxplot(aes_string(fill='LABEL'), outlier.shape=1, outlier.size=1.5)+
        scale_fill_manual(values=label.color, guide="none")+
        scale_x_discrete('MS runs', breaks=cumGroupAxis)+
        scale_y_continuous(yaxis.name, limits=c(y.limdown, y.limup))+
        geom_vline(xintercept=lineNameAxis+0.5, colour="grey", linetype="longdash")+
        labs(title = plot_title)+
        geom_text(data=groupName, aes(x=RUN, y=ABUNDANCE, label=Name), 
                  size=text.size, angle=text.angle, color="black")+
        theme(
            panel.background=element_rect(fill='white', colour="black"),
            legend.key=element_rect(fill='white', colour='white'),
            panel.grid.minor = element_blank(),
            strip.background=element_rect(fill='gray95'),	
            strip.text.x=element_text(colour=c("#00B0F6"), size=14),
            axis.text.x=element_text(size=x.axis.size, colour="black"),
            axis.text.y=element_text(size=y.axis.size, colour="black"),
            axis.ticks=element_line(colour="black"),
            axis.title.x=element_text(size=x.axis.size+5, vjust=-0.4),
            axis.title.y=element_text(size=y.axis.size+5, vjust=0.3),
            title=element_text(size=x.axis.size+8, vjust=1.5))
}

.getQCProtein = function(input, protein) {
    if (is.numeric(protein)) {
        protein = levels(input$PROTEIN)[protein]
    }
    if (!(protein %in% levels(input$PROTEIN))) {
        stop(paste0("Please check protein name. Data set does not have this protein. - ", 
                    toString(protein)))
    }
    protein
}