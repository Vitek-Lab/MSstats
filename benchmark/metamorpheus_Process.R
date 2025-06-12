calculate_Metrics(QuantData, Label, protein_mappings){

    # dataProcessPlots(QuantData, "QCPlot", which.Protein = "allonly", address = FALSE, isPlotly = TRUE)

    comparison <- matrix(c(-1,0,0,0,1, # 3x
                        -1,0,0,1,0, # 2.5x
                        -1,0,1,0,0, # 2x
                        -1,1,0,0,0),nrow=4,byrow = TRUE) # 1.5x
                        
    row.names(comparison) <- c("E-A", "D-A", "C-A", "B-A")
    groups = levels(QuantData$ProteinLevelData$GROUP)
    colnames(comparison) <- groups[order(as.numeric(groups))]
    model <- groupComparison(contrast.matrix=comparison, data=QuantData,
                                    use_log_file = FALSE)
                    
                                    
                                    
    ecoli_proteins = protein_mappings %>% filter(Organism == "Escherichia coli (strain K12)")
    model$ComparisonResult = model$ComparisonResult %>% mutate(ecoli = Protein %in% ecoli_proteins$`Protein Groups`)


    e_group = model$ComparisonResult %>% filter(Label == "B-A") %>% filter(is.na(issue))

    ecoli = e_group %>% filter(ecoli == TRUE)

    # hist(ecoli$log2FC)

    ecoli = e_group %>% filter(adj.pvalue < 0.05) %>% filter(ecoli == TRUE)
    human = e_group %>% filter(adj.pvalue < 0.05) %>% filter(ecoli == FALSE)
    FDR = nrow(human) / (nrow(ecoli) + nrow(human))

    cat(label, FDR, "\n")

    results <- data.frame(
    Label = label,
    FDR = fdr
  )
  
  return(results)

}