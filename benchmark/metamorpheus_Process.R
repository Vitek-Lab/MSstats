calculate_Metrics <- function(QuantData, protein_mappings, task_label, alpha = 0.05) {
  comparison <- matrix(
    c(-1,0,0,0,1,   # E-A
      -1,0,0,1,0,   # D-A
      -1,0,1,0,0,   # C-A
      -1,1,0,0,0),  # B-A
    nrow = 4, byrow = TRUE
  )
  rownames(comparison) <- c("E-A", "D-A", "C-A", "B-A")
  groups <- levels(QuantData$ProteinLevelData$GROUP)
  colnames(comparison) <- groups[order(as.numeric(groups))]

  model <- groupComparison(
    contrast.matrix = comparison,
    data = QuantData,
    use_log_file = FALSE
  )

  ecoli_ids <- protein_mappings %>%
    filter(Organism == "Escherichia coli (strain K12)") %>%
    pull(`Protein Groups`)

  filtered_results <- model$ComparisonResult %>%
    mutate(ecoli = Protein %in% ecoli_ids) %>%
    filter(is.na(issue))

  labels <- unique(filtered_results$Label)
  result_rows <- lapply(labels, function(lbl) {
    df <- filtered_results %>% filter(Label == lbl)
    sig <- df %>% filter(adj.pvalue < alpha)

    tp <- sig %>% filter(ecoli) %>% nrow()
    fp <- sig %>% filter(!ecoli) %>% nrow()
    tot <- tp + fp
    fdr <- if (tot > 0) fp / tot else NA_real_

    data.frame(
      Task       = task_label,
      Comparison = lbl,
      FDR        = fdr,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, result_rows)
  return(results)
}
