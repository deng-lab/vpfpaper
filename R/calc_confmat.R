#' Calculate confusion matrix
#'
#' @param df_collapse_mock Abundance profile of mock taxonomy
#' @param df_collapse_vpf Abundance profile of ViroProfiler taxonomy
#' @param df_collapse_brackenVPF Abundance profile of Bracken on ViroProfiler taxonomy
#' @param df_collapse_brackenSTD Abundance profile of Bracken on kraken2 STD taxonomy
#' @param vlevel Taxonomy level
#' @param minrel Min relative abundance
#'
#' @return DataFrame
#' @export
#'
#' @examples
calc_confmat <- function(df_collapse_mock, df_collapse_vpf, df_collapse_brackenVPF, df_collapse_brackenSTD, vlevel, minrel) {
  # Create a data frame with all combinations of taxa and sample names
  df_taxa <- expand.grid(taxonomy = unique(c(df_collapse_mock[[vlevel]], df_collapse_vpf[[vlevel]], df_collapse_brackenVPF[[vlevel]], df_collapse_brackenSTD[[vlevel]])),
                         sample_id = colnames(df_collapse_mock)[grep("Sample_", colnames(df_collapse_mock))]) %>%
    setnames("taxonomy", vlevel)

  prab_mock <- create_prab(df_collapse_mock, df_taxa, minrel, vlevel, "Mock")
  prab_vpf <- create_prab(df_collapse_vpf, df_taxa, minrel, vlevel, "ViroProfiler")
  prab_brackenVPF <- create_prab(df_collapse_brackenVPF, df_taxa, minrel, vlevel, "BrackenVPF")
  prab_brackenSTD <- create_prab(df_collapse_brackenSTD, df_taxa, minrel, vlevel, "BrackenSTD")

  df_prab <- rbind(prab_mock, prab_vpf, prab_brackenVPF, prab_brackenSTD) %>%
    # convert to factor, required by caret::confusionMatrix
    mutate(!!sym(vlevel) := factor(!!sym(vlevel), levels = unique(!!sym(vlevel)))) %>%
    mutate(prab = factor(.data$prab, levels = c(0, 1))) %>%
    pivot_wider(names_from = "method", values_from = "prab")

  metrics_list <- list()
  # Calculate confusion matrix for each sample + method combination
  for (sample_name in unique(df_prab$sample_id)) {
    for (method_name in c("ViroProfiler", "BrackenVPF", "BrackenSTD")) {
      df_sub <- df_prab %>%
        dplyr::filter(.data$sample_id == sample_name) %>%
        dplyr::select(-.data$sample_id) %>%
        dplyr::select(c(!!sym(vlevel), .data$Mock, !!sym(method_name)))

      confusion_matrix <- caret::confusionMatrix(df_sub[[method_name]], df_sub[["Mock"]], positive = "1")
      metrics_list[[paste0(sample_name, "_", method_name)]] <-
        data.frame(
          sample_id = sample_name,
          method = method_name,
          precision = confusion_matrix$byClass["Precision"],
          recall = confusion_matrix$byClass["Recall"],
          f1_score = confusion_matrix$byClass["F1"]
        )
    }
  }

  df_confmat <- do.call(rbind, metrics_list) %>%
    pivot_longer(cols = c("precision", "recall", "f1_score"), names_to = "metric", values_to = "value") %>%
    mutate(threshold = minrel)

  return(df_confmat)
}
