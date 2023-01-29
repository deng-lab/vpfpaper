#' Convert abundance table to prab (present/absent) table
#'
#' @param df_abundance Abundance table
#' @param df_taxa Table with all taxonomy name
#' @param minrel Min relative abundance
#' @param vlevel Taxonomy level
#' @param method_name Method name, Mock, ViroProfiler, etc
#'
#' @return DataFrame
#' @export
#'
#' @examples
create_prab <- function(df_abundance, df_taxa, minrel, vlevel, method_name) {
  df_prab <- df_abundance %>%
    dplyr::mutate_at(vars(starts_with("Sample_")), ~ ifelse(. < minrel, 0, 1)) %>%
    tidyr::pivot_longer(cols = starts_with("Sample_"), names_to = "sample_id", values_to = "prab") %>%
    dplyr::right_join(df_taxa, by = c(vlevel, "sample_id")) %>%
    dplyr::mutate(prab = replace_na(.data$prab, 0)) %>%
    dplyr::mutate(method = method_name)

  return(df_prab)
}
