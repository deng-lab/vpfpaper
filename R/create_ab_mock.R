#' Create taxonomy abundance of of mock data
#'
#' @param fin_t2lin_mock mapping file of gi <-> taxid <-> lineage
#' @param fins_ab_mock abundance files of mock data
#' @param vlevel taxonomy level
#' @param remove_unclassified default: 0
#'
#' @return A dataframe
#' @export
#'
#' @examples
create_ab_mock <- function(fin_t2lin_mock, fins_ab_mock, vlevel, remove_unclassified=0) {
  . <- NULL
  mock_taxid2lineage <- fread(fin_t2lin_mock) %>%
    dplyr::filter(!is.na(.data$taxid))
  mock_taxid2lineage <- mock_taxid2lineage[, c("gi_number", "taxid", "phylum", "class", "order", "family", "genus", "species", "subspecies")]

  df_mock_collapse <- lapply(fins_ab_mock, read_mock) %>%
    Reduce(function(x, y) merge(x, y, all = TRUE), .) %>%
    mutate(across(!tidyselect::matches("gi_number"), replace_na, 0)) %>%
    inner_join(mock_taxid2lineage, by = "gi_number") %>%
    dplyr::select(-.data$gi_number) %>%
    # mutate columns starting with "Sample_" to percentage, round to 4 digits
    mutate_at(vars(starts_with("Sample_")), ~ . / sum(.) * 100) %>%
    mutate_at(vars(starts_with("Sample_")), ~ round(., 3)) %>%
    dplyr::select(c(starts_with("Sample_"), !!sym(vlevel))) %>%
    dplyr::mutate(!!sym(vlevel) := ifelse(!!sym(vlevel) == "", "unclassified", !!sym(vlevel))) %>%
    group_by(!!sym(vlevel)) %>%
    summarise_all(sum) %>%
    # arrange by row sum of "Sample_" columns in descending order
    arrange(desc(rowSums(select(., starts_with("Sample_")))))

  # If remove unclassified
  if (remove_unclassified!=0) {
    df_mock_collapse <- df_mock_collapse %>%
      dplyr::filter(str_detect(!!sym(vlevel), "^unclassified", negate = TRUE))
  }

  return(df_mock_collapse)
}
