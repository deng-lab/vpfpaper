#' Calculate Bray-Curtis dissimilarity
#'
#' @param df_collapse_mock Abundance profile of mock taxonomy
#' @param df_collapse_vpf Abundance profile of ViroProfiler taxonomy
#' @param df_collapse_brackenVPF Abundance profile of Bracken on ViroProfiler taxonomy
#' @param df_collapse_brackenSTD Abundance profile of Bracken on kraken2 STD taxonomy
#' @param vlevel Taxonomy level
#'
#' @return DataFrame
#' @export
#'
#' @examples
calc_bc <- function(df_collapse_mock, df_collapse_vpf, df_collapse_brackenVPF, df_collapse_brackenSTD, vlevel) {
  dms <- list()
  sample_ids <- colnames(df_collapse_mock)[grep("Sample_", colnames(df_collapse_mock))]
  for (sample_sel in sample_ids) {
    sel_sample_mock <- df_collapse_mock %>%
      dplyr::mutate(abundance = !!sym(sample_sel)) %>%
      dplyr::select(c(!!sym(vlevel), .data$abundance)) %>%
      dplyr::mutate(method = "Mock") %>%
      dplyr::filter(.data$abundance > 0)

    sel_sample_vpf <- df_collapse_vpf %>%
      dplyr::mutate(abundance = !!sym(sample_sel)) %>%
      dplyr::select(c(!!sym(vlevel), .data$abundance)) %>%
      dplyr::mutate(method = "ViroProfiler") %>%
      dplyr::filter(.data$abundance > 0)

    sel_sample_bracken <- df_collapse_brackenVPF %>%
      dplyr::mutate(abundance = !!sym(sample_sel)) %>%
      dplyr::select(c(!!sym(vlevel), .data$abundance)) %>%
      dplyr::mutate(method = "BrackenVPF") %>%
      dplyr::filter(.data$abundance > 0)

    sel_sample_bracken_ref <- df_collapse_brackenSTD %>%
      dplyr::mutate(abundance = !!sym(sample_sel)) %>%
      dplyr::select(c(!!sym(vlevel), .data$abundance)) %>%
      dplyr::mutate(method = "BrackenSTD") %>%
      dplyr::filter(.data$abundance > 0)

    sel_sample_merge <- dplyr::bind_rows(sel_sample_mock, sel_sample_vpf, sel_sample_bracken, sel_sample_bracken_ref) %>%
      dplyr::mutate(!!sym(vlevel) := ifelse(!!sym(vlevel)=="", "unclassified", !!sym(vlevel)))

    # convert to wider format using pivot_wider, with method as rows and vlevel as columns
    sel_sample_merge_wide <- sel_sample_merge %>%
      tidyr::pivot_wider(names_from = vlevel, values_from = "abundance") %>%
      tibble::column_to_rownames("method") %>%
      dplyr::mutate(across(where(is.numeric), replace_na, 0))

    # Calculate Bray-Curtis dissimilarity
    dissimilarity_matrix <- vegan::vegdist(sel_sample_merge_wide, method = "bray")

    dm <- dissimilarity_matrix %>%
      as.matrix() %>%
      data.frame() %>%
      dplyr::select(.data$Mock) %>%
      data.table::setnames("Mock", sample_sel)

    dms[[sample_sel]] <- dm
  }

  dm <- do.call(cbind, dms)
  dm_long <- dm %>%
    tibble::rownames_to_column("method") %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("Sample_"), names_to = "sample_id", values_to = "bcs") %>%
    dplyr::filter(.data$method != "Mock") %>%
    dplyr::mutate(taxalevel = vlevel)

  return(dm_long)
}
