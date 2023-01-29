#' Read Bracken abundance report file
#'
#' @param fin_ab_bracken Bracken abundance report file
#' @param t2lin Mapping file of taxonomy ID <-> lineage
#' @param cnames Selected column names
#' @param vlevel Taxonomy level
#' @param remove_unclassified Merge unclassified taxa or not, default 0
#'
#' @return DataFrame
#' @export
#'
#' @examples
create_ab_bracken <- function(fin_ab_bracken, t2lin, cnames, vlevel, remove_unclassified=0) {
  . <- NULL
  df_collapse_bracken <- fread(fin_ab_bracken) %>%
    # selct columns ending with "_frac" and "taxonomy_id"
    dplyr::select(c(ends_with("_frac"), .data$taxonomy_id)) %>%
    # rename columns
    setnames(colnames(.), str_replace_all(colnames(.), ".tsv_frac", "")) %>%
    setnames(colnames(.), str_replace_all(colnames(.), "ds10Ms", "Sample_")) %>%
    setnames(colnames(.), str_replace_all(colnames(.), "Sample_0", "Sample_")) %>%
    setnames("taxonomy_id", "taxid") %>%
    # mutate value in columns starting with "Sample_" to multiply by 100
    mutate_at(vars(starts_with("Sample_")), ~ . * 100) %>%
    left_join(t2lin, by = "taxid") %>%
    dplyr::select(c(starts_with("Sample_"), !!sym(vlevel))) %>%
    dplyr::mutate(!!sym(vlevel) := ifelse(is.na(!!sym(vlevel)), "unclassified", !!sym(vlevel))) %>%
    group_by(!!sym(vlevel)) %>%
    summarise_all(sum)

  # reorder columns
  df_collapse_bracken <- df_collapse_bracken[, cnames]

  # If remove unclassified
  if (remove_unclassified!=0) {
    df_collapse_bracken <- df_collapse_bracken %>%
      dplyr::filter(str_detect(!!sym(vlevel), "^unclassified", negate = TRUE))
  }

  return(df_collapse_bracken)
}


#' Abundance profile of Bracken using standard Kraken2 database
#'
#' @param fin_t2lin_krakenDB Mapping file of Kraken2 STD reports taxid <-> lineage
#' @param fin_ab_brackenSTD Abundance profile of Bracken using STD
#' @param cnames Selected column names
#' @param vlevel taxonomy level
#' @param remove_unclassified Merge unclassified taxa or not, default 0
#'
#' @return DataFrame
#' @export
#'
#' @examples
create_ab_brackenSTD <- function(fin_t2lin_krakenDB, fin_ab_brackenSTD, cnames, vlevel, remove_unclassified=0) {
  t2lin_brackenSTD <- fread(fin_t2lin_krakenDB)
  df_brackenSTD_collapse <- create_ab_bracken(fin_ab_brackenSTD, t2lin_brackenSTD, cnames, vlevel, remove_unclassified)
  return(df_brackenSTD_collapse)
}


#' Abundance profile of Bracken using ViroProfiler taxonomy
#'
#' @param df_vpf_taxa dataframe of ViroProfiler taxonomy
#' @param fin_brackenVPF file of ViroProfiler bracken abundance table
#' @param cnames Selected column names
#' @param vlevel taxonomy level
#' @param remove_unclassified Merge unclassified taxa or not, default 0
#'
#' @return DataFrame
#' @export
#'
#' @examples
create_ab_brackenVPF <- function(df_vpf_taxa, fin_brackenVPF, cnames, vlevel, remove_unclassified=0) {
  t2lin_brackenVPF <- df_vpf_taxa[, c("taxid", "phylum", "class", "order", "family", "genus", "species", "subspecies")]
  df_brackenVPF_collapse <- create_ab_bracken(fin_brackenVPF, t2lin_brackenVPF, cnames, vlevel, remove_unclassified)
  return(df_brackenVPF_collapse)
}

