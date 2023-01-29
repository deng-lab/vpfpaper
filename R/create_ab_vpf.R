#' Abundance profile of ViroProfiler
#'
#' @param df_vpf_taxa dataframe of ViroProfiler taxonomy
#' @param fin_vpf_count file of ViroProfiler count table
#' @param fin_vpf_cov file of ViroProfiler coverage fraction table
#' @param cnames Selected column names
#' @param vlevel taxonomy level
#' @param mincov ViroProfiler minimum coverage fraction 0-1
#' @param bin2contig Mapping file of bins and contigs, default 0
#' @param remove_unclassified Merge unclassified taxa or not, default 0
#' @param fctglen Mapping file of contig length
#'
#' @return DataFrame
#' @export
#'
#' @examples
create_ab_vpf <- function(df_vpf_taxa, fin_vpf_count, fin_vpf_cov, cnames, vlevel, mincov = 0.1, bin2contig = 0, remove_unclassified=0, fctglen=0) {
  . <- NULL
  df_vpf_cov <- read_coverm(fin_vpf_cov, bin2contig, "cov") %>%
    column_to_rownames("genome_id") %>%
    # if value smaller than mincov, set to 0, otherwise set to 1
    mutate_at(vars(starts_with("Sample_")), ~ ifelse(. < mincov, 0, 1)) %>%
    as.matrix()

  df_vpf_count <- read_coverm(fin_vpf_count, bin2contig) %>%
    column_to_rownames("genome_id") %>%
    as.matrix()

  # Normalize by contig length
  if (fctglen!=0) {
    df_ctglen <- fread(fctglen, col.names = c("ctgid", "ctglen")) %>%
      column_to_rownames("ctgid")
    df_ctglen <- df_ctglen[rownames(df_vpf_count),]
    df_vpf_count <- df_vpf_count / df_ctglen
    df_vpf_count <- round(df_vpf_count, 6)
  }

  # multiply coverage matrix with count matrix
  df_vpf_cov_count <- df_vpf_cov * df_vpf_count %>%
    data.frame() %>%
    dplyr::mutate_at(vars(starts_with("Sample_")), ~ . / sum(.) * 100) %>%
    dplyr::mutate_at(vars(starts_with("Sample_")), ~ round(., 3))

  df_vpf_cov_count <- df_vpf_cov_count %>%
    rownames_to_column("genome_id")
  #
  df_collapse_vpf <- df_vpf_cov_count %>%
    left_join(df_vpf_taxa, by = "genome_id") %>%
    dplyr::select(c(starts_with("Sample_"), !!sym(vlevel))) %>%
    dplyr::mutate(!!sym(vlevel) := ifelse(!!sym(vlevel) == "", "unclassified", !!sym(vlevel))) %>%
    dplyr::mutate(!!sym(vlevel) := ifelse(is.na(!!sym(vlevel)), "unclassified", !!sym(vlevel))) %>%
    dplyr::group_by(!!sym(vlevel)) %>%
    summarise_all(sum) %>%
    arrange(desc(rowSums(select(., starts_with("Sample_")))))

  # If remove unclassified
  if (remove_unclassified!=0) {
    df_collapse_vpf <- df_collapse_vpf %>%
      dplyr::filter(str_detect(!!sym(vlevel), "^unclassified", negate = TRUE))
  }

  df_collapse_vpf <- df_collapse_vpf[, cnames]

  return(df_collapse_vpf)
}
