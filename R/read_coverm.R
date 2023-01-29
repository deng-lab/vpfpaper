#' Read CoverM output file
#'
#' @param fpath CoverM output file name
#' @param fbin2contig Contig binning to contig mapping
#' @param cov Coverage fraction file or not
#'
#' @return A dataframe
#' @export
#'
#' @examples
read_coverm <- function(fpath, fbin2contig = 0, cov = 0) {
  . <- NULL
  df_abundance <- fread(fpath) %>%
    setnames(colnames(.), str_replace_all(colnames(.), "ds10Ms", "Sample_")) %>%
    setnames(colnames(.), str_replace_all(colnames(.), "Sample_0", "Sample_")) %>%
    setnames("Contig", "genome_id")

  # check if bin2contig file is provided
  if (fbin2contig!=0) {
    df_bin2contig <- fread(fbin2contig, header = FALSE, col.names = c("bin_id", "genome_id"))
    df_abundance <- df_abundance %>%
      dplyr::left_join(df_bin2contig, by = "genome_id") %>%
      dplyr::mutate(bin_id = ifelse(is.na(.data$bin_id), .data$genome_id, .data$bin_id)) %>%
      dplyr::mutate(genome_id = .data$bin_id) %>%
      dplyr::select(-.data$bin_id) %>%
      dplyr::group_by(.data$genome_id)

    # If read coverage fraction, then choose the maximum value for each bin
    if (cov!=0) {
      df_abundance <- df_abundance %>% summarise_all(max)
    } else {
      df_abundance <- df_abundance %>% summarise_all(sum)
    }
  }
  return(df_abundance)
}
