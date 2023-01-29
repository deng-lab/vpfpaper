#' Read CheckV results
#'
#' @param fpath Path of CheckV tsv file
#'
#' @return DataFrame
#' @export
#'
#' @examples
read_checkv <- function(fpath) {
  . <- NULL
  df <- fread(fpath) %>%
    mutate(checkv_quality = factor(checkv_quality, levels=c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"))) %>%
    column_to_rownames("contig_id") %>%
    setnames(colnames(.), paste0("checkv_", colnames(.))) %>%
    setnames("checkv_checkv_quality", "checkv_quality") %>%
    rownames_to_column("Contig")
  return(df)
}
