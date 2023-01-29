#' Read mock abundance file
#'
#' @param fpath PATH that contains *_comp.tsv files
#'
#' @return A dataframe hello
#' @export
#'
#' @examples
read_mock <- function(fpath) {
  # get sample name
  sample_id <- str_replace_all(fpath, ".+/(.+)_comp.tsv", "\\1")
  df <- fread(fpath) %>%
    setnames(c("## Virus", "Coverage"), c("gi_number", sample_id)) %>%
    dplyr::select(-.data$Name)
  return(df)
}
