#' Taxonomy abundance profile in PRAB (present/absent) format
#'
#' Abundance profile of a specific taxonomy level
#' Report ...
#'
#' @format
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{method}{ViroProfiler, BrackenVPF, BrackenSTD, Mock}
#'   \item{metric}{precision, recall, f1_score}
#'   \item{sample_id}{Sample ID}
#'   \item{taxalevel}{Taxonomy level}
#'   \item{threshold}{Threshold of relative abundance}
#'   \item{value}{Value of precision, recall, f1_score}
#'   ...
#' }
#' @source {Created in-house PAB profile}
#' data(vpf_prab)
"vpf_prab"


#' Taxonomy abundance profile in RELAB (Relative abundance) format
#'
#' Abundance profile of a specific taxonomy level
#' Report ...
#'
#' @format
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{bcs}{Bray-Curtis dissimilarity}
#'   \item{method}{ViroProfiler, BrackenVPF, BrackenSTD}
#'   \item{sample_id}{Sample ID}
#'   \item{taxalevel}{Taxonomy level}
#'   ...
#' }
#' @source {Created in-house PAB profile}
#' data(vpf_relab)
"vpf_relab"

