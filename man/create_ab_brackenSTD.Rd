% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_ab_bracken.R
\name{create_ab_brackenSTD}
\alias{create_ab_brackenSTD}
\title{Abundance profile of Bracken using standard Kraken2 database}
\usage{
create_ab_brackenSTD(
  fin_t2lin_krakenDB,
  fin_ab_brackenSTD,
  cnames,
  vlevel,
  remove_unclassified = 0
)
}
\arguments{
\item{fin_t2lin_krakenDB}{Mapping file of Kraken2 STD reports taxid <-> lineage}

\item{fin_ab_brackenSTD}{Abundance profile of Bracken using STD}

\item{cnames}{Selected column names}

\item{vlevel}{taxonomy level}

\item{remove_unclassified}{Merge unclassified taxa or not, default 0}
}
\value{
DataFrame
}
\description{
Abundance profile of Bracken using standard Kraken2 database
}
