## code to prepare `vpf_prab` dataset goes here

library(here)
library(vpfpaper)

dpath <- "/Users/allen/github/rujinlong/viroprofiler-paper/workflow/data/00-rawdata/vpf_results/len500"

fin_t2lin_mock <- here(dpath, "taxid2lineage_mock.tsv")
fin_t2lin_VPF <- here(dpath, "taxid2lineage_vpf.tsv")
fin_t2lin_krakenDB <- here(dpath, "taxid2lineage_krakenSTD.tsv")

# input: abundance
fins_ab_mock <- list.files(here(dpath, "mockref"), pattern = "_comp.tsv$", full.names = T)
fin_ab_vpfcount <- here(dpath, "abundance_contigs_count.tsv.gz")
fin_ab_vpfcov <- here(dpath, "abundance_contigs_covered_fraction.tsv.gz")
fin_ab_brackenVPF <- here(dpath, "abundance_brackenVPF.txt")
fin_ab_brackenSTD <- here(dpath, "abundance_brackenSTD.txt")

# import data
df_taxa_vpf <- fread(fin_t2lin_VPF)

# calculate precision-recall metrics
mincov <- 0.1
remove_unclassified <- 1
dfmetrics = list()
vlevels <- c("phylum", "order", "family", "genus", "species", "subspecies")
for (vlevel in vlevels) {
  minrel_list <- seq(0.001, 0.5, 0.005)
  dfs <- lapply(minrel_list, function(minrel) {
    df_collapse_mock <- create_ab_mock(fin_t2lin_mock, fins_ab_mock, vlevel, remove_unclassified=remove_unclassified)
    df_collapse_vpf <- create_ab_vpf(df_taxa_vpf, fin_ab_vpfcount, fin_ab_vpfcov, colnames(df_collapse_mock), vlevel, mincov, remove_unclassified=remove_unclassified)
    df_collapse_brackenVPF <- create_ab_brackenVPF(df_taxa_vpf, fin_ab_brackenVPF, colnames(df_collapse_mock), vlevel, remove_unclassified=remove_unclassified)
    df_collapse_brackenSTD <- create_ab_brackenSTD(fin_t2lin_krakenDB, fin_ab_brackenSTD, colnames(df_collapse_mock), vlevel, remove_unclassified=remove_unclassified)
    calc_confmat(df_collapse_mock, df_collapse_vpf, df_collapse_brackenVPF, df_collapse_brackenSTD, vlevel, minrel)
  })

  df <- bind_rows(dfs) %>% mutate(taxalevel = vlevel)
  dfmetrics[[vlevel]] <- df
}

vpf_prab <- do.call(rbind, dfmetrics) %>%
  data.frame() %>%
  mutate(value = 100 * value) %>%
  mutate(taxalevel = factor(taxalevel, levels = vlevels))

usethis::use_data(vpf_prab, overwrite = TRUE)
