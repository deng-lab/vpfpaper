library(tidyverse)
library(data.table)
library(ggrepel)

# -------------------- functions ------------------
read_CheckV <- function(fpath) {
  dtbl <- fread(fpath) %>%
    mutate(checkv_quality = factor(checkv_quality, levels=c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"))) %>%
    column_to_rownames("contig_id") %>%
    setnames(colnames(.), paste0("checkv_", colnames(.))) %>%
    setnames("checkv_checkv_quality", "checkv_quality") %>%
    rownames_to_column("Contig")
  return(dtbl)
}


read_VirSorter2 <- function(fpath) {
  dtbl <- fread(fpath) %>%
    column_to_rownames("seqname") %>%
    setnames(colnames(.), paste0("virsorter2_", colnames(.))) %>%
    rownames_to_column("Contig")
  return(dtbl)
}


read_CATBAT <- function(fpath) {
  dtbl <- fread(fpath, fill = TRUE, sep = "\t") %>%
    column_to_rownames("# contig") %>%
    setnames("lineage scores", "lineage_score") %>%
    setnames(colnames(.), paste0("CATBAT_", colnames(.))) %>%
    rownames_to_column("virsorter2_contig_id")
  return(dtbl)
}


read_vConTACT2 <- function(fpath, assembler_label="_NODE_") {
  # TODO: need to tested using a test dataset
  df_vcontact <- fread(fpath) %>%
    mutate(source=ifelse(str_detect(Genome, assembler_label), "queryseq", "refseq")) %>%
    mutate(cluster_status=ifelse(str_detect(`VC Status`, "Overlap"), "Overlap", `VC Status`)) %>%
    mutate(cluster_status=factor(cluster_status)) %>%
    setnames(colnames(.), paste0("vConTACT_", str_replace_all(colnames(.), " ", "_")))

  vc2_refs <- df_vcontact %>%
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label, negate = TRUE)) %>%
    dplyr::filter(vConTACT_VC != "")

  # Get VC stats
  vcontact_stats <- df_vcontact %>%
    group_by(vConTACT_source, vConTACT_cluster_status) %>%
    dplyr::summarise(seqs_with_VC = sum(vConTACT_VC != ""),
                     seqs_without_VC = sum(vConTACT_VC == "")) %>%
    gather(seq_status, num_seqs, -vConTACT_source, -vConTACT_cluster_status)

  plot_vc_stats <- ggplot(vcontact_stats, aes(y=num_seqs, x=seq_status, fill=vConTACT_source, label=num_seqs)) +
    geom_point(aes(color=vConTACT_source), alpha=0.5, size=3) +
    facet_grid(vConTACT_cluster_status ~ .) +
    # theme(text = element_text(size = 12)) +
    geom_text_repel()

  # Annotate clusters using reference genomes
  vc2_contigs_vclst_anno <- df_vcontact %>%
    # Get distinct cluster ID of contigs in samples
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label)) %>%
    dplyr::filter(vConTACT_VC != "") %>%
    dplyr::select(vConTACT_VC) %>%
    distinct() %>%
    # Annotate cluster with reference taxonomy
    inner_join(vc2_refs, by = "vConTACT_VC")

  # Whether contig clusters were annotated by reference (1) or not (0)
  df_vcontact <- df_vcontact %>%
    mutate(vConTACT_classified=ifelse(vConTACT_VC %in% vc2_contigs_vclst_anno$vConTACT_VC, 1, 0)) %>%
    # only choose assemblies
    dplyr::filter(str_detect(vConTACT_Genome, assembler_label)) %>%
    mutate(vConTACT_VC2=ifelse(vConTACT_VC_Status %in% c("Outlier", "Singleton"), vConTACT_Genome, vConTACT_VC)) %>%
    mutate(vConTACT_VC2=ifelse(str_detect(vConTACT_VC_Status, "Overlap"), vConTACT_VC_Status, vConTACT_VC2)) %>%
    setnames("vConTACT_Genome", "virsorter2_contig_id") %>%
    mutate(Contig=str_replace(virsorter2_contig_id, "-cat_[1-6]", "")) %>%
    mutate(virsorter2_category=as.integer(str_replace(virsorter2_contig_id, ".*-cat_", "")))

  vc2 <- list("vc_tbl" = df_vcontact,
              "vc_stats" = vcontact_stats,
              "vc_plot" = plot_vc_stats,
              "vc_annotated"= vc2_contigs_vclst_anno,
              "vc_refs" = vc2_refs)
  return(vc2)
}


read_vtaxonomy <- function(fpath) {
  dtbl <- fread(fpath) %>%
    setnames(colnames(.), paste0("vtaxa_", colnames(.))) %>%
    setnames("vtaxa_contig_id", "virsorter2_contig_id")

  return(dtbl)
}

read_vtaxonomy2 <- function(fpath) {
  dtbl <- fread(fpath) %>%
    mutate(Contig=str_replace(contig_id, "-cat_[1-6]", "")) %>%
    dplyr::select(-length) %>%
    setnames(c("contig_id", "Superkingdom"), c("virsorter2_contig_id", "Kingdom")) %>%
    mutate(virsorter2_category=as.integer(str_replace(virsorter2_contig_id, ".*-cat_", "")))
  return(dtbl)
}


read_lifestyle <- function(fpath) {
  dtbl <- fread(fpath) %>%
    column_to_rownames("V1") %>%
    setnames(colnames(.), paste0("lifestyle_", colnames(.))) %>%
    mutate(lifestyle=ifelse(lifestyle_Temperate>=lifestyle_Virulent, "temperate", "virulent")) %>%
    rownames_to_column("virsorter2_contig_id")
  return(dtbl)
}


read_virhost <- function(fpath) {
  dtbl <- fread(fpath) %>%
    setnames(colnames(.), paste0("vhost_", colnames(.))) %>%
    setnames("vhost_contig_id", "virsorter2_contig_id")
  return(dtbl)
}


read_abundance <- function(fpath, anno_viral, sample_ids) {
  # for Phyloseq
  df_abundance <- fread(fpath) %>%
    inner_join(anno_viral[c("Contig", "checkv_contig_length")]) %>%
    arrange(desc(checkv_contig_length), Contig) %>%
    dplyr::select(-checkv_contig_length) %>%
    column_to_rownames("Contig") %>%
    dplyr::select(sample_ids) %>%
    rownames_to_column("Contig")
  return(df_abundance)
}

read_abundance2 <- function(fpath, sample_metadata, rownames="Contig", covfrac_threshold=0.2, ...) {
  # for TreeSummarizedExperiment
  myArgs <- list(...)
  if (!is.null(myArgs[['df_covfrac']])) {
    df_filter <- df_covfrac
    df_filter[df_filter<covfrac_threshold] <- 0
    df_filter[df_filter>=covfrac_threshold] <- 1
  } else {
    df_filter <- 1
  }

  df <- fread(fpath) %>%
    column_to_rownames(rownames)
  smeta_vir <- sample_metadata %>%
    dplyr::filter(id_virome %in% colnames(df))
  df2 <- df %>%
    setnames(smeta_vir$id_virome, smeta_vir$sample_id) %>%
    dplyr::select(smeta_vir$sample_id)
  return(as.matrix(df2 * df_filter))
}

read_abundance3 <- function(fpath, sample_metadata, rownames="Contig", id_col = "id_virome") {
  df <- fread(fpath) %>%
    column_to_rownames(rownames)
  smeta_vir <- sample_metadata %>%
    dplyr::filter(.data[[id_col]] %in% colnames(df))
  df2 <- df %>%
    setnames(smeta_vir[[id_col]], smeta_vir$sample_id) %>%
    dplyr::select(smeta_vir$sample_id) %>%
    as.matrix()
  return(df2)
}

create_viral_abundance <- function(df_abundance, df_covfrac, covfrac_threshold=0.2) {
  df_abundance <- df_abundance %>%
    column_to_rownames("Contig")
  df_filter <- df_covfrac %>%
    column_to_rownames("Contig")

  df_filter[df_filter<covfrac_threshold] <- 0
  df_filter[df_filter>=covfrac_threshold] <- 1
  df_abundance_filtered <- df_abundance*df_filter
  return(df_abundance_filtered)
}


create_vir_phyloseq <- function(smeta, df_abundance, contig_anno) {
  smeta_sel_vir <- smeta %>%
    dplyr::filter(id_virome %in% colnames(df_abundance))
  phy_abundance <- df_abundance %>%
    setnames(smeta_sel_vir$id_virome, smeta_sel_vir$sample_id) %>%
    as.matrix() %>%
    otu_table(taxa_are_rows = TRUE)

  phy_anno <- contig_anno %>%
    mutate(ctgid = Contig) %>%
    column_to_rownames("ctgid") %>%
    as.matrix() %>%
    tax_table()

  phy_smeta <- smeta %>%
    dplyr::filter(sample_id %in% colnames(df_abundance)) %>%
    column_to_rownames("sample_id") %>%
    dplyr::select(-c("id_virome", "id_16s", "id_18s")) %>%
    sample_data()

  vir_pseq <- phyloseq(phy_abundance, phy_anno, phy_smeta)
  return(vir_pseq)
}


read_gene_annotation <- function(fpath, contig_anno) {
  df <- fread(fpath) %>%
    dplyr::select(-c(fasta)) %>%
    setnames(c("scaffold"), c("virsorter2_contig_id")) %>%
    left_join(contig_anno[c("Contig", "virsorter2_contig_id")], by = "virsorter2_contig_id") %>%
    filter(!is.na(Contig))
  return(df)
}


tax_glom_by_column <- function(pseq, column, bad_empty = c(NA, "", " ")) {
  df <- data.frame(pseq@tax_table) %>%
    rownames_to_column("feature_id") %>%
    mutate(cluster_id = ifelse(.data[[column]] %in% bad_empty, feature_id, .data[[column]]))
  cluster_id <- setNames(df[["cluster_id"]], df[["feature_id"]])
  pseq2 <- merge_taxa_vec(pseq, group = cluster_id, tax_adjust = 0)
  return(pseq2)
}

clean_taxa <- function(pseq) {
  tax.clean <- data.frame(pseq@tax_table)
  tax.clean[is.na(tax.clean)] <- ""
  for (i in 1:nrow(tax.clean)){
    if (tax.clean[i,7] != ""){
      tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
    } else if (tax.clean[i,2] == ""){
      kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
      tax.clean[i, 2:7] <- kingdom
    } else if (tax.clean[i,3] == ""){
      phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
      tax.clean[i, 3:7] <- phylum
    } else if (tax.clean[i,4] == ""){
      class <- paste("Unclassified", tax.clean[i,3], sep = " ")
      tax.clean[i, 4:7] <- class
    } else if (tax.clean[i,5] == ""){
      order <- paste("Unclassified", tax.clean[i,4], sep = " ")
      tax.clean[i, 5:7] <- order
    } else if (tax.clean[i,6] == ""){
      family <- paste("Unclassified", tax.clean[i,5], sep = " ")
      tax.clean[i, 6:7] <- family
    } else if (tax.clean[i,7] == ""){
      tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
    }
  }
  pseq@tax_table <- as.matrix(tax.clean) %>% tax_table()
  return(pseq)
}


create_sampleMap <- function(tse_obj, assay_name) {
  df <- data.frame(assay = assay_name,
                   primary = colnames(tse_obj),
                   colname = colnames(tse_obj))
  return(df)
}


abundance_filter_by_coverage <- function(tse, covfrac_thr=0.2, abundances=c("counts", "tpm_raw", "trimmed_mean")) {
  df_covfrac <- assay(tse, "covfrac")
  df_covfrac[df_covfrac<covfrac_thr] <- 0
  df_covfrac[df_covfrac>=covfrac_thr] <- 1
  for (abundance in abundances) {
    assay(tse, abundance) <- assay(tse, abundance) * df_covfrac
  }
  return(tse)
}


filter_feature_by_prevalance <- function(tse_obj, min_rel, min_prev) {
  ft_prev <- getPrevalence(tse_obj, detection = min_rel, sort = T, as_relative = T, include_lowest = T) * ncol(tse_obj)
  tse_sel <- tse_obj[ft_prev >= min_prev,]
  return(tse_sel)
}



get_alpha_diversity <- function(tse) {
  # tse <- tse_obj
  tse <- transformCounts(tse, method = "relabundance")

  # colData(tse) <- colData(tse) %>%
  #   data.frame() %>%
  #   dplyr::select(colnames(.)[str_detect(colnames(.), "richness_|diversity_|even_|dominance_|divergence", negate = T)]) %>%
  #   DataFrame()

  # Richness
  index_richness <- c("observed", "chao1", "ace", "hill")
  tse <- estimateRichness(tse,
                          abund_values = "counts",
                          index = index_richness,
                          name = paste0("richness_", index_richness))

  # Diversity
  index_diversity = c("coverage", "fisher", "gini_simpson", "inverse_simpson", "shannon", "log_modulo_skewness")
  tse <- estimateDiversity(tse,
                           abund_values = "counts",
                           index = index_diversity,
                           name = paste0("diversity_", index_diversity))

  # Evenness
  index_even <- c("pielou", "camargo", "simpson_evenness", "evar", "bulla")
  tse <- estimateEvenness(tse,
                          abund_values = "counts",
                          index = index_even,
                          name = paste0("even_", index_even))

  # Dominance
  index_dominance <- c("absolute", "dbp", "core_abundance", "gini", "dmn", "relative", "simpson_lambda")
  tse <- estimateDominance(tse,
                           abund_values = "counts",
                           index = index_dominance,
                           name = paste0("dominance_", index_dominance))


  # Divergence
  tse <- estimateDivergence(tse,
                            abund_values = "counts",
                            reference = "median",
                            FUN = vegan::vegdist,
                            method = "bray",
                            name = "divergence")

  coldata <- colData(tse) %>% data.frame() %>% rownames_to_column("sample_id")
  return(coldata)
}
