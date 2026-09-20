## =========================================================
## GO enrichment + global semantic module pipeline
## version: filter generic terms + filter oversized GO terms
## =========================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GOSemSim)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

## -----------------------------
## user parameters
## -----------------------------

## highly generic GO terms to remove before semantic clustering
generic_terms <- c(
  "cell fate commitment",
  "cell differentiation",
  "cell development",
  "developmental process",
  "multicellular organism development",
  "anatomical structure development",
  "system development",
  "cellular developmental process",
  "regulation of developmental process",
  "regulation of cell differentiation"
)

## remove GO terms with very large background annotation size
## you can try 300 / 500 / 800 depending on strictness
max_term_size <- 200

## semantic similarity cutoffs
sim_cut_BP <- 0.5
sim_cut_MF <- 0.5
sim_cut_CC <- 0.5

## -----------------------------
## helper functions
## -----------------------------

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

parse_ratio_num <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) < 1) return(NA_real_)
    suppressWarnings(as.numeric(z[1]))
  })
}

parse_ratio_den <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) < 2) return(NA_real_)
    suppressWarnings(as.numeric(z[2]))
  })
}

## -----------------------------
## GO enrichment for one cluster
## -----------------------------

cluP.go <- function(input.tg, ref.tg, keyType = "SYMBOL", ont,
                    qvalue = 0.05,
                    pvalue = 0.05,
                    minGSSize = 10,
                    maxGSSize = 300) {

  input.tg <- unique(input.tg)
  ref.tg <- unique(ref.tg)

  if (length(input.tg) < minGSSize) {
    return(data.frame(
      ID = NA_character_,
      Description = NA_character_,
      Ontology = ont,
      qvalue = NA_real_,
      Count = NA_real_,
      BgRatio = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  ego <- tryCatch({
    clusterProfiler::enrichGO(
      gene = input.tg,
      universe = ref.tg,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = keyType,
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue,
      qvalueCutoff = qvalue,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      readable = TRUE
    )
  }, error = function(e) {
    NULL
  })

  if (is.null(ego)) {
    return(data.frame(
      ID = NA_character_,
      Description = NA_character_,
      Ontology = ont,
      qvalue = NA_real_,
      Count = NA_real_,
      BgRatio = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  ego.df <- as.data.frame(ego)

  if (nrow(ego.df) == 0) {
    return(data.frame(
      ID = NA_character_,
      Description = NA_character_,
      Ontology = ont,
      qvalue = NA_real_,
      Count = NA_real_,
      BgRatio = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    ID = ego.df$ID,
    Description = ego.df$Description,
    Ontology = ont,
    qvalue = ego.df$qvalue,
    Count = ego.df$Count,
    BgRatio = ego.df$BgRatio,
    stringsAsFactors = FALSE
  )
}

## -----------------------------
## run GO for all clusters
## -----------------------------

run_go_for_clusters <- function(direction = c("up", "down"),
                                res1,
                                ref,
                                roi.proximal.tg,
                                tg.keeped,
                                cluster_ids,
                                keyType = "SYMBOL",
                                q_cutoff = 0.05,
                                p_cutoff = 0.05,
                                minGSSize = 10,
                                maxGSSize = 300) {
  direction <- match.arg(direction)

  if (direction == "up") {
    target_gene_id <- rownames(res1[res1$change == "Up", , drop = FALSE])
  } else {
    target_gene_id <- rownames(res1[res1$change == "Down", , drop = FALSE])
  }

  target_gene_name <- unique(ref[ref$gene_id %in% target_gene_id, "gene_name"])
  universe <- unique(ref[ref$gene_id %in% tg.keeped, "gene_name"])

  all.go <- data.frame()

  for (clu.id in cluster_ids) {
    message("Running ", direction, " GO for ", clu.id, "...")

    clu.tg.id <- intersect(
      roi.proximal.tg[roi.proximal.tg$cluster == clu.id, "id"],
      tg.keeped
    )

    clu.tg <- unique(ref[ref$gene_id %in% clu.tg.id, "gene_name"])
    clu.target.tg <- intersect(clu.tg, target_gene_name)

    out_bp <- cluP.go(
      input.tg = clu.target.tg,
      ref.tg = universe,
      keyType = keyType,
      ont = "BP",
      qvalue = q_cutoff,
      pvalue = p_cutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )

    out_mf <- cluP.go(
      input.tg = clu.target.tg,
      ref.tg = universe,
      keyType = keyType,
      ont = "MF",
      qvalue = q_cutoff,
      pvalue = p_cutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )

    out_cc <- cluP.go(
      input.tg = clu.target.tg,
      ref.tg = universe,
      keyType = keyType,
      ont = "CC",
      qvalue = q_cutoff,
      pvalue = p_cutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )

    output <- dplyr::bind_rows(out_bp, out_mf, out_cc) %>%
      dplyr::mutate(cluster = clu.id)

    all.go <- dplyr::bind_rows(all.go, output)
  }

  all.go
}

## -----------------------------
## clean and filter GO table
## -----------------------------

prepare_go_df <- function(all.go.df,
                          generic_terms = NULL,
                          max_term_size = Inf) {

  go_df <- all.go.df %>%
    dplyr::filter(!is.na(ID), !is.na(Description), !is.na(qvalue)) %>%
    dplyr::distinct(cluster, ID, Description, Ontology, .keep_all = TRUE) %>%
    dplyr::mutate(
      score = -log10(qvalue),
      term_size = parse_ratio_num(BgRatio),
      bg_universe_size = parse_ratio_den(BgRatio)
    )

  if (!is.null(generic_terms) && length(generic_terms) > 0) {
    go_df <- go_df %>%
      dplyr::filter(!Description %in% generic_terms)
  }

  if (is.finite(max_term_size)) {
    go_df <- go_df %>%
      dplyr::filter(is.na(term_size) | term_size <= max_term_size)
  }

  go_df
}

## -----------------------------
## build GO semantic modules
## -----------------------------

build_go_modules <- function(df, ont = "BP", sim_cut = 0.5) {
  subdf <- df %>%
    dplyr::filter(Ontology == ont) %>%
    dplyr::distinct(ID, Description, Ontology)

  if (nrow(subdf) == 0) {
    return(data.frame())
  }

  if (nrow(subdf) == 1) {
    return(
      subdf %>%
        dplyr::mutate(module = paste0(ont, "_M1")) %>%
        dplyr::select(ID, Description, Ontology, module)
    )
  }

  semData <- GOSemSim::godata(
    annoDb = "org.Hs.eg.db",
    ont = ont
  )

  sim_mat <- GOSemSim::mgoSim(
    GO1 = subdf$ID,
    GO2 = subdf$ID,
    semData = semData,
    measure = "Wang",
    combine = NULL
  )

  sim_mat <- as.matrix(sim_mat)
  rownames(sim_mat) <- subdf$ID
  colnames(sim_mat) <- subdf$ID
  sim_mat[is.na(sim_mat)] <- 0
  diag(sim_mat) <- 1

  dist_mat <- stats::as.dist(1 - sim_mat)
  hc <- stats::hclust(dist_mat, method = "average")
  module_id <- stats::cutree(hc, h = 1 - sim_cut)

  data.frame(
    ID = names(module_id),
    module = paste0(ont, "_M", module_id),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(subdf, by = "ID") %>%
    dplyr::select(ID, Description, Ontology, module)
}

## -----------------------------
## choose representative term per module
## -----------------------------

pick_module_representative <- function(go_df, term2module_df) {
  go_df %>%
    dplyr::left_join(term2module_df, by = c("ID", "Description", "Ontology")) %>%
    dplyr::filter(!is.na(module)) %>%
    dplyr::group_by(module, ID, Description, Ontology) %>%
    dplyr::summarise(
      n_cluster = dplyr::n_distinct(cluster),
      mean_score = safe_mean(score),
      max_score = safe_max(score),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      module,
      dplyr::desc(n_cluster),
      dplyr::desc(mean_score),
      dplyr::desc(max_score),
      Description
    ) %>%
    dplyr::group_by(module) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(module_label = Description)
}

## -----------------------------
## build matrices
## -----------------------------

build_term_level_matrix <- function(go_df) {
  go_df %>%
    dplyr::select(ID, Description, Ontology, cluster, score) %>%
    tidyr::pivot_wider(
      id_cols = c(ID, Description, Ontology),
      names_from = cluster,
      values_from = score
    ) %>%
    dplyr::arrange(Ontology, Description)
}

build_module_matrix <- function(go_df, term2module_df, module_label_df) {
  go_df %>%
    dplyr::left_join(term2module_df, by = c("ID", "Description", "Ontology")) %>%
    dplyr::filter(!is.na(module)) %>%
    dplyr::group_by(module, cluster) %>%
    dplyr::summarise(score = safe_max(score), .groups = "drop") %>%
    dplyr::left_join(
      module_label_df %>% dplyr::select(module, module_label, Ontology),
      by = "module"
    ) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      id_cols = c(module, module_label, Ontology),
      names_from = cluster,
      values_from = score
    ) %>%
    dplyr::arrange(Ontology, module_label)
}

build_module_members <- function(term2module_df) {
  term2module_df %>%
    dplyr::group_by(module, Ontology) %>%
    dplyr::summarise(
      member_terms = paste(Description, collapse = " ; "),
      member_ids = paste(ID, collapse = " ; "),
      n_terms = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(Ontology, module)
}

get_module_presence <- function(module_mat_df, cluster_ids) {
  module_mat_df %>%
    dplyr::mutate(
      n_cluster_present = rowSums(
        !is.na(as.matrix(dplyr::select(., dplyr::any_of(cluster_ids))))
      )
    )
}

## -----------------------------
## plot heatmap
## -----------------------------

plot_module_heatmap <- function(module_mat_df, outfile, cluster_ids,
                                main = "GO module heatmap") {
  if (nrow(module_mat_df) == 0) return(NULL)

  mat <- module_mat_df %>%
    dplyr::select(dplyr::any_of(cluster_ids)) %>%
    as.data.frame()

  rownames(mat) <- paste(module_mat_df$Ontology, module_mat_df$module_label, sep = " | ")
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0

  if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)

  grDevices::pdf(outfile, width = 8, height = max(6, nrow(mat) * 0.18))
  pheatmap::pheatmap(
    mat,
    scale = "row",
    clustering_method = "ward.D2",
    border_color = NA,
    fontsize_row = 7,
    fontsize_col = 10,
    main = main
  )
  grDevices::dev.off()
}

## -----------------------------
## full pipeline
## -----------------------------

run_global_module_pipeline <- function(all.go.df,
                                       sim_cut_BP = 0.5,
                                       sim_cut_MF = 0.5,
                                       sim_cut_CC = 0.5,
                                       generic_terms = NULL,
                                       max_term_size = Inf) {
  go_df <- prepare_go_df(
    all.go.df = all.go.df,
    generic_terms = generic_terms,
    max_term_size = max_term_size
  )

  term_mat <- build_term_level_matrix(go_df)

  term2module_bp <- build_go_modules(go_df, ont = "BP", sim_cut = sim_cut_BP)
  term2module_mf <- build_go_modules(go_df, ont = "MF", sim_cut = sim_cut_MF)
  term2module_cc <- build_go_modules(go_df, ont = "CC", sim_cut = sim_cut_CC)

  term2module <- dplyr::bind_rows(term2module_bp, term2module_mf, term2module_cc)

  module_label <- pick_module_representative(go_df, term2module)

  module_mat <- build_module_matrix(go_df, term2module, module_label)

  module_members <- build_module_members(term2module)

  list(
    go_df = go_df,
    term_mat = term_mat,
    term2module = term2module,
    module_label = module_label,
    module_mat = module_mat,
    module_members = module_members
  )
}

############################
## 1. Parameters
############################
outdir <- "GO_global_module_results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

keyType_use <- "SYMBOL"
q_cutoff <- 0.05
p_cutoff <- 0.05
minGSSize_use <- 10
maxGSSize_use <- 200

# GO semantic similarity cutoff
# larger -> stricter merge; smaller -> more aggressive merge
sim_cut_BP <- 0.5
sim_cut_MF <- 0.5
sim_cut_CC <- 0.5

cluster_ids <- c(1:7)


all.go.up <- run_go_for_clusters(
  direction = "up",
  res1 = res1,
  ref = ref,
  roi.proximal.tg = roi.proximal.tg,
  tg.keeped = tg.keeped,
  cluster_ids = cluster_ids,
  keyType = keyType_use,
  q_cutoff = q_cutoff,
  p_cutoff = p_cutoff,
  minGSSize = minGSSize_use,
  maxGSSize = maxGSSize_use
)

all.go.down <- run_go_for_clusters(
  direction = "down",
  res1 = res1,
  ref = ref,
  roi.proximal.tg = roi.proximal.tg,
  tg.keeped = tg.keeped,
  cluster_ids = cluster_ids,
  keyType = keyType_use,
  q_cutoff = q_cutoff,
  p_cutoff = p_cutoff,
  minGSSize = minGSSize_use,
  maxGSSize = maxGSSize_use
)

## run pipeline for up-regulated genes
res.up <- run_global_module_pipeline(
  all.go.df = all.go.up,
  sim_cut_BP = sim_cut_BP,
  sim_cut_MF = sim_cut_MF,
  sim_cut_CC = sim_cut_CC,
  generic_terms = generic_terms,
  max_term_size = max_term_size
)

## run pipeline for down-regulated genes
res.down <- run_global_module_pipeline(
  all.go.df = all.go.down,
  sim_cut_BP = sim_cut_BP,
  sim_cut_MF = sim_cut_MF,
  sim_cut_CC = sim_cut_CC,
  generic_terms = generic_terms,
  max_term_size = max_term_size
)


write.csv(res.up$go_df, "GO/d7_d0_res_up_go_df_filtered.csv", row.names = FALSE)
write.csv(res.up$term2module, "GO/d7_d0_res_up_term2module.csv", row.names = FALSE)
write.csv(res.up$module_mat, "GO/d7_d0_res_up_module_mat.csv", row.names = FALSE)
write.csv(res.up$module_members, "GO/d7_d0_res_up_module_members.csv", row.names = FALSE)

write.csv(res.down$go_df, "GO/d7_d0_res_down_go_df_filtered.csv", row.names = FALSE)
write.csv(res.down$term2module, "GO/d7_d0_res_down_term2module.csv", row.names = FALSE)
write.csv(res.down$module_mat, "GO/d7_d0_res_down_module_mat.csv", row.names = FALSE)
write.csv(res.down$module_members, "GO/d7_d0_res_down_module_members.csv", row.names = FALSE)




all.go.up <- run_go_for_clusters(
  direction = "up",
  res1 = res2,
  ref = ref,
  roi.proximal.tg = roi.proximal.tg,
  tg.keeped = tg.keeped,
  cluster_ids = cluster_ids,
  keyType = keyType_use,
  q_cutoff = q_cutoff,
  p_cutoff = p_cutoff,
  minGSSize = minGSSize_use,
  maxGSSize = maxGSSize_use
)

all.go.down <- run_go_for_clusters(
  direction = "down",
  res1 = res2,
  ref = ref,
  roi.proximal.tg = roi.proximal.tg,
  tg.keeped = tg.keeped,
  cluster_ids = cluster_ids,
  keyType = keyType_use,
  q_cutoff = q_cutoff,
  p_cutoff = p_cutoff,
  minGSSize = minGSSize_use,
  maxGSSize = maxGSSize_use
)

## run pipeline for up-regulated genes
res.up <- run_global_module_pipeline(
  all.go.df = all.go.up,
  sim_cut_BP = sim_cut_BP,
  sim_cut_MF = sim_cut_MF,
  sim_cut_CC = sim_cut_CC,
  generic_terms = generic_terms,
  max_term_size = max_term_size
)

## run pipeline for down-regulated genes
res.down <- run_global_module_pipeline(
  all.go.df = all.go.down,
  sim_cut_BP = sim_cut_BP,
  sim_cut_MF = sim_cut_MF,
  sim_cut_CC = sim_cut_CC,
  generic_terms = generic_terms,
  max_term_size = max_term_size
)


write.csv(res.up$go_df, "GO/d14-d7_res_up_go_df_filtered.csv", row.names = FALSE)
write.csv(res.up$term2module, "GO/d14-d7_res_up_term2module.csv", row.names = FALSE)
write.csv(res.up$module_mat, "GO/d14-d7_res_up_module_mat.csv", row.names = FALSE)
write.csv(res.up$module_members, "GO/d14-d7_res_up_module_members.csv", row.names = FALSE)

write.csv(res.down$go_df, "GO/d14-d7_res_down_go_df_filtered.csv", row.names = FALSE)
write.csv(res.down$term2module, "GO/d14-d7_res_down_term2module.csv", row.names = FALSE)
write.csv(res.down$module_mat, "GO/d14-d7_res_down_module_mat.csv", row.names = FALSE)
write.csv(res.down$module_members, "GO/d14-d7_res_down_module_members.csv", row.names = FALSE)



selected_modules1 <- c(
  "BP_M2",
  "BP_M16",
  "BP_M10",
  "BP_M97",
  "BP_M62",
  "BP_M8",
  "BP_M31",
  "BP_M102",
  "MF_M11",
  "MF_M10",
  "BP_M105",
  "BP_M85",
  "BP_M89",
  "BP_M106",
  "BP_M43",
  "BP_M54",
  "BP_M100",
  "BP_M101",
  "BP_M27",
  "BP_M57",
  "BP_M92",
  "BP_M95",
  "BP_M41",
  "BP_M38",
  "MF_M5",
  "BP_M21",
  "BP_M80",
  "MF_M7",
  "BP_M46",
  "BP_M84"
)

d0_d7_up <- read.table(file = "GO/d7_d0_res_up_module_mat.csv",sep = ",",header = T)

out1 <- d0_d7_up[d0_d7_up$module %in% selected_modules1,]


selected_modules2 <- c(
  "BP_M61",
  "BP_M9",
  "BP_M31",
  "BP_M14",
  "BP_M81",
  "BP_M42",
  "BP_M79",
  "MF_M8",
  "MF_M3",
  "CC_M10",
  "BP_M4",
  "BP_M7",
  "BP_M76",
  "BP_M50",
  "BP_M64",
  "BP_M44",
  "BP_M59",
  "BP_M25",
  "BP_M13",
  "BP_M68",
  "BP_M78",
  "BP_M82",
  "BP_M5",
  "BP_M15",
  "BP_M54",
  "BP_M33",
  "BP_M55",
  "BP_M63",
  "BP_M36",
  "MF_M2"
)

d7_d14_up <- read.table(file = "GO/d14-d7_res_up_module_mat.csv",sep = ",",header = T)

out2 <- d7_d14_up[d7_d14_up$module %in% selected_modules2,]


intersect(d0_d7_up$module_label,d7_d14_up$module_label)


d0_d7_down <- read.table(file = "GO/d7_d0_res_down_module_mat.csv",sep = ",",header = T)
selected_modules3 <- c(
  "BP_M10",  # somatic stem cell population maintenance
  "BP_M29",  # developmental induction
  "BP_M35",  # regeneration
  "BP_M30",  # positive regulation of animal organ morphogenesis
  "BP_M32",  # osteoblast differentiation
  "BP_M18",  # neural precursor cell proliferation
  "BP_M24",  # regulation of neural precursor cell proliferation
  "BP_M20",  # oligodendrocyte differentiation
  "BP_M22",  # negative regulation of neurogenesis
  "BP_M36",  # columnar/cuboidal epithelial cell differentiation
  "BP_M14",  # regulation of morphogenesis of an epithelium
  "BP_M37",  # digestive tract development
  "BP_M31",  # cardiac muscle tissue growth
  "BP_M15",  # tissue migration
  "BP_M26",  # actomyosin structure organization
  "BP_M33",  # cell-cell adhesion mediated by cadherin
  "MF_M1",   # extracellular matrix binding
  "CC_M5",   # catenin complex
  "CC_M7",   # lamellipodium
  "CC_M3",   # site of polarized growth
  "BP_M17",  # regulation of Notch signaling pathway
  "BP_M27",  # regulation of JNK cascade
  "BP_M1",   # ephrin receptor signaling pathway
  "MF_M5",   # transmembrane receptor protein serine/threonine kinase binding
  "MF_M2",   # protein tyrosine kinase activity
  "MF_M6",   # phosphoprotein binding
  "MF_M4",   # G protein-coupled receptor binding
  "BP_M11",  # calcium ion homeostasis
  "BP_M12",  # negative regulation of inflammatory response
  "BP_M9"    # phagocytosis
)


out3 <- d0_d7_down[d0_d7_down$module %in% selected_modules3,]

write.table(out3,file = "GO/d0_d7_down.mat2.tsv",row.names = F,col.names = T,sep = "\t",quote = F)
