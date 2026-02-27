# ===========================
# Strategy B: Infer sex from expression (XIST vs Y genes) and run pseudobulk DE
# Dataset: GSE248728 (scRNA-seq processed MTX/TSV per capture)
# ===========================

suppressPackageStartupMessages({
  library(GEOquery)
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(patchwork)
})

# ---------------------------
# 0) Paths and settings
# ---------------------------
gse_id <- "GSE248728"
out_dir <- "data/raw_geo"
geo_dir <- file.path(out_dir, gse_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)


# QC thresholds (tune after inspecting violin plots)
qc_min_features <- 300
qc_max_features <- 6000
qc_max_percent_mt <- 20

# Clustering parameters
pca_dims <- 1:30
cluster_resolution <- 0.6

# Sex inference genes
# female_marker <- "XIST"
# y_genes <- c("RPS4Y1", "DDX3Y", "KDM5D", "EIF1AY")

female_marker <- "XIST"
y_genes <- c("RPS4Y1","DDX3Y","KDM5D","EIF1AY","USP9Y","UTY","ZFY","NLGN4Y","TMSB4Y","PRKY")



# Sex inference threshold on log2CPM difference (conservative default)
sex_delta <- 1.0

# ---------------------------
# 1) Download GEO series matrix + supplementary files
# ---------------------------
message("Downloading GEO series matrix and supplementary files ...")
gse <- getGEO(gse_id, GSEMatrix = TRUE)[[1]]
pheno <- pData(gse)
write.csv(pheno, file.path(out_dir, paste0(gse_id, "_pData.csv")), row.names = FALSE)

getGEOSuppFiles(GEO = gse_id, makeDirectory = TRUE, baseDir = out_dir)

# ---------------------------
# 2) Read all captures and merge into one Seurat object
# ---------------------------
message("Detecting capture prefixes ...")
caps <- list.files(geo_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = FALSE) |>
  str_replace("_matrix\\.mtx\\.gz$", "") |>
  sort()

stopifnot(length(caps) > 0)
message("Found captures: ", paste(caps, collapse = ", "))

read_capture <- function(capture_prefix, dir_path) {
  mtx_path  <- file.path(dir_path, paste0(capture_prefix, "_matrix.mtx.gz"))
  feat_path <- file.path(dir_path, paste0(capture_prefix, "_features.tsv.gz"))
  barc_path <- file.path(dir_path, paste0(capture_prefix, "_barcodes.tsv.gz"))
  
  stopifnot(file.exists(mtx_path), file.exists(feat_path), file.exists(barc_path))
  
  mtx <- Matrix::readMM(gzfile(mtx_path))
  
  # Windows-safe reading of gz TSVs
  features <- data.table::fread(feat_path, header = FALSE)
  barcodes <- data.table::fread(barc_path, header = FALSE)
  
  gene_names <- features$V2
  barcode_names <- barcodes$V1
  
  rownames(mtx) <- make.unique(gene_names)
  colnames(mtx) <- barcode_names
  
  obj <- CreateSeuratObject(counts = mtx, project = capture_prefix, min.cells = 0, min.features = 0)
  obj$capture <- capture_prefix
  obj
}

message("Reading captures and merging ...")
objs <- lapply(caps, function(cp) read_capture(cp, geo_dir))
names(objs) <- caps
bm <- Reduce(function(x, y) merge(x, y), objs)

saveRDS(bm, "results/bm_raw_merged.rds")
message("Merged object: ", ncol(bm), " cells; ", nrow(bm), " genes.")

# ---------------------------
# 3) QC, normalization, clustering
# ---------------------------
bm <- readRDS("results/bm_raw_merged.rds")
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^MT-")

p_qc <- VlnPlot(bm, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) +
  plot_annotation(title = "QC distributions (pre-filter)")
ggsave("results/plots/QC_violin_pre_filter.png", p_qc, width = 12, height = 4, dpi = 200)

p_feat <- VlnPlot(bm, features = "nFeature_RNA", pt.size = 0.1) +
  ggtitle("Number of detected genes per cell") +
  theme_classic()

p_count <- VlnPlot(bm, features = "nCount_RNA", pt.size = 0.1) +
  ggtitle("Total UMI counts per cell") +
  theme_classic()

p_mt <- VlnPlot(bm, features = "percent.mt", pt.size = 0.1) +
  ggtitle("Percentage mitochondrial transcripts") +
  theme_classic()

# Save each plot separately
ggsave("results/plots/QC_nFeature_RNA.png", p_feat, width = 15, height = 5, dpi = 300)
ggsave("results/plots/QC_nCount_RNA.png", p_count, width = 15, height = 5, dpi = 300)
ggsave("results/plots/QC_percent_mt.png", p_mt, width = 15, height = 5, dpi = 300)

bm <- subset(
  bm,
  subset = nFeature_RNA >= qc_min_features &
    nFeature_RNA <= qc_max_features &
    percent.mt <= qc_max_percent_mt
)

### QC plots after filtering ##############################

p_qc <- VlnPlot(bm, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) +
  plot_annotation(title = "QC distributions (pre-filter)")
ggsave("results/plots/QC_violin_post_filter.png", p_qc, width = 12, height = 4, dpi = 200)

p_feat <- VlnPlot(bm, features = "nFeature_RNA", pt.size = 0.1) +
  ggtitle("Number of detected genes per cell") +
  theme_classic()

p_count <- VlnPlot(bm, features = "nCount_RNA", pt.size = 0.1) +
  ggtitle("Total UMI counts per cell") +
  theme_classic()

p_mt <- VlnPlot(bm, features = "percent.mt", pt.size = 0.1) +
  ggtitle("Percentage mitochondrial transcripts") +
  theme_classic()

# Save each plot separately
ggsave("results/plots/postQC_nFeature_RNA.png", p_feat, width = 15, height = 5, dpi = 300)
ggsave("results/plots/postQC_nCount_RNA.png", p_count, width = 15, height = 5, dpi = 300)
ggsave("results/plots/postQC_percent_mt.png", p_mt, width = 15, height = 5, dpi = 300)

#########################################################################################

message("After QC: ", ncol(bm), " cells retained.")

bm <- NormalizeData(bm)
bm <- FindVariableFeatures(bm, selection.method = "vst", nfeatures = 3000)
bm <- ScaleData(bm, features = VariableFeatures(bm))
bm <- RunPCA(bm, npcs = 50)

p_elbow <- ElbowPlot(bm, ndims = 50) + ggtitle("Elbow plot")
ggsave("results/plots/PCA_elbow.png", p_elbow, width = 7, height = 5, dpi = 200)

bm <- FindNeighbors(bm, dims = pca_dims)
bm <- FindClusters(bm, resolution = cluster_resolution)
bm <- RunUMAP(bm, dims = pca_dims)

p_umap_cluster <- DimPlot(bm, group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP (clusters)")
ggsave("results/plots/UMAP_clusters.png", p_umap_cluster, width = 7, height = 5, dpi = 200)

saveRDS(bm, "results/bm_qc_clustered.rds")


# ---------------------------
# 4) Infer sex per capture from expression (pseudobulk per capture)
# ---------------------------
bm <- readRDS("results/bm_qc_clustered.rds")

DefaultAssay(bm) <- "RNA"
Layers(bm[["RNA"]])

bm <- JoinLayers(bm, assay = "RNA")
Layers(bm[["RNA"]])

counts <- GetAssayData(bm, assay = "RNA", layer = "counts")

DefaultAssay(bm) <- "RNA"
bm <- JoinLayers(bm, assay = "RNA")
counts <- GetAssayData(bm, assay = "RNA", layer = "counts")

female_marker <- "XIST"
y_genes <- c("RPS4Y1","DDX3Y","KDM5D","EIF1AY","USP9Y","UTY","ZFY","NLGN4Y","TMSB4Y","PRKY")

y_present <- intersect(rownames(counts), y_genes)
stopifnot(length(y_present) >= 2)

libsize <- Matrix::colSums(counts)

logcpm_gene <- function(g) log2((counts[g, ] / libsize) * 1e6 + 1)

xist <- if (female_marker %in% rownames(counts)) as.numeric(logcpm_gene(female_marker)) else rep(NA_real_, ncol(counts))

Ymat <- sapply(y_present, function(g) as.numeric(logcpm_gene(g)))
y_score <- rowMeans(Ymat)

y_score
# Conservative classification rules
# (tune if needed by inspecting distributions)
sex_cell <- dplyr::case_when(
  !is.na(xist) & y_score >= 2.0 & xist <= 0.5 ~ "M",
  !is.na(xist) & xist >= 2.0 & y_score <= 0.5 ~ "F",
  TRUE ~ "Ambiguous"
)

bm$sex_cell <- sex_cell

mix_wide <- bm@meta.data %>%
  dplyr::count(capture, sex_cell) %>%
  dplyr::group_by(capture) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::select(-n) %>%
  tidyr::pivot_wider(names_from = sex_cell, values_from = prop, values_fill = 0) %>%
  dplyr::ungroup()

mix_wide ## sample sex information extraction 


# Pseudobulk counts per capture = sum counts over all cells in that capture

cap <- bm$capture
cap_levels <- sort(unique(cap))
cap_levels

pb <- sapply(cap_levels, function(cn) {
  Matrix::rowSums(counts[, cap == cn, drop = FALSE])
})
pb <- as.matrix(pb)  # genes x captures


# # log2CPM transform
# libsize <- colSums(pb)
# logcpm <- log2(t(t(pb) / libsize) * 1e6 + 1)
# 
# # Score captures by XIST vs mean(Y genes)
# present_y <- intersect(rownames(logcpm), y_genes)
# has_xist <- female_marker %in% rownames(logcpm)
# 
# score_tbl <- tibble(
#   capture = colnames(logcpm),
#   XIST = if (has_xist) as.numeric(logcpm[female_marker, ]) else NA_real_,
#   Ymean = if (length(present_y) > 0) as.numeric(colMeans(logcpm[present_y, , drop = FALSE])) else NA_real_
# ) %>%
#   mutate(
#     # Conservative classification; keep ambiguous if not clearly separated
#     sex_inferred = case_when(
#       !is.na(XIST) & !is.na(Ymean) & (XIST > Ymean + sex_delta) ~ "F",
#       !is.na(XIST) & !is.na(Ymean) & (Ymean > XIST + sex_delta) ~ "M",
#       TRUE ~ "Ambiguous"
#     )
#   )
# 
# write.csv(score_tbl, "results/tables/capture_sex_inference_scores.csv", row.names = FALSE)
# print(score_tbl %>% count(sex_inferred))

# # Plot capture-level sex separation
# p_sex_scatter <- ggplot(score_tbl, aes(x = XIST, y = Ymean, label = capture, color = sex_inferred)) +
#   geom_point(size = 3) +
#   geom_text(vjust = 1.2, size = 3, show.legend = FALSE) +
#   theme_minimal() +
#   labs(
#     title = "Capture-level sex inference (log2CPM)",
#     x = "log2CPM(XIST)",
#     y = "mean log2CPM(Y-linked genes)"
#   )
# ggsave("results/plots/capture_sex_inference_scatter.png", p_sex_scatter, width = 9, height = 6, dpi = 200)
# 
# # Attach inferred sex to each cell based on capture
# bm$sex_inferred <- score_tbl$sex_inferred[match(bm$capture, score_tbl$capture)]
# saveRDS(bm, "results/bm_with_inferred_sex.rds")

# ---------------------------
# 5) Sanity-check at cell level: marker behavior by inferred sex
# ---------------------------
bm <- readRDS("results/bm_with_inferred_sex.rds")

# UMAP by inferred sex
p_umap_sex <- DimPlot(bm, group.by = "sex_inferred") + ggtitle("UMAP colored by inferred sex")
ggsave("results/plots/UMAP_inferred_sex.png", p_umap_sex, width = 7, height = 5, dpi = 200)

# Feature/Vln plots for sex markers (only for present genes)
genes_to_plot <- c(female_marker, y_genes)
genes_present <- intersect(genes_to_plot, rownames(bm))

if (length(genes_present) > 0) {
  p_feat <- FeaturePlot(bm, features = genes_present[1:min(4,length(genes_present))], ncol = 2) +
    plot_annotation(title = "Sex-marker FeaturePlots")
  ggsave("results/plots/sex_marker_featureplots.png", p_feat, width = 10, height = 8, dpi = 200)
  
  p_vln <- VlnPlot(bm, features = genes_present[1:min(4,length(genes_present))],
                   group.by = "sex_inferred", pt.size = 0.1, ncol = 2) +
    plot_annotation(title = "Sex-marker violin plots by inferred sex")
  ggsave("results/plots/sex_marker_violin.png", p_vln, width = 10, height = 8, dpi = 200)
}

# ---------------------------
# 6) Sex differential expression using pseudobulk per capture within each cluster
# ---------------------------
bm <- readRDS("results/bm_with_inferred_sex.rds")

# We use clusters as "cell types" proxy unless you have annotated labels.
celltype_col <- "seurat_clusters"

md <- bm@meta.data %>%
  mutate(celltype = .data[[celltype_col]]) %>%
  filter(!is.na(sex_inferred), sex_inferred %in% c("F","M"))

DefaultAssay(bm) <- "RNA"
bm <- JoinLayers(bm, assay = "RNA") 
# Pseudobulk per capture within each celltype

pseudobulk_counts <- function(seurat_obj, meta_df, group_var, assay = "RNA", layer = "counts") {
  DefaultAssay(seurat_obj) <- assay
  
  # Ensure single layer (safe if already joined)
  if (length(Layers(seurat_obj[[assay]])) > 1) {
    seurat_obj <- JoinLayers(seurat_obj, assay = assay)
  }
  
  counts <- GetAssayData(seurat_obj, assay = assay, layer = layer)
  
  groups <- factor(meta_df[[group_var]])
  mm <- Matrix::sparse.model.matrix(~ 0 + groups)
  pb <- counts %*% mm
  colnames(pb) <- levels(groups)
  pb
}


out_de <- "results/tables/sex_DE_pseudobulk"
dir.create(out_de, recursive = TRUE, showWarnings = FALSE)

celltypes <- sort(unique(md$celltype))
de_summary <- list()

for (ct in celltypes) {
  message("DE for cluster/celltype: ", ct)
  cells_ct <- rownames(md)[md$celltype == ct]
  if (length(cells_ct) < 200) next
  
  obj_ct <- subset(bm, cells = cells_ct)
  md_ct <- md[cells_ct, , drop = FALSE]
  
  # pseudobulk per capture (sample unit)
  pb <- pseudobulk_counts(obj_ct, md_ct, group_var = "capture")
  
  pb_meta <- md_ct %>%
    group_by(capture) %>%
    summarise(
      sex_inferred = first(sex_inferred),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_cells >= 20)  # consistent minimum-cells rule used in similar pseudobulk workflows
  
  if (nrow(pb_meta) < 4) next
  pb <- pb[, pb_meta$capture, drop = FALSE]
  
  # edgeR + voom + limma (no donor random effect possible here without donor IDs)
  dge <- DGEList(counts = pb)
  keep <- filterByExpr(dge, group = pb_meta$sex_inferred)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  
  design <- model.matrix(~ sex_inferred, data = pb_meta)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # coefficient: sex_inferredM (if F is reference)
  coef_name <- grep("^sex_inferred", colnames(fit$coefficients), value = TRUE)[1]
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  
  out_csv <- file.path(out_de, paste0("DE_sex_cluster_", ct, ".csv"))
  write.csv(tt, out_csv, row.names = TRUE)
  
  # Positive-control extraction
  pc_genes <- intersect(rownames(tt), c("XIST", y_genes))
  pc <- tt[pc_genes, c("logFC","P.Value","adj.P.Val"), drop = FALSE]
  write.csv(pc, file.path(out_de, paste0("posctrl_cluster_", ct, ".csv")))
  
  de_summary[[as.character(ct)]] <- tibble(
    celltype = as.character(ct),
    n_genes_tested = nrow(tt),
    n_fdr_0_05 = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
    n_captures_used = nrow(pb_meta),
    min_cells_per_capture = min(pb_meta$n_cells)
  )
}

de_summary <- bind_rows(de_summary)
write.csv(de_summary, "results/tables/sex_DE_summary.csv", row.names = FALSE)

# Plot summary
p_de_bar <- ggplot(de_summary, aes(x = celltype, y = n_fdr_0_05)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Sex DE genes per cluster (FDR < 0.05)", x = "Cluster", y = "Count")
ggsave("results/plots/sex_DE_counts_per_cluster.png", p_de_bar, width = 10, height = 5, dpi = 200)

message("Done. Outputs in results/plots and results/tables.")