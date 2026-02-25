library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(SeuratWrappers)
library(glue)
library(data.table)
library(ggsci)
library(patchwork)

## 批次校正相关包
library(harmony)       # Harmony
library(batchelor)     # FastMNN
library(rliger)        # LIGER (iNMF)
library(conos)         # Conos
library(sva)           # ComBat
library(STACAS)        # STACAS
library(sceasy)        # R <-> Python h5ad 转换

seu = qs::qread("/home/data/tanglei/project/prostate_altas/output/02/Seurat_Integration.qs")
seu@reductions = list()

#### ========== 关键参数设定 ========== ####
batch_var   <- "sample.ID"   # <-- 手动修改批次变量（可选 "orig.ident" / "GSE.ID" 或 metadata 中的任意列）
celltype_var <- "celltype"     # <-- 细胞类型列名，scANVI 需要
n_hvg       <- 3000            # HVG 数量
n_pcs       <- 50              # PCA 维度
n_dims      <- 1:n_pcs         # 下游统一使用 50 维用于 benchmark

output_dir <- "/home/data/tanglei/project/prostate_altas/output/02"

## !! 661k 细胞数据量大，必须解除 future 全局变量大小限制
## 否则 IntegrateLayers (CCA/RPCA 等) 会报 future.globals.maxSize 错误
options(future.globals.maxSize = Inf)

cat("批次变量:", batch_var, "\n")
cat("细胞数:", ncol(seu), "\n")
cat("批次数:", length(unique(seu[[batch_var, drop = TRUE]])), "\n")

#### ========== Part 1: 公共预处理 ========== ####
## Seurat v5: split layers -> Normalize -> HVG -> Scale -> PCA

# 按批次拆分 layer（IntegrateLayers 要求）
batch_vec <- seu[[batch_var, drop = TRUE]]
if (is.null(batch_vec)) stop(paste0("metadata 中不存在批次列: ", batch_var))
batch_vec <- as.character(batch_vec)
names(batch_vec) <- colnames(seu)
if (anyNA(batch_vec)) stop(paste0("批次列含 NA，请先处理: ", batch_var))

# 如果之前已经 split 过，先 Join 再按当前 batch 重新 split
layer_names <- Layers(seu[["RNA"]])
already_split <- any(grepl("^(counts|data)\\..+$", layer_names))
if (already_split) {
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
}
seu[["RNA"]] <- split(seu[["RNA"]], f = batch_vec)

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = n_hvg, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = n_pcs, verbose = FALSE)

# 未校正 PCA 的 UMAP（作为 baseline）
seu <- RunUMAP(seu, reduction = "pca", dims = n_dims,
               reduction.name = "umap.unintegrated", verbose = FALSE)

#### ========== Part 2a: CCA 整合 ========== ####
## 注意：CCA 在大数据集（>100k 细胞）上非常耗时，如需跳过可注释本 cell
cat(">>> CCA Integration 开始...\n")
t_cca <- system.time({
  seu <- IntegrateLayers(
    object = seu,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "cca",
    verbose = FALSE
  )
})
cat("CCA 耗时:", t_cca["elapsed"], "秒\n")

qs::qsave(seu, "/home/data/tanglei/project/prostate_altas/output/02/Seurat_Integration.qs")

#### ========== Part 2b: RPCA 整合 ========== ####
cat(">>> RPCA Integration 开始...\n")
t_rpca <- system.time({
  seu <- IntegrateLayers(
    object = seu,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "rpca",
    verbose = FALSE
  )
})
cat("RPCA 耗时:", t_rpca["elapsed"], "秒\n")

qs::qsave(seu, "/home/data/tanglei/project/prostate_altas/output/02/Seurat_Integration.qs")

#### ========== Part 2c: Harmony 整合 ========== ####
cat(">>> Harmony Integration 开始...\n")
t_harmony <- system.time({
  seu <- IntegrateLayers(
    object = seu,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = FALSE
  )
})
cat("Harmony 耗时:", t_harmony["elapsed"], "秒\n")

qs::qsave(seu, "/home/data/tanglei/project/prostate_altas/output/02/Seurat_Integration.qs")

#### ========== Part 2d: FastMNN (MNN) 整合 ========== ####
cat(">>> FastMNN Integration 开始...\n")
t_mnn <- system.time({
  seu <- IntegrateLayers(
    object = seu,
    method = FastMNNIntegration,
    orig.reduction = "pca",
    new.reduction = "mnn",
    verbose = FALSE
  )
})
cat("FastMNN 耗时:", t_mnn["elapsed"], "秒\n")

qs::qsave(seu, "/home/data/tanglei/project/prostate_altas/output/02/Seurat_Integration.qs")

#### ========== Part 3: LIGER (iNMF) 整合 ========== ####
## rliger v2 workflow: 从 Seurat 提取 counts -> 构建 liger 对象 -> iNMF -> 存回 Seurat
cat(">>> LIGER Integration 开始...\n")

## 先 JoinLayers 以获取完整 counts 矩阵
seu_tmp <- JoinLayers(seu)

t_liger <- system.time({
  ## 提取 counts 和 batch 信息
  counts_mat <- GetAssayData(seu_tmp, assay = "RNA", layer = "counts")
  batch_labels <- seu_tmp[[batch_var, drop = TRUE]]
  
  ## 按 batch 拆分 counts 为 list
  batches <- unique(batch_labels)
  raw_data_list <- lapply(setNames(batches, batches), function(b) {
    idx <- which(batch_labels == b)
    counts_mat[, idx, drop = FALSE]
  })
  
  ## 创建 liger 对象
  liger_obj <- createLiger(raw_data_list)
  
  ## 标准 LIGER 流程
  liger_obj <- rliger::normalize(liger_obj)
  liger_obj <- selectGenes(liger_obj, var.thresh = 0.1)
  liger_obj <- scaleNotCenter(liger_obj)
  
  ## 运行 iNMF（online iNMF 适用于大数据集）
  k_val <- 50
  liger_obj <- runIntegration(liger_obj, k = k_val, method = "iNMF")
  liger_obj <- quantileNorm(liger_obj)
  
  ## 提取 H 矩阵（cell embedding）
  H_norm <- getMatrix(liger_obj, "H.norm")
  H_norm <- as.matrix(H_norm)

  ## 对齐 cell 顺序（兼容 rliger 自动添加 batch 前缀）
  cells_seu <- colnames(seu)
  batches_chr <- as.character(batches)

  strip_batch_prefix <- function(x, batch_names) {
    out <- as.character(x)
    for (b in batch_names) {
      pref <- paste0(b, "_")
      hit <- startsWith(out, pref)
      out[hit] <- substring(out[hit], nchar(pref) + 1L)
    }
    out
  }

  if (!is.null(rownames(H_norm)) && all(cells_seu %in% rownames(H_norm))) {
    H_norm <- H_norm[cells_seu, , drop = FALSE]
  } else if (!is.null(colnames(H_norm)) && all(cells_seu %in% colnames(H_norm))) {
    ## 如果 cells 在列，转置成 cells x factors
    H_norm <- t(H_norm[, cells_seu, drop = FALSE])
  } else if (!is.null(rownames(H_norm)) && all(cells_seu %in% strip_batch_prefix(rownames(H_norm), batches_chr))) {
    rn_norm <- strip_batch_prefix(rownames(H_norm), batches_chr)
    idx <- match(cells_seu, rn_norm)
    H_norm <- H_norm[idx, , drop = FALSE]
    rownames(H_norm) <- cells_seu
  } else if (!is.null(colnames(H_norm)) && all(cells_seu %in% strip_batch_prefix(colnames(H_norm), batches_chr))) {
    cn_norm <- strip_batch_prefix(colnames(H_norm), batches_chr)
    idx <- match(cells_seu, cn_norm)
    H_norm <- t(H_norm[, idx, drop = FALSE])
    rownames(H_norm) <- cells_seu
  } else {
    ## 名字仍无法匹配时，按构建 raw_data_list 时的顺序重建映射
    cells_by_batch <- unlist(lapply(batches_chr, function(b) {
      colnames(counts_mat)[which(batch_labels == b)]
    }), use.names = FALSE)

    if (nrow(H_norm) == length(cells_by_batch)) {
      rownames(H_norm) <- cells_by_batch
      H_norm <- H_norm[cells_seu, , drop = FALSE]
    } else if (ncol(H_norm) == length(cells_by_batch)) {
      H_norm <- t(H_norm)
      rownames(H_norm) <- cells_by_batch
      H_norm <- H_norm[cells_seu, , drop = FALSE]
    } else {
      common_r <- if (is.null(rownames(H_norm))) NA_integer_ else length(intersect(cells_seu, rownames(H_norm)))
      common_c <- if (is.null(colnames(H_norm))) NA_integer_ else length(intersect(cells_seu, colnames(H_norm)))
      common_r_strip <- if (is.null(rownames(H_norm))) NA_integer_ else length(intersect(cells_seu, strip_batch_prefix(rownames(H_norm), batches_chr)))
      common_c_strip <- if (is.null(colnames(H_norm))) NA_integer_ else length(intersect(cells_seu, strip_batch_prefix(colnames(H_norm), batches_chr)))
      stop(paste0(
        "H.norm 无法与 Seurat 细胞名对齐。\n",
        "- Seurat cells: ", length(cells_seu), "\n",
        "- H.norm dim: ", paste(dim(H_norm), collapse = " x "), "\n",
        "- H.norm rownames: ", if (is.null(rownames(H_norm))) "NULL" else "present", " (common=", common_r, ", strip_common=", common_r_strip, ")\n",
        "- H.norm colnames: ", if (is.null(colnames(H_norm))) "NULL" else "present", " (common=", common_c, ", strip_common=", common_c_strip, ")\n"
      ))
    }
  }

  colnames(H_norm) <- paste0("liger_", seq_len(ncol(H_norm)))
})

## 存入 Seurat
seu[["liger"]] <- CreateDimReducObject(
  embeddings = as.matrix(H_norm),
  key = "liger_",
  assay = "RNA"
)
cat("LIGER 耗时:", t_liger["elapsed"], "秒\n")

## 清理
rm(seu_tmp, counts_mat, batch_labels, batches, raw_data_list, liger_obj, H_norm)
gc()

qs::qsave(seu, "/home/data/tanglei/project/prostate_altas/output/02/Seurat_Integration.qs")


#### ========== Part 5: Combat 整合 (sva) ========== ####
## ComBat 在 PCA embedding 上做批次校正（比在表达矩阵上更实用且快速）
cat(">>> Combat Integration 开始...\n")

t_combat <- system.time({
  ## 提取 PCA embedding
  pca_emb <- Embeddings(seu, reduction = "pca")  # cells x n_pcs
  batch_labels <- seu[[batch_var, drop = TRUE]]
  
  ## ComBat 要求输入为 features x samples 矩阵
  pca_t <- t(pca_emb)  # n_pcs x cells
  
  ## 过滤掉只有 1 个细胞的 batch（ComBat 无法处理）
  batch_tab <- table(batch_labels)
  valid_batches <- names(batch_tab[batch_tab > 1])
  valid_idx <- which(batch_labels %in% valid_batches)
  
  combat_corrected_t <- sva::ComBat(
    dat = pca_t[, valid_idx],
    batch = batch_labels[valid_idx],
    mod = NULL,
    par.prior = TRUE
  )
  
  ## 转回 cells x n_pcs
  combat_corrected <- t(combat_corrected_t)
  
  ## 如果有被过滤的细胞，用原始 PCA 填充
  if (length(valid_idx) < ncol(seu)) {
    full_combat <- pca_emb
    full_combat[valid_idx, ] <- combat_corrected
    combat_corrected <- full_combat
  }
  
  colnames(combat_corrected) <- paste0("combat_", seq_len(ncol(combat_corrected)))
})

## 存入 Seurat
seu[["combat"]] <- CreateDimReducObject(
  embeddings = combat_corrected,
  key = "combat_",
  assay = "RNA"
)
cat("Combat 耗时:", t_combat["elapsed"], "秒\n")

rm(pca_emb, pca_t, combat_corrected_t, combat_corrected, batch_labels, batch_tab, valid_batches, valid_idx)

qs::qsave(seu, "/home/data/tanglei/project/prostate_altas/output/02/Seurat_Integration.qs")

#### ========== Part 6: STACAS 整合 ========== ####
## STACAS 专为大规模 atlas 整合设计，基于 anchor-based 但优化了锚点选择
## 注意：STACAS 在大数据集上也较慢
cat(">>> STACAS Integration 开始...\n")

t_stacas <- system.time({
  ## STACAS 需要 per-sample Seurat list（已预处理）
  seu_joined <- JoinLayers(seu)
  seu_list <- SplitObject(seu_joined, split.by = batch_var)
  rm(seu_joined); gc()
  
  ## 每个样本做预处理
  seu_list <- lapply(seu_list, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures = n_hvg, verbose = FALSE)
    return(x)
  })
  
  ## 找锚点（STACAS 的核心优势：自动过滤不可靠锚点）
  anchors <- FindAnchors.STACAS(
    seu_list,
    dims = n_dims,
    anchor.features = n_hvg
  )
  
  ## 整合
  seu_stacas <- IntegrateData(
    anchorset = anchors,
    dims = n_dims
  )
  
  ## 在 integrated assay 上做 PCA
  DefaultAssay(seu_stacas) <- "integrated"
  seu_stacas <- ScaleData(seu_stacas, verbose = FALSE)
  seu_stacas <- RunPCA(seu_stacas, npcs = n_pcs, verbose = FALSE)
  
  ## 提取 corrected PCA embedding
  stacas_emb <- Embeddings(seu_stacas, reduction = "pca")
  stacas_emb <- stacas_emb[colnames(seu), ]
  colnames(stacas_emb) <- paste0("stacas_", seq_len(ncol(stacas_emb)))
})

## 存回主对象
seu[["stacas"]] <- CreateDimReducObject(
  embeddings = stacas_emb,
  key = "stacas_",
  assay = "RNA"
)
cat("STACAS 耗时:", t_stacas["elapsed"], "秒\n")

rm(seu_list, anchors, seu_stacas, stacas_emb); gc()

qs::qsave(seu, "/home/data/tanglei/project/prostate_altas/output/02/Seurat_Integration.qs")

#### ========== Part 7: Seurat V5 转换为 h5ad (reticulate 方式) ========== ####
## 通过 reticulate 直接构建 AnnData 对象，避免 SeuratDisk 的兼容性问题
## 保留: counts、normalized data、metadata、所有 reductions、HVG 信息
## 用途: scVI / scANVI / Scanorama / BBKNN 等 Python 方法 + benchmark

library(reticulate)
use_condaenv("scvi-env", conda = "/home/data/tanglei/miniconda3/bin/conda", required = TRUE)
sc  <- import("scanpy")
sp  <- import("scipy.sparse")
np  <- import("numpy")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Part 7: Seurat V5 --> h5ad（reticulate 直接构建 AnnData）\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

## ---- 1. JoinLayers（Seurat V5 split layers 必须先合并） ---- ##
cat("[1/7] JoinLayers...\n")
seu_export <- JoinLayers(seu)

## 确定要导出的 assay
available_assays <- names(seu_export@assays)
assay_to_export <- if ("RNA" %in% available_assays) "RNA" else if ("sketch" %in% available_assays) "sketch" else available_assays[1]
cat(sprintf("  导出 assay: %s\n", assay_to_export))

## ---- 2. 提取 counts 矩阵（稀疏格式, genes x cells） ---- ##
cat("[2/7] 提取 counts 矩阵...\n")
counts <- GetAssayData(seu_export, assay = assay_to_export, layer = "counts")
cat(sprintf("  Counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))

## ---- 3. 提取 normalized data ---- ##
cat("[3/7] 提取 normalized data...\n")
data_norm <- tryCatch(
  GetAssayData(seu_export, assay = assay_to_export, layer = "data"),
  error = function(e) { cat("  警告: 无法提取 data layer\n"); NULL }
)

## ---- 4. 提取 obs（metadata） ---- ##
cat("[4/7] 提取 metadata (obs)...\n")
obs <- seu_export@meta.data

## 清理 obs 列名中的特殊字符（Python 兼容）
colnames(obs) <- gsub("[^A-Za-z0-9_.]", "_", colnames(obs))

## 清理每一列的数据类型，确保 h5py 能正确写入
for (col in colnames(obs)) {
  x <- obs[[col]]
  if (is.factor(x)) {
    ## factor -> character（避免 h5py 写入 mixed type）
    obs[[col]] <- as.character(x)
  } else if (is.logical(x)) {
    obs[[col]] <- as.character(x)
  } else if (is.numeric(x)) {
    ## numeric 列保持原样，NA 会被 anndata 正确处理
    next
  } else if (is.character(x)) {
    ## character 列中的 NA 替换为 "NA"（避免 h5py mixed type 报错）
    obs[[col]][is.na(obs[[col]])] <- "NA"
  } else {
    ## 其他类型强制转 character
    obs[[col]] <- as.character(x)
    obs[[col]][is.na(obs[[col]])] <- "NA"
  }
}

cat(sprintf("  Metadata: %d cells x %d 列\n", nrow(obs), ncol(obs)))
cat(sprintf("  包含列: %s\n", paste(colnames(obs), collapse = ", ")))

## ---- 5. 提取 var（基因信息 + HVG 标记） ---- ##
cat("[5/7] 构建 var (基因信息 + HVG)...\n")

var <- data.frame(
  row.names    = rownames(counts),
  gene_symbols = rownames(counts)
)

## 标记 HVG
hvg_genes <- tryCatch(VariableFeatures(seu_export), error = function(e) character(0))
var$highly_variable      <- rownames(var) %in% hvg_genes
var$highly_variable_rank <- NA_integer_
var$highly_variable_rank[match(hvg_genes, rownames(var))] <- seq_along(hvg_genes)
cat(sprintf("  HVG 数量: %d\n", sum(var$highly_variable)))

## 尝试把 Seurat var-level 统计信息也带上（mean, variance 等）
tryCatch({
  hvf_info <- HVFInfo(seu_export, assay = assay_to_export)
  for (col in colnames(hvf_info)) {
    safe_col <- gsub("[^A-Za-z0-9_.]", "_", col)
    var[[safe_col]] <- hvf_info[[col]][match(rownames(var), rownames(hvf_info))]
  }
  cat(sprintf("  附加 HVFInfo 列: %s\n", paste(colnames(hvf_info), collapse = ", ")))
}, error = function(e) {
  cat("  HVFInfo 无法提取（跳过）\n")
})

## ---- 6. 构建 AnnData 对象 ---- ##
cat("[6/7] 构建 AnnData 对象...\n")

## 将 counts 转为 Python scipy.sparse.csc_matrix（AnnData X 存储为 cells x genes）
counts_t <- Matrix::t(counts)  # cells x genes
adata <- sc$AnnData(
  X   = sp$csc_matrix(counts_t),  # X = raw counts (cells x genes, sparse)
  obs = obs,
  var = var
)
cat(sprintf("  AnnData shape: %d cells x %d genes\n", adata$shape[[1]], adata$shape[[2]]))

## 存储 counts 到 layers["counts"]（scVI 要求）
adata$layers["counts"] <- sp$csc_matrix(counts_t)

## 存储 normalized data 到 layers["data"]（可选, 供 scanpy 等直接使用）
if (!is.null(data_norm)) {
  data_norm_t <- Matrix::t(data_norm)
  adata$layers["data"] <- sp$csc_matrix(data_norm_t)
  cat("  已存储 layers: counts, data\n")
} else {
  cat("  已存储 layers: counts\n")
}

## ---- 7. 转移所有 reductions 到 obsm ---- ##
cat("[7/7] 转移 reductions 到 obsm...\n")
all_reductions <- Reductions(seu_export)
cat(sprintf("  共 %d 个 reductions: %s\n",
            length(all_reductions),
            paste(all_reductions, collapse = ", ")))

for (red in all_reductions) {
  tryCatch({
    emb <- Embeddings(seu_export, reduction = red)
    ## AnnData obsm key 格式: X_pca, X_harmony, ...
    ## 如果名字已经以 X_ 开头就不再加
    obsm_key <- if (startsWith(red, "X_")) red else paste0("X_", red)
    ## 用 numpy array 存储
    adata$obsm[obsm_key] <- np$array(emb)
    cat(sprintf("    [OK] %s -> obsm['%s'] (%d dims)\n", red, obsm_key, ncol(emb)))
  }, error = function(e) {
    cat(sprintf("    [FAIL] %s: %s\n", red, e$message))
  })
}

## ---- 额外: 存储 PCA loadings 到 varm（可选, 用于投影新数据） ---- ##
tryCatch({
  pca_loadings <- Loadings(seu_export, reduction = "pca")
  if (nrow(pca_loadings) > 0) {
    ## loadings: genes x n_pcs -> 存到 varm['PCs']
    ## 需要对齐到 var 的基因顺序
    common_genes <- intersect(rownames(pca_loadings), rownames(var))
    loadings_aligned <- matrix(0, nrow = nrow(var), ncol = ncol(pca_loadings))
    rownames(loadings_aligned) <- rownames(var)
    loadings_aligned[common_genes, ] <- pca_loadings[common_genes, ]
    adata$varm["PCs"] <- np$array(loadings_aligned)
    cat(sprintf("  PCA loadings -> varm['PCs'] (%d genes x %d PCs)\n",
                nrow(loadings_aligned), ncol(loadings_aligned)))
  }
}, error = function(e) {
  cat(sprintf("  PCA loadings 跳过: %s\n", e$message))
})

## ---- 额外: 记录 uns 信息（便于 scVI / benchmark 使用） ---- ##
adata$uns["batch_key"]    <- batch_var
adata$uns["n_hvg"]        <- as.integer(sum(var$highly_variable))
adata$uns["n_pcs"]        <- as.integer(n_pcs)
adata$uns["assay_source"] <- assay_to_export
adata$uns["reductions"]   <- paste(all_reductions, collapse = ",")
cat("  uns 信息已记录: batch_key, n_hvg, n_pcs, assay_source, reductions\n")

## ---- 保存 h5ad ---- ##
h5ad_path <- file.path(output_dir, "seurat_for_python.h5ad")
if (file.exists(h5ad_path)) file.remove(h5ad_path)

cat(sprintf("\n保存 h5ad: %s\n", h5ad_path))
adata$write_h5ad(h5ad_path)
cat(sprintf("  文件大小: %.2f GB\n", file.size(h5ad_path) / 1024^3))

## ---- 清理内存 ---- ##
rm(seu_export, counts, counts_t, data_norm, obs, var, hvg_genes)
if (exists("data_norm_t")) rm(data_norm_t)
if (exists("pca_loadings")) rm(pca_loadings)
if (exists("loadings_aligned")) rm(loadings_aligned)
gc()