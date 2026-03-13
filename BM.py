#### ========== Part 4: Python - 加载 h5ad & 准备 Benchmark ========== ####
## 【Python 内核 (scib-env)】

import os

# 避免 OpenBLAS/OMP 在线程嵌套并行时触发线程上限和内存区域分配错误
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors
from scib_metrics.nearest_neighbors import NeighborsResults
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

## 参数配置
input_dir    = "/home/data/tanglei/project/prostate_altas/output/02"
h5ad_path    = os.path.join(input_dir, "benchmark_GSM.h5ad")
batch_key    = "sample.ID"
celltype_key = "celltype"

## 加载 h5ad
print(f"加载 h5ad: {h5ad_path}")
adata = sc.read_h5ad(h5ad_path)
print(f"数据维度: {adata.shape[0]} cells x {adata.shape[1]} genes")
print(f"批次数 ({batch_key}): {adata.obs[batch_key].nunique()}")
print(f"细胞类型数 ({celltype_key}): {adata.obs[celltype_key].nunique()}")
print(f"\n可用 obsm keys: {sorted(adata.obsm.keys())}")
print(f"可用 layers: {list(adata.layers.keys())}")

## 确认必要的 embedding 存在
required_keys = ["X_pca", "X_cca", "X_rpca", "X_harmony", "X_mnn",
                 "X_liger", "X_combat", "X_scvi", "X_scanvi",
                 "X_scanorama", "X_bbknn","X_stacas"]

print("\n=== Embedding 检查 ===")
available_keys = []
for key in required_keys:
    if key in adata.obsm:
        shape = adata.obsm[key].shape
        print(f"  [OK]   {key:<25} shape: {shape}")
        available_keys.append(key)
    else:
        print(f"  [缺失] {key}")

print(f"\n共 {len(available_keys)} 个 embedding 可用于 Benchmark")

#### ========== Part 5: Python - 运行 scib_metrics Benchmark ========== ####
## 【Python 内核 (scib-env)】
##
## 指标说明:
##   Bio Conservation: Isolated labels, KMeans NMI/ARI, Silhouette label, cLISI
##   Batch Correction:  Silhouette batch, iLISI, kBET, Graph connectivity, PCR comparison

from pathlib import Path
import pickle

def sklearn_nn(X: np.ndarray, k: int) -> NeighborsResults:
    """CPU brute-force KNN（避免 GPU/pynndescent 依赖问题）"""
    X = np.ascontiguousarray(X, dtype=np.float32)
    nbrs = NearestNeighbors(
        n_neighbors=k,
        algorithm="brute",
        metric="euclidean",
        n_jobs=1
    ).fit(X)
    distances, indices = nbrs.kneighbors(X)
    return NeighborsResults(indices=indices, distances=distances)


## 整合方法列表（排除未校正的 pca 和 umap 类）
integration_keys = [k for k in available_keys if k != "X_pca"]
print(f"Benchmark 评测方法 ({len(integration_keys)} 个):\n  " +
      ", ".join(integration_keys))
print(f"未校正基线: X_pca")
print(f"Batch key:  {batch_key}")
print(f"Label key:  {celltype_key}")

## Benchmarker 缓存文件（保存到 input_dir 下，可在 tmux 先跑再回 notebook 复用）
bm_cache_path = Path(input_dir) / "bm_scib_metrics_GSM.pkl"
bm_cache_path.parent.mkdir(parents=True, exist_ok=True)

if bm_cache_path.exists():
    print(f"\n检测到缓存，直接加载 bm: {bm_cache_path.resolve()}")
    with bm_cache_path.open("rb") as f:
        bm = pickle.load(f)
    print("bm 加载完成，可直接运行下面结果展示单元格。")
else:
    print("\n未检测到 bm 缓存，开始完整 benchmark 计算...")

    biocons = BioConservation(isolated_labels=True)
    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=celltype_key,
        embedding_obsm_keys=integration_keys,
        pre_integrated_embedding_obsm_key="X_pca",
        bio_conservation_metrics=biocons,
        batch_correction_metrics=BatchCorrection(),
        n_jobs=4,
    )

    print("\n开始 prepare（计算 KNN）...")
    bm.prepare(neighbor_computer=sklearn_nn)

    print("\n开始 benchmark（计算各项指标）...")
    bm.benchmark()

    with bm_cache_path.open("wb") as f:
        pickle.dump(bm, f)

    print(f"\nBenchmark 完成，bm 已保存到: {bm_cache_path.resolve()}")