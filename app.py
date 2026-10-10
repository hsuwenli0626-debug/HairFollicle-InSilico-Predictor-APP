import streamlit as st
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from hf_predictor import HairFolliclePredictor
import os
import gzip
import shutil
import gdown

st.set_page_config(page_title="Hair Follicle In Silico Predictor", layout="wide")
st.title("🧬 毛囊細胞狀態轉移預測系統 (In Silico Predictor)")
st.markdown("本系統基於 scRNA-seq 資料，可模擬特定基因擾動後，毛囊細胞向 Anagen、Catagen 或 Telogen 轉移的機率。")

GDRIVE_FOLDER_ID = "1EOZUHP638H6nZSUHjFDhMcK6qPTVClej"

@st.cache_resource
def load_system():
    file_path = "adata_full_genes_clean.h5ad"
    
    if not os.path.exists(file_path):
        with st.spinner("首次啟動：正在自 Google Drive 下載 9 個原始數據檔案... (約需 1~2 分鐘)"):
            gdown.download_folder(id=GDRIVE_FOLDER_ID, output="gdrive_data", quiet=False)
            
        with st.spinner("正在自動重組 10x 矩陣並計算特徵統計量..."):
            stages = ['Anagen', 'Catagen', 'Telogen']
            adatas = []
            
            os.makedirs("processed_10x", exist_ok=True)
            
            # 尋找下載下來的所有檔案
            downloaded_files = []
            for root, dirs, files in os.walk("gdrive_data"):
                for f in files:
                    downloaded_files.append(os.path.join(root, f))
                    
            for stage in stages:
                stage_folder = os.path.join("processed_10x", stage)
                os.makedirs(stage_folder, exist_ok=True)
                
                # 尋找對應時期的 barcodes, features/genes, matrix 檔案
                b_file, f_file, m_file = None, None, None
                for fp in downloaded_files:
                    fname = os.path.basename(fp).lower()
                    if stage.lower() in fname:
                        if 'barcodes' in fname:
                            b_file = fp
                        elif 'features' in fname or 'genes' in fname:
                            f_file = fp
                        elif 'matrix' in fname:
                            m_file = fp
                            
                if b_file and f_file and m_file:
                    # 複製並標準化命名為 10x 讀取看得懂的檔名
                    shutil.copy(b_file, os.path.join(stage_folder, "barcodes.tsv.gz"))
                    shutil.copy(f_file, os.path.join(stage_folder, "features.tsv.gz"))
                    shutil.copy(m_file, os.path.join(stage_folder, "matrix.mtx.gz"))
                    
                    try:
                        temp_adata = sc.read_10x_mtx(stage_folder, var_names='gene_symbols', cache=True)
                        temp_adata.obs['stage'] = stage
                        adatas.append(temp_adata)
                    except Exception as e:
                        print(f"讀取 {stage} 矩陣失敗: {e}")
                        
            if not adatas:
                raise FileNotFoundError("無法從雲端檔案中對應出 Anagen、Catagen、Telogen 的 10x 檔案結構。")
                
            adata = ad.concat(adatas, join='inner', batch_key='sample_batch')
            adata.var_names_make_unique()
            
            sc.pp.filter_cells(adata, min_genes=200)
            sc.pp.filter_genes(adata, min_cells=3)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=2000)
            
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, svd_solver='arpack')
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=0.8)
            
            cluster2celltype_updated = {
                '0': 'DP_Mesenchymal', '8': 'IRS', '10': 'Bulge_KRT15_High',
                '9': 'HG_Hair_Germ', '6': 'Bulge_Progenitor',
                '13': 'SOX9_High_Progenitor', '3': 'SOX9_High_Progenitor',
                '5': 'SOX9_High_Progenitor', '4': 'SOX9_High_Progenitor',
                '1': 'Matrix_Candidate',
                '18': 'Endothelial_Noise', '19': 'Fibroblast_VIM_High',
                '2': 'ORS_Epithelial', '7': 'ORS_Epithelial', '11': 'ORS_Epithelial',
                '12': 'ORS_Epithelial', '14': 'ORS_Epithelial', '15': 'ORS_Epithelial',
                '16': 'ORS_Epithelial', '17': 'ORS_Epithelial', '20': 'ORS_Epithelial',
                '21': 'ORS_Epithelial'
            }
            
            adata.obs['cell_type'] = adata.obs['leiden'].map(cluster2celltype_updated).astype('category')
            cells_to_keep = adata.obs['cell_type'] != 'Endothelial_Noise'
            adata_hf = adata[cells_to_keep].copy()
            
            # 計算基因均值與標準差供擾動模擬使用
            if sparse.issparse(adata_hf.X):
                adata_hf.var['mean'] = np.array(adata_hf.X.mean(axis=0)).flatten()
                var_sq = np.array(adata_hf.X.power(2).mean(axis=0)).flatten() - np.square(adata_hf.var['mean'])
                adata_hf.var['std'] = np.sqrt(np.maximum(var_sq, 1e-8))
            else:
                adata_hf.var['mean'] = np.mean(adata_hf.X, axis=0)
                adata_hf.var['std'] = np.std(adata_hf.X, axis=0)
                adata_hf.var['std'] = np.maximum(adata_hf.var['std'], 1e-8)

            sc.pp.scale(adata_hf, max_value=10)
            sc.tl.pca(adata_hf, svd_solver='arpack')
            
            adata_hf.write_h5ad(file_path)

    with st.spinner("正在初始化預測模型並訓練分類器..."):
        adata = sc.read_h5ad(file_path)
        return HairFolliclePredictor(adata)

try:
    hf_system = load_system()
    st.success("✅ 系統載入完成，請於左側設定模擬參數！")
except Exception as e:
    st.error(f"資料載入失敗: {e}")
    st.stop()

# --- 側邊欄參數設定 ---
st.sidebar.header("⚙️ 模擬參數設定")

target_celltype = st.sidebar.selectbox(
    "1. 目標細胞群 (Target Celltype)", 
    options=[
        'ORS_Epithelial',
        'SOX9_High_Progenitor',
        'Bulge_Progenitor',
        'Bulge_KRT15_High',
        'HG_Hair_Germ',
        'DP_Mesenchymal',
        'IRS',
        'Matrix_Candidate',
        'Fibroblast_VIM_High'
    ]
)

current_stage = st.sidebar.selectbox(
    "2. 當前毛囊時期 (Current Stage)", 
    options=['Anagen', 'Catagen', 'Telogen'], index=1
)

target_gene_input = st.sidebar.text_input("3. 擾動基因 (輸入一個或多個基因，以逗號分隔)", value="PER1")
target_gene = [g.strip().upper() for g in target_gene_input.split(',')] if ',' in target_gene_input else target_gene_input.strip().upper()

perturbation_type = st.sidebar.selectbox(
    "4. 擾動方式 (Perturbation Type)", 
    options=['overexpress', 'suppress']
)

perturbation_value = st.sidebar.slider("5. 擾動強度 (Perturbation Value)", 1, 100, 100)
loading_multiplier = st.sidebar.slider("6. PCA 權重放大 (Loading Multiplier)", 1.0, 50.0, 10.0, step=1.0)
temperature = st.sidebar.slider("7. 決策邊界軟化 (Temperature)", 1.0, 20.0, 10.0, step=0.5)

if st.button("🚀 執行預測模擬 (Run Simulation)", type="primary"):
    with st.spinner("正在進行高維度特徵投影與機率運算..."):
        fig = hf_system.simulate(
            target_celltype=target_celltype, 
            current_stage=current_stage, 
            target_gene=target_gene, 
            perturbation_type=perturbation_type, 
            perturbation_value=perturbation_value,
            loading_multiplier=loading_multiplier,
            temperature=temperature
        )
        
        if fig:
            st.pyplot(fig)
        else:
            st.warning("⚠️ 模擬失敗，請確認該細胞群是否在該時期存在，或基因名稱是否正確。")
