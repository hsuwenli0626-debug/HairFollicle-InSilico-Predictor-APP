import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
from hf_predictor import HairFolliclePredictor
import os
import gdown

st.set_page_config(page_title="Hair Follicle In Silico Predictor", layout="wide")
st.title("🧬 毛囊細胞狀態轉移預測系統 (In Silico Predictor)")
st.markdown("本系統基於 scRNA-seq 資料，可模擬特定基因擾動後，毛囊細胞向 Anagen、Catagen 或 Telogen 轉移的機率。")

# Google Drive 設定：
# 預設使用先前建立的 Google Drive 資料夾/檔案 ID
GDRIVE_FOLDER_ID = "1B0QFE8Z231glvLGqZ5ld3x9NnB9IABPP"

@st.cache_resource
def load_system():
    file_path = "adata_full_genes_clean.h5ad"
    
    # 若本地尚未存在資料庫，從 Google Drive 下載
    if not os.path.exists(file_path):
        with st.spinner("首次啟動：正在自 Google Drive 下載單細胞數據庫... (約需 1~3 分鐘)"):
            try:
                # 優先嘗試以單檔形式下載
                url = f"https://drive.google.com/uc?id={GDRIVE_FOLDER_ID}"
                gdown.download(url, file_path, quiet=False)
            except Exception:
                # 若為資料夾則下載整個資料夾
                gdown.download_folder(id=GDRIVE_FOLDER_ID, output="gdrive_data", quiet=False)
                # 搜尋資料夾內的 .h5ad 檔案
                for root, dirs, files in os.walk("gdrive_data"):
                    for f in files:
                        if f.endswith(".h5ad"):
                            file_path = os.path.join(root, f)
                            break

    if not os.path.exists(file_path):
        raise FileNotFoundError("未能找到 adata_full_genes_clean.h5ad，請確認 Google Drive 連結與權限設定。")

    with st.spinner("正在初始化預測模型並訓練分類器..."):
        adata = sc.read_h5ad(file_path)
        return HairFolliclePredictor(adata)

try:
    hf_system = load_system()
    st.success("✅ 系統載入完成，請於左側設定模擬參數！")
except Exception as e:
    st.error(f"資料載入失敗: {e}")
    st.info("提示：若上傳至 Google Drive 的是單一 .h5ad 檔案，請將 app.py 中的 GDRIVE_FOLDER_ID 替換為該檔案之專屬 File ID。")
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
