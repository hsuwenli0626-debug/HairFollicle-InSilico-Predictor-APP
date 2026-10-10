import streamlit as st
import os
import gdown
from hf_predictor import HairFolliclePredictor

# 頁面設定
st.set_page_config(page_title="毛囊細胞狀態轉移預測系統", layout="wide")
st.title("🧬 毛囊細胞狀態轉移預測系統 (In Silico Predictor)")
st.markdown("基於 scRNA-seq 資料，模擬特定基因擾動後，毛囊細胞向 Anagen、Catagen 或 Telogen 轉移的機率。")

DATA_PATH = "adata_full_genes.h5ad"

# ==========================================
# 你的 Google Drive 檔案 ID (已自動帶入)
# ==========================================
FILE_ID = "12Ry2S9HozStCvRqdd7ZlBYv3NQ9YLiWt" 

@st.cache_resource
def load_predictor():
    # 如果檔案不存在，則從 Google Drive 自動下載
    if not os.path.exists(DATA_PATH):
        st.info("首次啟動：正在從雲端下載模型資料庫，檔案較大，請耐心等候幾分鐘...")
        try:
            url = f'https://drive.google.com/uc?id={FILE_ID}'
            gdown.download(url, DATA_PATH, quiet=False)
            st.success("資料庫下載完成！")
        except Exception as e:
            st.error(f"下載失敗，請確認 Google Drive 連結是否已設定為「知道連結的人均可讀取」。錯誤訊息: {e}")
            st.stop()
            
    return HairFolliclePredictor(DATA_PATH)

hf_system = load_predictor()

# 依照你指定的細胞清單更新
cell_types = [
    "ORS_Epithelial", "SOX9_High_Progenitor", "Bulge_Progenitor",
    "Bulge_KRT15_High", "HG_Hair_Germ", "DP_Mesenchymal",
    "IRS", "Matrix_Candidate", "Endothelial_Noise", "Fibroblast_VIM_High"
]

# 側邊欄輸入介面
st.sidebar.header("設定基因擾動參數")
target_celltype = st.sidebar.selectbox("目標細胞種類 (Target Celltype)", cell_types)
current_stage = st.sidebar.selectbox("當前毛囊時期 (Current Stage)", ['Anagen', 'Catagen', 'Telogen'], index=1)
target_gene = st.sidebar.text_input("目標基因 (Target Gene)", value="PER1")
perturbation_type = st.sidebar.selectbox("擾動類型 (Perturbation Type)", ['overexpress', 'suppress'])
perturbation_value = st.sidebar.number_input("擾動倍率 (Perturbation Value)", value=100.0, step=10.0)
loading_multiplier = st.sidebar.number_input("PCA 影響力乘數 (Loading Multiplier)", value=10.0, step=1.0)
temperature = st.sidebar.number_input("決策邊界軟化溫度 (Temperature)", value=10.0, step=1.0)

# 執行預測
if st.sidebar.button("開始模擬預測", type="primary"):
    with st.spinner("模型運算中..."):
        results = hf_system.simulate(
            target_celltype=target_celltype, current_stage=current_stage, target_gene=target_gene,
            perturbation_type=perturbation_type, perturbation_value=perturbation_value,
            loading_multiplier=loading_multiplier, temperature=temperature
        )
        
        if results.get("error"):
            st.error(results["error"])
        else:
            st.subheader(f"預測結果: {target_gene.upper()} 基因擾動 ({perturbation_type})")
            probs = results["probabilities"]
            col1, col2, col3 = st.columns(3)
            col1.metric("Anagen 機率", f"{probs.get('Anagen', 0):.2%}")
            col2.metric("Catagen 機率", f"{probs.get('Catagen', 0):.2%}")
            col3.metric("Telogen 機率", f"{probs.get('Telogen', 0):.2%}")
            st.pyplot(results["figure"])
