# Hair Follicle Cell State Transition Predictor

這是一個基於 scRNA-seq 資料的 In Silico 預測系統。本系統可模擬特定基因在不同毛囊細胞類型中受到擾動（過度表現或抑制）後，細胞向 Anagen、Catagen 或 Telogen 時期轉移的機率。

## 功能特色
* **互動式介面**：使用 Streamlit 建立，支援下拉選單與自訂基因輸入。
* **極大化擾動演算法 (Max Impact)**：藉由強制放大 PCA 空間的權重影響力 (Loading Multiplier) 與柔化分類器邊界 (Temperature)，精確模擬基因對毛囊細胞路徑的推力。
* **視覺化輸出**：自動計算轉移機率，並產出擾動前後的長條圖對比分析。

## 資料庫準備 (重要)
啟動本系統需要預先處理好的單細胞資料庫 (`adata_full_genes.h5ad`)。
請至以下 Google Drive 連結下載資料庫檔案，並將其放置於本專案的根目錄下（與 `app.py` 放在同一層）：
🔗 **[點此下載模型資料庫 (Google Drive)](https://drive.google.com/drive/folders/1EOZUHP638H6nZSUHjFDhMcK6qPTVClej?usp=drive_link)**

## 安裝與執行環境

1. 複製此專案到本地端：
   ```bash
   git clone <你的 GitHub Repo 網址>
   cd <你的資料夾名稱>
