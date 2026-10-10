# Hair Follicle Cell State Transition Predictor

這是一個基於 scRNA-seq 資料的 In Silico 預測系統。本系統可模擬特定基因在不同毛囊細胞類型中受到擾動（過度表現或抑制）後，細胞向 Anagen、Catagen 或 Telogen 時期轉移的機率。

## 功能特色
* **互動式介面**：使用 Streamlit 建立，支援下拉選單與自訂基因輸入。
* **極大化擾動演算法 (Max Impact)**：藉由強制放大 PCA 空間的權重影響力 (Loading Multiplier) 與柔化分類器邊界 (Temperature)，精確模擬基因對毛囊細胞路徑的推力。
* **視覺化輸出**：自動計算轉移機率，並產出擾動前後的長條圖對比分析。
* **雲端自動下載**：首次啟動時，系統會自動從 Google Drive 下載預先訓練好的 `h5ad` 大型資料庫，無須手動配置。

## 資料庫準備
系統啟動時會自動抓取模型資料庫。若在特殊環境下自動下載失敗，您也可以手動點擊 [此 Google Drive 連結](https://drive.google.com/file/d/12Ry2S9HozStCvRqdd7ZlBYv3NQ9YLiWt/view?usp=sharing) 下載 `adata_full_genes.h5ad`，並將其放置於本專案的根目錄下（與 `app.py` 放在同一層）。

## 安裝與執行環境
   1. 複製此專案到本地端：
   ```bash
   git clone https://github.com/hsuwenli0626-debug/HairFollicle-InSilico-Predictor-APP.git
   cd HairFollicle-InSilico-Predictor-APP
   ```
   2.安裝必要的 Python 套件：
   ```bash
   pip install -r requirements.txt
   ```
   
   3.啟動 Streamlit 應用程式：
   ```bash
   streamlit run app.py
   ```

   
