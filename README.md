# Hair Follicle In Silico Predictor

本專案是一個基於單細胞 RNA 定序（scRNA-seq）技術開發的 *in silico* 毛囊細胞狀態預測與模擬系統。

系統的底層訓練資料源自於公開數據庫 **GEO: GSE205456**。透過訓練多項式羅吉斯迴歸分類器（Multinomial Logistic Regression），本系統能夠以量化的方式模擬並預測：當特定的毛囊細胞群體遭受基因的虛擬擾動時，其細胞命運發生跨週期狀態轉移（State Transition）的機率變化。

## 🌐 雲端互動式網頁 (Web App)
本專案支援透過 **Streamlit Community Cloud** 進行雲端部署，並透過 `gdown` 自動掛載 Google Drive 之數據矩陣。
- **主程式入口**：`app.py`
- **佈署平台**：[Streamlit Community Cloud](https://share.streamlit.io/)

## 資料庫配置說明
受限於 GitHub 單一檔案 100MB 限制，矩陣資料庫託管於 Google Drive：
- **Google Drive 資料夾連結**：[點此前往](https://drive.google.com/drive/folders/1EOZUHP638H6nZSUHjFDhMcK6qPTVClej?usp=drive_link)
