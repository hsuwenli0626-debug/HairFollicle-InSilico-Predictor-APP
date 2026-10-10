import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)

class HairFolliclePredictor:
    def __init__(self, adata):
        print("正在初始化系統與訓練狀態分類器...")
        self.adata = adata.copy()
        
        if 'X_pca' not in self.adata.obsm:
            raise ValueError("找不到 PCA 降維資料。請先確保資料包含 .obsm['X_pca']")
        if 'stage' not in self.adata.obs.columns:
            raise ValueError("找不到 'stage' 欄位。")
            
        # 定義給 GitHub 使用者的前台簡化標籤映射
        self.simplified_map = {
            'DP': ['DP_Mesenchymal'],
            'Bulge_Stem_Cells': ['Bulge_KRT15_High', 'Bulge_Progenitor'],
            'Matrix_HG': ['Matrix_Candidate', 'HG_Hair_Germ'],
            'ORS': [c for c in self.adata.obs['cell_type'].unique() if 'ORS' in str(c)]
        }
            
        self.X = self.adata.obsm['X_pca']
        self.y = self.adata.obs['stage'].values
        self.stages = np.unique(self.y)
        
        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(self.X)
        
        X_train, X_test, y_train, y_test = train_test_split(X_scaled, self.y, test_size=0.2, random_state=42)
        
        self.classifier = LogisticRegression(multi_class='multinomial', solver='lbfgs', max_iter=2000, random_state=42)
        self.classifier.fit(X_train, y_train)
        
        y_pred = self.classifier.predict(X_test)
        acc = accuracy_score(y_test, y_pred)
        print(f"✅ 狀態分類器訓練完成！Anagen, Catagen, Telogen 分類準確度: {acc*100:.2f}%\n" + "-"*50)

    def simulate(self, target_celltype, current_stage, target_gene, perturbation_type, perturbation_value, 
                 loading_multiplier=3.0, temperature=1.5):
                 
        if isinstance(target_gene, str):
            target_genes = [target_gene.upper()]
        elif isinstance(target_gene, list):
            target_genes = [g.upper() for g in target_gene]
        else:
            raise ValueError("target_gene 必須是字串或字串列表。")
            
        valid_genes = [g for g in target_genes if g in self.adata.var_names]
        if not valid_genes:
            print(f"❌ 錯誤：在資料庫中找不到指定的基因。")
            return None
            
        genes_str = ", ".join(valid_genes)
        
        if target_celltype in self.simplified_map:
            target_subtypes = self.simplified_map[target_celltype]
            mask = (self.adata.obs['cell_type'].isin(target_subtypes)) & (self.adata.obs['stage'] == current_stage)
            display_celltype = f"{target_celltype} (Merged Group)"
        else:
            mask = (self.adata.obs['cell_type'] == target_celltype) & (self.adata.obs['stage'] == current_stage)
            display_celltype = target_celltype
            
        target_indices = np.where(mask)[0]
        
        if len(target_indices) == 0:
            print(f"❌ 錯誤：找不到符合條件的細胞 ({current_stage} 的 {display_celltype})。")
            return None
            
        print(f"🚀 啟動模擬：針對 {len(target_indices)} 顆細胞")
        print(f"🧬 基因: [{genes_str}] | 擾動: {perturbation_type} ({perturbation_value}x) | 權重槓桿: {loading_multiplier}x")
        
        target_X = self.adata.X[target_indices, :].copy()
        if type(target_X) != np.ndarray:
            target_X = target_X.toarray()
            
        modified_loadings = self.adata.varm['PCs'].copy()
            
        for g in valid_genes:
            gene_idx = self.adata.var_names.get_loc(g)
            
            modified_loadings[gene_idx, :] = modified_loadings[gene_idx, :] * loading_multiplier
            
            g_mean = self.adata.var['mean'].values[gene_idx]
            g_std = self.adata.var['std'].values[gene_idx]
            
            log1p_vals = (target_X[:, gene_idx] * g_std) + g_mean
            log1p_vals = np.maximum(0, log1p_vals) 
            original_exp = np.expm1(log1p_vals) 
            
            if perturbation_type == 'suppress':
                alpha = max(0.0, 1.0 - (perturbation_value / 10.0))
                new_exp = original_exp * alpha
            elif perturbation_type == 'overexpress':
                new_exp = original_exp * perturbation_value
                
            new_log1p_vals = np.log1p(new_exp)
            new_scaled_vals = (new_log1p_vals - g_mean) / g_std
            new_scaled_vals = np.clip(new_scaled_vals, a_min=None, a_max=100.0)
            
            target_X[:, gene_idx] = new_scaled_vals
            
        perturbed_pca = np.dot(target_X, modified_loadings)
        perturbed_pca_scaled = self.scaler.transform(perturbed_pca)
        
        original_pca = self.X[target_indices, :]
        original_pca_scaled = self.scaler.transform(original_pca)
        
        decision_perturbed = self.classifier.decision_function(perturbed_pca_scaled)
        decision_original = self.classifier.decision_function(original_pca_scaled)
        
        def temperature_softmax(logits, T):
            logits = logits / T
            exp_logits = np.exp(logits - np.max(logits, axis=1, keepdims=True))
            return exp_logits / np.sum(exp_logits, axis=1, keepdims=True)
        
        if decision_perturbed.ndim == 1:
             decision_perturbed = np.c_[-decision_perturbed, decision_perturbed]
             decision_original = np.c_[-decision_original, decision_original]
             
        probs = temperature_softmax(decision_perturbed, temperature)
        orig_probs = temperature_softmax(decision_original, temperature)
        
        avg_probs = np.mean(probs, axis=0)
        avg_orig_probs = np.mean(orig_probs, axis=0)
        
        classes = self.classifier.classes_
        fig, ax = plt.subplots(figsize=(8, 5))
        x = np.arange(len(classes))
        width = 0.35
        
        rects1 = ax.bar(x - width/2, avg_orig_probs * 100, width, label='Before Perturbation', color='#D0D0D0')
        rects2 = ax.bar(x + width/2, avg_probs * 100, width, label='After Perturbation', color='#E64B35')
        
        display_genes = genes_str if len(genes_str) < 30 else f"{len(valid_genes)} Target Genes"
        ax.set_ylabel('Probability (%)', fontsize=12)
        ax.set_title(f'Pathway Transition Probability\nCell: {display_celltype} | Stage: {current_stage}\nPathway: [{display_genes}] ({perturbation_type}: {perturbation_value})', fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels(classes, fontsize=12)
        ax.legend(fontsize=11)
        
        ax.bar_label(rects1, fmt='%.1f%%', padding=3, fontsize=10, color='#505050')
        ax.bar_label(rects2, fmt='%.1f%%', padding=3, fontsize=10, color='black')
        
        for i in range(len(classes)):
            diff = (avg_probs[i] - avg_orig_probs[i]) * 100
            if abs(diff) >= 0.1:
                sign = "+" if diff > 0 else ""
                color = "red" if diff > 0 else "blue"
                ax.text(x[i] + width/2, (avg_probs[i]*100) + 7, f'{sign}{diff:.1f}%', ha='center', va='bottom', color=color, fontweight='bold')

        max_y = max(max(avg_orig_probs), max(avg_probs)) * 100
        ax.set_ylim(0, min(120, max_y + 25))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout()
        
        return fig
