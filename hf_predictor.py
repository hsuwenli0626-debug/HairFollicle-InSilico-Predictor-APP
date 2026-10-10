import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
import warnings

warnings.filterwarnings('ignore', category=FutureWarning)

class HairFolliclePredictor:
    def __init__(self, adata_path):
        print(f"正在載入模型資料庫: {adata_path} ...")
        self.adata = sc.read_h5ad(adata_path)
        
        if 'X_pca' not in self.adata.obsm:
            raise ValueError("找不到 PCA 降維資料。")
        if 'stage' not in self.adata.obs.columns:
            raise ValueError("找不到 'stage' 欄位。")
            
        self.X = self.adata.obsm['X_pca']
        self.y = self.adata.obs['stage'].values
        self.stages = np.unique(self.y)
        
        print("正在初始化系統與訓練狀態分類器...")
        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(self.X)
        
        self.classifier = LogisticRegression(solver='lbfgs', max_iter=2000, random_state=42)
        self.classifier.fit(X_scaled, self.y)
        
        y_pred = self.classifier.predict(X_scaled)
        acc = accuracy_score(self.y, y_pred)
        print(f"狀態分類器訓練完成！整體分類準確度: {acc*100:.2f}%")

    def simulate(self, target_celltype, current_stage, target_gene, perturbation_type, perturbation_value, 
                 loading_multiplier=10.0, temperature=10.0):
        target_gene = target_gene.upper()
        if target_gene not in self.adata.var_names:
            return {"error": f"在資料庫中找不到基因 '{target_gene}'。"}
            
        mask = (self.adata.obs['cell_type'] == target_celltype) & (self.adata.obs['stage'] == current_stage)
        target_indices = np.where(mask)[0]
        
        if len(target_indices) == 0:
            return {"error": f"找不到符合條件的細胞 ({current_stage} 的 {target_celltype})。"}
            
        target_X = self.adata.X[target_indices, :].copy()
        if type(target_X) != np.ndarray:
            target_X = target_X.toarray()
            
        modified_loadings = self.adata.varm['PCs'].copy()
        gene_idx = self.adata.var_names.get_loc(target_gene)
        
        modified_loadings[gene_idx, :] = modified_loadings[gene_idx, :] * loading_multiplier
        
        g_mean = self.adata.var['mean'].values[gene_idx]
        g_std = self.adata.var['std'].values[gene_idx]
        
        log1p_vals = (target_X[:, gene_idx] * g_std) + g_mean
        log1p_vals = np.maximum(0, log1p_vals) 
        original_exp = np.expm1(log1p_vals) 
        
        if perturbation_type == 'suppress':
            alpha = max(0.0, 1.0 - (perturbation_value / 10.0))
            new_exp = original_exp * alpha
        else: 
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
        
        ax.set_ylabel('Probability (%)', fontsize=12)
        ax.set_title(f'State Transition Probability\nCell: {target_celltype} | Stage: {current_stage} | Gene: {target_gene} ({perturbation_type})', fontsize=11)
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
        
        result_probs = {cls: prob for cls, prob in zip(classes, avg_probs)}
        
        return {
            "error": None,
            "figure": fig,
            "probabilities": result_probs
        }
