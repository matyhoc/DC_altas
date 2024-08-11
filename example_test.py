##################examples: predict labels
from dc_predict import *

adata_DC=sc.read("adata_test.h5ad") #make sure the data has been processed by log1p-normalized with size-factor 1e4
adata_DC=pred_dc_ensemble(adata_DC,n_models=10,n_jobs=40,copy=True,seed=0)#anndata: predicted labels saved in .obs['predicted_label']
#df_pred=pred_dc_ensemble(adata_DC,n_models=10,n_jobs=10,copy=False,seed=0)#series:predicted labels saved as pd.Series


#################examples: extract highly expressed genes and lowly expressed genes
marker_list=[]
all_celltypes=['AS DC','DC_DPYD','DC_LAMP3','DC_MKI67','DC_pro_CD1C','DC_pro_PRDM16','LC','cDC1','cDC2','pDC']
high_markers,low_markers=extract_marker_genes(n_models=10,copy=True,seed=0,top_genes_each_model=20)

#format high_markers and low_markers into pd.DataFrame

df_high_markers = pd.DataFrame(dict([(key, pd.Series(value)) for key, value in high_markers.items()]))
df_high_markers=df_high_markers.replace(np.nan, '')
df_low_markers=pd.DataFrame(dict([(key, pd.Series(value)) for key, value in low_markers.items()]))
df_low_markers=df_low_markers.replace(np.nan, '')

