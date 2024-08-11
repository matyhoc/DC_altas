import celltypist
from celltypist import models #pip install celltypist
import scanpy as sc
import numpy as np
import pandas as pd
import os
from functools import reduce
from joblib import Parallel, delayed
model_abspath=os.path.dirname(os.path.abspath(__file__))
models.models_path=model_abspath+"/models/"

def pred_dc_single(adata,model_cur):
    predictions = celltypist.annotate(adata, model =model_cur,majority_voting=False)
    return predictions.predicted_labels.copy()

#load models saved in "./models"
models_list=[models.models_path+model0 for model0 in os.listdir(models.models_path) if str(model0).startswith("dc_models")]

def pred_dc_ensemble(adata,n_models=10,n_jobs=1,copy=True,seed=0):
    """
    adata_DC: processed anndata format; Make sure adata.X. has been log1p-normalized, that is: adata_DC.X.expm1().sum(1) return all 1e4
    n_modesl: can be set integer 1-40;
    n_jobs: The maximum number of concurrently running jobs,If -1 all CPUs are used.
 |      If 1 is given, no parallel computing code is used at all, and the
 |      behavior amounts to a simple python `for` loop. This mode is not
 |      compatible with `timeout`.
 |      For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for
 |      n_jobs = -2, all CPUs but one are used.
    copy: if true, return a anndata object that .obs['predicted_labels']; if false return a pd.DataFrame with a column names "predicted_labels"
    seed=0: if n_models<40, we will randomly choose n_models with seed
    """
    np.random.seed(seed)
    if n_models <1 or n_models>40:
        n_models=40 #
        print("Using all models to predict")

    df_pred_list=Parallel(n_jobs=n_jobs)(delayed(pred_dc_single)(adata,model0) for model0 in [models_list[i] for i in np.random.choice(40,n_models)])
    for id0 in range(len(df_pred_list)):
        df_pred_tmp=df_pred_list[id0]
        df_pred_tmp.columns=['predicted_labels_'+'{0:02}'.format(id0+1)]
        df_pred_list[id0]=df_pred_tmp.copy()

    def row_max_index(x):
        tmp0=pd.Series(x).value_counts(sort=True)
        return tmp0.index.values[0]

    df_pred=reduce(lambda x,y: pd.merge(x,y, left_index=True, right_index=True,how='left'), df_pred_list)
    df_final=df_pred.apply(lambda x: row_max_index(x) ,axis=1)

    if not copy:
        return df_final 
    else:
        adata.obs['predicted_label']=df_final.values
        return adata


def extract_marker_genes(n_models=10,copy=True,seed=0,top_genes_each_model=20):
    np.random.seed(seed)
    if n_models <1 or n_models>40:
        n_models=40 #
        print("Using all models to predict")
    high_degenes={}
    low_degenes={}
    all_celltypes=['AS DC','DC_DPYD','DC_LAMP3','DC_MKI67','DC_pro_CD1C', 'DC_pro_PRDM16','LC','cDC1','cDC2','pDC']
    for id0,model0 in enumerate([models_list[i] for i in np.random.choice(40,n_models)]):
        model=models.Model.load(model0)
        for i in all_celltypes:
            if id0==0:
                high_degenes[i]=[]
                low_degenes[i]=[]
            high_degenes[i]=high_degenes[i]+list(model.extract_top_markers(i,top_genes_each_model,only_positive=True))
            low_degenes[i]=low_degenes[i]+list(set(model.extract_top_markers(i,top_genes_each_model*2,only_positive=False)).difference(set(high_degenes[i])))

    high_markers={}
    low_markers={}
    for i in all_celltypes:
        tmp0=pd.Series(high_degenes[i])
        high_markers[i]=list(tmp0.value_counts().index[tmp0.value_counts()>=(tmp0.shape[0]/top_genes_each_model*0.4)])
        tmp0=pd.Series(low_degenes[i])
        low_markers[i]=list(tmp0.value_counts().index[tmp0.value_counts()>=(tmp0.shape[0]/top_genes_each_model*0.4)])

    return high_markers,low_markers

if __name__=="__main__":
    ##################examples: predict labels
    adata_DC=sc.read("adata_test.h5ad")
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
