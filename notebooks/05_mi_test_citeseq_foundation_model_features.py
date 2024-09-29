#!/usr/bin/env python
# coding: utf-8

# The aim of this notebook is to do an assessment of mutual information between the CITE-seq ADT (surface protein expression) and cytoplasmic RNA expression (SCT) in the PBMC cells. The MI will be calculated once again by the LatentMI method, but in this case, we'll test the method with different features - including adding scGPT, Orthrus, and ESM2 embeddings to the cell and protein spaces. 
# 
# scGPT features will be added to the cell, and protein spaces, while Orthrus and ESM2 embeddings will be added to the protein space. The embeddings will be added to the protein space by concatenating the embeddings to the protein expression matrix. Similarly, the scGPT embeddings will be added to the cell space by concatenating the embeddings to the SCT matrix.
# 
# Using the ceiling of the number of features for 2000 HVGs + (scGPT + Orthrus + ESM2 embeddings), we'll also include a baseline with an increased number of HVGs. E.g. if the Orthrus features have the highest number of dimensions, we'll use the top 2000 HVGs + n HVG features (n total orthrus dimensions).
# 
# Each iteration will be ran 10 times on different random samples of cells in the PBMC cite-seq dataset.

# In[1]:


import pandas as pd 
import numpy as np 
import scanpy as sc 
import anndata as ann 
import muon as mu
import mudata as md 
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


##Load the MuData object with scgpt, orthrus and esm2 embeddings
mdata = md.read("/fs01/projects/isoclr/cite_seq_with_seq_embed_with_cell_embed.h5mu")


# In[3]:


# Extract all of the scgpt embeddings 
scgpt_cols = [col for col in mdata["SCT"].obs.columns if "scgpt" in col]
scgpt_col_indices = [mdata["SCT"].obs.columns.get_loc(col) for col in scgpt_cols]
scgpt_embeddings = np.array(mdata["SCT"].obs[scgpt_cols].values)

# Extract all of the orthrus embeddings
orthrus_cols = [col for col in mdata["ADT"].var.columns if "orthrus" in col]
orthrus_col_indices = [mdata["ADT"].var.columns.get_loc(col) for col in orthrus_cols]
orthrus_embeddings = np.array(mdata["ADT"].var[orthrus_cols].values)

# Extract all of the esm2 embeddings
esm2_cols = [col for col in mdata["ADT"].var.columns if "esm" in col]
esm2_col_indices = [mdata["ADT"].var.columns.get_loc(col) for col in esm2_cols]
esm2_embeddings = np.array(mdata["ADT"].var[esm2_cols].values)


# In[4]:


# Within mdata, drop any vars that have nan values for either the esm2 or orthrus scolumns 
orthrus_cols_subset = mdata["ADT"].var[orthrus_cols]
esm2_cols_subset = mdata["ADT"].var[esm2_cols]

nan_mask_orthrus = orthrus_cols_subset.isna().any(axis=1)
nan_mask_esm2 = esm2_cols_subset.isna().any(axis=1)

mdata_adt_subset = mdata["ADT"][:, ~nan_mask_orthrus & ~nan_mask_esm2]


# In[5]:


# Using the same mask, subset the esm2 and orthrus embeddings
orthrus_embeddings_subset = orthrus_embeddings[~nan_mask_orthrus & ~nan_mask_esm2]
esm2_embeddings_subset = esm2_embeddings[~nan_mask_orthrus & ~nan_mask_esm2]


# In[6]:


# Determine embedding dimensions
print(f"orthrus embeddings shape: {orthrus_embeddings_subset.shape}")
print(f"esm2 embeddings shape: {esm2_embeddings_subset.shape}")


# In[6]:


# Start with the HVG baseline 
mdata_sct_sub = mdata["SCT"]

# Get HVG subset of SCT
sc.pp.highly_variable_genes(mdata_sct_sub, n_top_genes=2000)
sct_counts = mdata["SCT"].X[:, mdata_sct_sub.var.highly_variable]
adt_counts = mdata_adt_subset.X

# Randomly sample 1000 cells with a seed - repeat 10 times 
from latentmi import lmi

n_mi_calc = 10
#MI_real = np.zeros(n_mi_calc)
#for i in range(n_mi_calc):
#    np.random.seed(i)
#    indices = np.random.choice(mdata["SCT"].shape[0], 1000, replace=True)
#    sct_counts_sub = sct_counts[indices, :].todense()
#    adt_counts_sub = adt_counts[indices, :]
#    pmis, _, _ = lmi.estimate(sct_counts_sub, adt_counts_sub)
#    MI_real[i] = np.nanmean(pmis)
    
# Save the MI estimate from the calculation and the MI estimates from the simulations
#np.save("../data/MI_estimate_1k_10_times_adt_sct_exp_5.npy", MI_real)


# In[7]:


# Now test ESM2 - get 2000 HVGs and test by appending the esm2 embeddings to the ADT data

# Randomly sample 1000 cells with a seed - repeat 10 times

MI_real_esm2 = np.zeros(n_mi_calc)
for i in range(n_mi_calc):
    np.random.seed(i)
    indices = np.random.choice(mdata["SCT"].shape[0], 1000, replace=True)
    sct_counts_sub = sct_counts[indices, :].todense()
    adt_counts_sub = adt_counts[indices, :]
    # Dot product of adt_counts_sub and esm2_embeddings
    esm2_embeddings_subset = np.abs(esm2_embeddings_subset)
    adt_esm_product = np.dot(adt_counts_sub, esm2_embeddings_subset)
    adt_esm_product_adt = np.concatenate((adt_counts_sub, adt_esm_product), axis=1)
    pmis, _, _ = lmi.estimate(sct_counts_sub, adt_esm_product_adt)
    MI_real_esm2[i] = np.nanmean(pmis)
    
# Save the MI estimate from the calculation and the MI estimates from the simulations
np.save("../data/MI_estimate_1k_10_times_adt_sct_exp_5_esm2.npy", MI_real_esm2)


# In[8]:


# Now test Orthrus - get 2000 HVGs and test by appending the orthrus embeddings to the ADT data

# Randomly sample 1000 cells with a seed - repeat 10 times
MI_real_orthrus = np.zeros(n_mi_calc)
for i in range(n_mi_calc):
    np.random.seed(i)
    indices = np.random.choice(mdata["SCT"].shape[0], 1000, replace=True)
    sct_counts_sub = sct_counts[indices, :].todense()
    adt_counts_sub = adt_counts[indices, :]
    # Dot product of adt_counts_sub and orthrus_embeddings
    orthrus_embeddings_subset = np.abs(orthrus_embeddings_subset)
    adt_orthrus_product = np.dot(adt_counts_sub, orthrus_embeddings_subset)
    adt_orthrus_product_adt = np.concatenate((adt_counts_sub, adt_orthrus_product), axis=1)
    pmis, _, _ = lmi.estimate(sct_counts_sub, adt_orthrus_product_adt)
    MI_real_orthrus[i] = np.nanmean(pmis)
    
# Save the MI estimate from the calculation and the MI estimates from the simulations
np.save("../data/MI_estimate_1k_10_times_adt_sct_exp_5_orthrus.npy", MI_real_orthrus)


# In[9]:


# Now test SCGPT - substitute the hvg counts for the scgpt embeddings
#sct_counts_scgpt = scgpt_embeddings
#adt_counts = mdata_adt_subset.X

# Randomly sample 1000 cells with a seed - repeat 10 times
# MI_real_scgpt = np.zeros(n_mi_calc)
#for i in range(n_mi_calc):
#    np.random.seed(i)
#    indices = np.random.choice(mdata["SCT"].shape[0], 1000, replace=True)
#    sct_counts_sub = sct_counts_scgpt[indices, :]
#    adt_counts_sub = adt_counts[indices, :]
#    pmis, _, _ = lmi.estimate(sct_counts_sub, adt_counts_sub)
#    MI_real_scgpt[i] = np.nanmean(pmis)

# Save the MI estimate from the calculation and the MI estimates from the simulations
#np.save("../data/MI_estimate_1k_10_times_adt_sct_exp_5_scgpt.npy", MI_real_scgpt)


# In[ ]:


# Now test including the scgpt and hvg counts together
sct_counts_scgpt = np.concatenate((scgpt_embeddings, sct_counts), axis=1) 
adt_counts = mdata_adt_subset.X

# Randomly sample 1000 cells with a seed - repeat 10 times
MI_real_scgpt_hvg = np.zeros(n_mi_calc)
for i in range(n_mi_calc):
    np.random.seed(i)
    indices = np.random.choice(mdata["SCT"].shape[0], 1000, replace=True)
    sct_counts_sub = sct_counts_scgpt[indices, :]
    adt_counts_sub = adt_counts[indices, :]
    pmis, _, _ = lmi.estimate(sct_counts_sub, adt_counts_sub)
    MI_real_scgpt_hvg[i] = np.nanmean(pmis)
    
# Save the MI estimate from the calculation and the MI estimates from the simulations
np.save("../data/MI_estimate_1k_10_times_adt_sct_exp_5_scgpt_hvg.npy", MI_real_scgpt_hvg)

