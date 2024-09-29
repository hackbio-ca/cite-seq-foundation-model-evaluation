#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd 
import numpy as np 
import scanpy as sc 
import anndata as ann 
import muon as mu
import mudata as md 
import matplotlib.pyplot as plt
import seaborn as sns


# The aim of this notebook is to do an assessment of mutual information between the CITE-seq ADT (surface protein expression) and cytoplasmic RNA expression (SCT) in the PBMC cells. 
# 
# The MI will be estimated using LMI (https://arxiv.org/abs/2409.02732). First install the LMI library:

# In[2]:

# In[3]:


##Load the MuData object (This is a big file - should have at last 16GB of RAM to prevent memory issues)
mdata = md.read("../data/multi.h5mu")


# In[4]:


# Extract and subset the counts for SCT and ADT. We'll sample 1000 cells to start with to speed up the analysis


# In[5]:


#indices = np.random.choice(mdata["SCT"].shape[0], 1000, replace=False)
#sct_counts = mdata["SCT"].X[indices, :].todense()
#adt_counts = mdata["ADT"].X[indices, :]


# In[6]:


# Test the latentmi function in this respect 
# from latentmi import lmi

#pmis, embedding, model = lmi.estimate(sct_counts, adt_counts)


# 10000 cells take too long of a time- test with 1000 instead - this is also taking a long time. 
# 
# Try with 2500 HVGs instead and 1000 cells

# In[9]:


# Get HVG subset of SCT
sc.pp.highly_variable_genes(mdata["SCT"], n_top_genes=2000)
sct_counts = mdata["SCT"].X[:, mdata["SCT"].var.highly_variable]
adt_counts = mdata["ADT"].X

# Randomly sample 1000 cells with a seed - repeat 10 times 
from latentmi import lmi

n_mi_calc = 10
# MI_real = np.zeros(n_mi_calc)
# for i in range(n_mi_calc):
#     np.random.seed(i)
#     indices = np.random.choice(mdata["SCT"].shape[0], 1000, replace=True)
#     sct_counts_sub = sct_counts[indices, :].todense()
#     adt_counts_sub = adt_counts[indices, :]
#     pmis, _, _ = lmi.estimate(sct_counts_sub, adt_counts_sub)
#     MI_real[i] = np.nanmean(pmis)


# # In[ ]:


# # Now simulate a bunch of random gaussian distributions with the same dimensions as the adt and sct counts and calculate the MI for each of them
# # This will give us a null distribution to compare the MI estimate to

# n_simulations = 20
# MI_simulations = np.zeros(n_simulations)
# for i in range(n_simulations):
#     seed = np.random.randint(0, 2**32 - 1)
#     np.random.seed(seed)
#     sct_sim = np.random.normal(size=sct_counts.shape)
#     adt_sim = np.random.normal(size=adt_counts.shape)
#     pmis, _, _ = lmi.estimate(sct_sim, adt_sim)
#     MI_simulations[i] = np.nanmean(pmis)


# # In[ ]:


# # Save the MI estimate from the calculation and the MI estimates from the simulations
# np.save("../data/MI_estimate_10k_10_times_adt_sct.npy", MI_real)
# np.save("../data/MI_simulations_100_gaussian.npy", MI_simulations)


# In[ ]:


# Next see if MI changes based on the celltypes 

# Extract list of all celltypes from the mdata
celltypes = mdata.obs["celltype.l2"].values

sct_counts = mdata["SCT"].X[:, mdata["SCT"].var.highly_variable].todense()
adt_counts = mdata["ADT"].X

# Sample 500 cells from each celltype (that has at least 500 cells)
celltype_counts_sct = []
celltype_counts_adt = []
celltype_names = []
for celltype in np.unique(celltypes):
    if np.sum(celltypes == celltype) >= 500:
        celltype_names.append(celltype)
        indices = np.random.choice(np.where(celltypes == celltype)[0], 500, replace=False)
        celltype_counts_sct.append(sct_counts[indices, :])
        celltype_counts_adt.append(adt_counts[indices, :])


# In[ ]:


# Calculate the MI for each celltype
MI_celltypes = np.zeros(len(celltype_counts_sct))
for i in range(len(celltype_counts_sct)):
    pmis, _, _ = lmi.estimate(celltype_counts_sct[i], celltype_counts_adt[i])
    MI_celltypes[i] = np.nanmean(pmis) 


# In[ ]:


# Save the MI estimates for each celltype
np.save("../data/MI_celltypes_500.npy", MI_celltypes)
np.save("../data/celltype_names.npy", celltype_names)


# In[ ]:




