import scanpy as sc 

def hvg_getter(data, n_hvgs = 2500, mudata = False, lognormed = False):
    if mudata:
        data = data["SCT"]
    else:
        data = data
        
    data = data.copy() # Prevent overwrites
    
    if lognormed:
        sc.pp.highly_variable_genes(data, n_top_genes = n_hvgs, flavor = "seurat")
    else:
        sc.pp.highly_variable_genes(data, n_top_genes = n_hvgs, flavor = "seurat_v3")
        
    # Convert the boolean to a list of 0 (False) and 1 (True)
    data.var["highly_variable"] = data.var.highly_variable.astype(int)
        
    return data.var.highly_variable